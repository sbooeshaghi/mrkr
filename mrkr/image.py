"""Image processing utilities for mrkr."""

import base64
import io
from pathlib import Path
from typing import Optional, Tuple

from PIL import Image


def resize_and_encode_image(
    image_path: Path,
    max_size: int = 1568,
    quality: int = 85,
) -> Optional[Tuple[str, str]]:
    """Resize an image and encode as base64 for Claude API.

    Claude API has image size limits. This function:
    - Resizes to fit within max_size x max_size (maintains aspect ratio)
    - Converts to RGB (handles RGBA, grayscale, etc.)
    - Saves as JPEG and encodes as base64
    - Uses 1568px default (safe for Claude with multiple images)

    Args:
        image_path: Path to input image (any format PIL supports: TIF, PNG, JPG, etc.)
        max_size: Maximum dimension in pixels (default: 1568)
        quality: JPEG quality 1-100 (default: 85)

    Returns:
        Tuple of (base64_data, media_type) or None if processing fails

    Example:
        >>> from pathlib import Path
        >>> data, media_type = resize_and_encode_image(Path("figure.tif"))
        >>> print(media_type)
        'image/jpeg'
    """
    try:
        # Open image
        image = Image.open(image_path)

        # Convert to RGB if necessary (handles RGBA, P, L, etc.)
        if image.mode not in ('RGB', 'L'):
            image = image.convert('RGB')
        elif image.mode == 'L':
            # Convert grayscale to RGB for consistency
            image = image.convert('RGB')

        # Calculate new size maintaining aspect ratio
        width, height = image.size
        if width > max_size or height > max_size:
            # Resize to fit within max_size x max_size
            if width > height:
                new_width = max_size
                new_height = int(height * (max_size / width))
            else:
                new_height = max_size
                new_width = int(width * (max_size / height))

            image = image.resize((new_width, new_height), Image.Resampling.LANCZOS)

        # Convert to JPEG and encode as base64
        buffer = io.BytesIO()
        image.save(buffer, format='JPEG', quality=quality, optimize=True)
        buffer.seek(0)

        base64_data = base64.b64encode(buffer.read()).decode('utf-8')

        return base64_data, 'image/jpeg'

    except Exception as e:
        # Return None on error - caller will handle
        return None


def process_image_for_api(
    image_path: Path,
    max_size: int = 1568,
    verbose: bool = False,
) -> Optional[Tuple[str, str]]:
    """Process a single image file for API consumption.

    Wraps resize_and_encode_image with error handling and optional logging.

    Args:
        image_path: Path to image file
        max_size: Maximum dimension in pixels
        verbose: Print processing information

    Returns:
        Tuple of (base64_data, media_type) or None if processing fails
    """
    if verbose:
        print(f"  📷 Processing image: {image_path.name}")

    result = resize_and_encode_image(image_path, max_size=max_size)

    if result is None:
        if verbose:
            print(f"  ⚠️  Failed to process: {image_path.name}")
        return None

    if verbose:
        base64_data, media_type = result
        size_kb = len(base64_data) / 1024
        print(f"     → Resized to ≤{max_size}px, {size_kb:.1f} KB")

    return result
