"""Verification of LLM extractions using taln alignment and exact string matching.

Verification is a three-level check for each text extraction:

1. source_rationale: Must be an exact contiguous substring of the manuscript
   (after Unicode normalization + ASCII transliteration via norm_text).
   If exact match fails, the extraction is discarded.

2. group_label: Must align within source_rationale via taln (supports
   non-contiguous token matches, e.g. abbreviations split across tokens).

3. feature_label: Same as group_label.

Text normalization (norm_text from taln) applies:
- Unicode NFC normalization
- ASCII transliteration (unidecode)
- Common punctuation/symbol mapping to ASCII equivalents
- Whitespace/newline collapsing
"""

from taln.taln_aln import align_ng, align_ng_casefold, norm_text, reconstruct_target_by_token


def verify_source_rationale(manuscript_text: str, source_rationale: str) -> dict:
    """Check that source_rationale appears in the manuscript text.

    Uses exact substring match on normalized text. If not found, the
    extraction should be discarded.

    Args:
        manuscript_text: Full manuscript text
        source_rationale: Extracted quote to verify

    Returns:
        dict with keys: found (bool), method (str), normalized_rationale (str)
    """
    norm_ms = norm_text(manuscript_text)
    norm_sr = norm_text(source_rationale)

    if norm_sr in norm_ms:
        return {
            "found": True,
            "method": "exact",
            "normalized_rationale": norm_sr,
        }

    return {
        "found": False,
        "method": "none",
        "normalized_rationale": norm_sr,
    }


def verify_label_in_rationale(source_rationale: str, label: str) -> dict:
    """Check that a label (group_label or feature_label) aligns within source_rationale.

    Uses taln's n-gram alignment which supports non-contiguous matches
    (e.g. "CD4" appearing as part of "CD4+" in the rationale).

    Tries alignment with the label as-is first, then with a leading space prepended.
    Subword tokenizers (tiktoken/BPE) tokenize differently depending on whether
    a token appears at the start of text vs mid-text (e.g. "SPG1" vs " SPG1"),
    so trying both variants improves alignment recall.

    Args:
        source_rationale: The source text to search within
        label: The label to find (group_label or feature_label)

    Returns:
        dict with keys: found (bool), n_alignments (int), reconstructed (str)
    """
    # Case-sensitive alignment (exact token match + space-prefix variant)
    alignments = align_ng(source_rationale, label)
    best_recon = ""
    if alignments:
        best = max(alignments, key=len)
        best_recon = reconstruct_target_by_token(source_rationale, best)

    # Case-insensitive fallback when case-sensitive alignment is missing or
    # incomplete (e.g. LLM extracted "regenerative" but the sentence has
    # "Regenerative" at sentence start — partial match covers fewer tokens)
    if len(best_recon.strip()) < len(label.strip()):
        cf_alignments = align_ng_casefold(source_rationale, label)
        if cf_alignments:
            cf_best = max(cf_alignments, key=len)
            cf_recon = reconstruct_target_by_token(source_rationale, cf_best)
            if len(cf_recon.strip()) > len(best_recon.strip()):
                best_recon = cf_recon

    if best_recon:
        return {
            "found": True,
            "n_alignments": len(alignments or []),
            "reconstructed": best_recon,
        }

    return {
        "found": False,
        "n_alignments": 0,
        "reconstructed": "",
    }


def verify_extraction(manuscript_text: str, record: dict) -> dict:
    """Run all verification checks on a single extraction record.

    Checks:
    1. source_rationale is found in manuscript_text (exact normalized match)
    2. group_label aligns within source_rationale via taln
    3. feature_label aligns within source_rationale via taln

    Args:
        manuscript_text: Full manuscript text
        record: Evidence dict with source_rationale, group_label, feature_label

    Returns:
        dict with verification results for each check
    """
    source_rationale = record.get("source_rationale", "")
    group_label = record.get("group_label", "")
    feature_label = record.get("feature_label", "")

    sr_check = verify_source_rationale(manuscript_text, source_rationale)
    gl_check = verify_label_in_rationale(source_rationale, group_label)
    fl_check = verify_label_in_rationale(source_rationale, feature_label)

    return {
        "source_rationale_found": sr_check["found"],
        "source_rationale_method": sr_check["method"],
        "group_label_found": gl_check["found"],
        "group_label_reconstructed": gl_check["reconstructed"],
        "feature_label_found": fl_check["found"],
        "feature_label_reconstructed": fl_check["reconstructed"],
        "all_verified": sr_check["found"] and gl_check["found"] and fl_check["found"],
    }


def verify_extractions(manuscript_text: str, records: list[dict], verbose: bool = False) -> list[dict]:
    """Verify all extraction records against manuscript text.

    Text records that fail source_rationale verification are discarded.
    Records that pass source_rationale but fail label checks are kept
    with '_verification' metadata. Non-text records (deg, image) are
    passed through unchanged.

    Args:
        manuscript_text: Full manuscript text
        records: List of evidence dicts
        verbose: Print progress

    Returns:
        List of verified records (failed source_rationale extractions removed).
    """
    verified_records = []
    n_text = 0
    n_verified = 0
    n_sr_fail = 0
    n_gl_fail = 0
    n_fl_fail = 0

    for record in records:
        if record.get("source_type") in ("deg", "image"):
            verified_records.append(record)
            continue

        n_text += 1
        result = verify_extraction(manuscript_text, record)

        # Discard extractions where source_rationale not found in manuscript
        if not result["source_rationale_found"]:
            n_sr_fail += 1
            continue

        record["_verification"] = result
        verified_records.append(record)

        if result["all_verified"]:
            n_verified += 1
        else:
            if not result["group_label_found"]:
                n_gl_fail += 1
            if not result["feature_label_found"]:
                n_fl_fail += 1

    if verbose and n_text > 0:
        n_kept = n_text - n_sr_fail
        print(f"   Verification: {n_kept}/{n_text} text extractions kept ({n_sr_fail} discarded)")
        print(f"   Of kept: {n_verified}/{n_kept} fully verified")
        if n_gl_fail:
            print(f"   - {n_gl_fail} group_label not found in source_rationale")
        if n_fl_fail:
            print(f"   - {n_fl_fail} feature_label not found in source_rationale")

    return verified_records
