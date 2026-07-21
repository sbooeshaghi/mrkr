"""Tests for the production mrkr command surface."""

from click.testing import CliRunner

from mrkr.cli import cli


def test_help_exposes_only_the_core_workflow():
    result = CliRunner().invoke(cli, ["--help"])

    assert result.exit_code == 0
    assert "extract" in result.output
    assert "ground" in result.output
    assert "validate" in result.output
    assert "\n  verify " not in result.output
    assert "\n  generate " not in result.output


def test_ground_requires_the_source_manuscript(tmp_path):
    claims = tmp_path / "claims.json"
    claims.write_text("{}", encoding="utf-8")

    result = CliRunner().invoke(
        cli, ["ground", str(claims), "--output", str(tmp_path / "onto.json")]
    )

    assert result.exit_code == 2
    assert "Missing option '--manuscript'" in result.output


def test_ground_requires_organism(tmp_path):
    manuscript = tmp_path / "paper.txt"
    manuscript.write_text("paper", encoding="utf-8")
    claims = tmp_path / "claims.json"
    claims.write_text("{}", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        [
            "ground",
            str(claims),
            "--manuscript",
            str(manuscript),
            "--output",
            str(tmp_path / "onto.json"),
        ],
    )

    assert result.exit_code == 2
    assert "Missing option '--organism'" in result.output


def test_validate_reports_malformed_json_without_a_traceback(tmp_path):
    manuscript = tmp_path / "paper.txt"
    manuscript.write_text("paper", encoding="utf-8")
    malformed = tmp_path / "broken.json"
    malformed.write_text("{", encoding="utf-8")

    result = CliRunner().invoke(
        cli,
        ["validate", str(malformed), "--manuscript", str(manuscript)],
    )

    assert result.exit_code == 1
    assert "Error:" in result.output
    assert "Traceback" not in result.output
