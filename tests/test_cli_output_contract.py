"""CLI contract tests for machine-readable output and error messaging."""

from __future__ import annotations

import json
from types import SimpleNamespace

from typer.testing import CliRunner

from scalpel.cli import app
from scalpel.models.data_classes import EfficiencyScore, SpacerCandidate
from scalpel.models.enums import Strand
from scalpel.offtarget.searcher import OffTargetIndexNotFoundError


def _patch_design_pipeline(monkeypatch) -> None:
    spacer = SpacerCandidate(
        spacer_sequence="ATCGATCGATCGATCGATCG",
        pam_sequence="AGG",
        strand=Strand.PLUS,
        genomic_start=100,
        genomic_end=120,
        cut_site=117,
        context_sequence="NNNNNATCGATCGATCGATCGATCGAGGNNNNN",
    )

    class _FakeResolver:
        def __init__(self, *_args, **_kwargs):
            pass

        def resolve(self, _spec):
            return SimpleNamespace(
                sequence="A" * 200,
                chromosome="chr1",
                start=100,
                end=300,
                strand=Strand.PLUS,
                gene_info=SimpleNamespace(symbol="TP53", gene_id="ENSG00000141510"),
                transcript=None,
            )

    class _FakeExtractor:
        def __init__(self, *_args, **_kwargs):
            pass

        def extract_spacers(self, *_args, **_kwargs):
            return [spacer]

    class _FakeScoredGuide:
        def __init__(self, candidate):
            self.spacer = candidate
            self.efficiency = EfficiencyScore(
                overall_score=0.77,
                interpretation="Good efficiency",
            )

    class _FakeScorer:
        def score_batch(self, spacers, _modality, n_top=None):
            limit = n_top if n_top is not None else len(spacers)
            return [_FakeScoredGuide(s) for s in spacers[:limit]]

    monkeypatch.setattr("scalpel.genome.target_resolver.TargetResolver", _FakeResolver)
    monkeypatch.setattr("scalpel.design.SpacerExtractor", _FakeExtractor)
    monkeypatch.setattr("scalpel.design.efficiency.EnsembleScorer", _FakeScorer)
    monkeypatch.setattr("scalpel.core.red_flags.detect_red_flags", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        "scalpel.core.red_flags.summarize_red_flags",
        lambda _flags: {"total_flags": 0, "interpretation": "No critical red flags"},
    )


def test_design_json_output_is_parseable_in_noninteractive_mode(monkeypatch) -> None:
    _patch_design_pipeline(monkeypatch)
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "design",
            "--sequence",
            "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCG",
            "--format",
            "json",
            "--n-guides",
            "1",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["status"] == "success"
    assert len(payload["guides"]) == 1
    assert "Resolving target" not in result.output


def test_design_tsv_output_has_no_progress_noise(monkeypatch) -> None:
    _patch_design_pipeline(monkeypatch)
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "design",
            "--sequence",
            "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCG",
            "--format",
            "tsv",
            "--n-guides",
            "1",
        ],
    )

    assert result.exit_code == 0
    assert result.stdout.splitlines()[0] == "rank\tspacer\tpam\tscore"
    assert "Resolving target" not in result.output


def test_offtarget_missing_index_message_is_actionable(monkeypatch) -> None:
    class _MissingIndexSearcher:
        def __init__(self, *_args, **_kwargs):
            pass

        def search(self, _spacer):
            raise OffTargetIndexNotFoundError(
                "Off-target index not found for GRCh38. "
                "Build it with: python -m scalpel.offtarget.build_index --genome GRCh38"
            )

    class _FakeRiskCalculator:
        def __init__(self, *_args, **_kwargs):
            pass

    monkeypatch.setattr("scalpel.offtarget.OffTargetSearcher", _MissingIndexSearcher)
    monkeypatch.setattr("scalpel.offtarget.RiskCalculator", _FakeRiskCalculator)

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "offtarget",
            "--spacer",
            "ATCGATCGATCGATCGATCG",
            "--genome",
            "GRCh38",
        ],
    )

    assert result.exit_code == 1
    assert "build_index" in result.output
