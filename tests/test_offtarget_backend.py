"""Tests for training-free off-target DuckDB backend."""

from __future__ import annotations

from pathlib import Path

import pytest

from scalpel.models.enums import Genome
from scalpel.offtarget.searcher import OffTargetIndexNotFoundError, OffTargetSearcher

duckdb = pytest.importorskip("duckdb")


def _write_test_index(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with duckdb.connect(str(path)) as conn:
        conn.execute(
            """
            CREATE TABLE pam_sites (
                chromosome VARCHAR,
                position INTEGER,
                strand VARCHAR(1),
                spacer VARCHAR(23),
                pam VARCHAR(3)
            )
            """
        )
        conn.executemany(
            "INSERT INTO pam_sites VALUES (?, ?, ?, ?, ?)",
            [
                ("chr1", 100, "+", "AAAAAAAAAAAAAAAAAAAA", "AGG"),  # exact on-target
                ("chr1", 200, "+", "CAAAAAAAAAAAAAAAAAAA", "AGG"),  # 1 mismatch
                ("chr1", 300, "+", "CCCCCCCCCCCCCCCCCCCC", "AGG"),  # many mismatches
            ],
        )


def test_search_database_returns_filtered_hits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db_path = tmp_path / "test.duckdb"
    _write_test_index(db_path)
    monkeypatch.setenv("SCALPEL_OFFTARGET_DB", str(db_path))

    searcher = OffTargetSearcher(Genome.HUMAN_GRCH38, max_mismatches=2)
    analysis = searcher.search("AAAAAAAAAAAAAAAAAAAA")

    assert analysis.total_sites >= 1
    assert all(site.mismatch_count > 0 for site in analysis.sites)
    assert any(site.mismatch_count == 1 for site in analysis.sites)
    assert all(site.cutting_probability > 0.01 for site in analysis.sites)


def test_missing_index_raises_actionable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    searcher = OffTargetSearcher(Genome.HUMAN_GRCH38, max_mismatches=2)
    monkeypatch.delenv("SCALPEL_OFFTARGET_DB", raising=False)
    monkeypatch.setattr(searcher, "_resolve_database_path", lambda: None)

    with pytest.raises(OffTargetIndexNotFoundError, match="build_index"):
        searcher.search("AAAAAAAAAAAAAAAAAAAA")


def test_resolve_database_path_uses_env_var_first(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    env_db = tmp_path / "explicit.duckdb"
    env_db.write_bytes(b"")
    monkeypatch.setenv("SCALPEL_OFFTARGET_DB", str(env_db))

    # Even if a default path could exist, env var should win.
    monkeypatch.setattr("scalpel.offtarget.searcher.Path.home", lambda *_args, **_kwargs: tmp_path)
    default_db = tmp_path / ".scalpel" / "offtargets" / "GRCh38.duckdb"
    default_db.parent.mkdir(parents=True, exist_ok=True)
    default_db.write_bytes(b"")

    searcher = OffTargetSearcher(Genome.HUMAN_GRCH38)
    resolved = searcher._resolve_database_path()

    assert resolved == env_db
