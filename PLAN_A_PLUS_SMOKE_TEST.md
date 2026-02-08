# Plan A + Smoke Test

Complete Plan A sprint by committing core files and validating API endpoints.

---

## Part 1: Plan A Commit

### Files to Stage

**Core Modules**
- `scalpel/design/efficiency/ensemble.py` - Efficiency scoring
- `scalpel/cli.py` - CLI
- `scalpel/api/server.py` - FastAPI server

**Off-target Module**
- `scalpel/offtarget/__init__.py`
- `scalpel/offtarget/searcher.py`
- `scalpel/offtarget/build_index.py`
- `scalpel/offtarget/cfd_scorer.py`

**Tests**
- `tests/test_efficiency.py`
- `tests/test_modalities.py`
- `tests/test_spacer.py`
- `tests/test_cli_output_contract.py`
- `tests/test_ensemble_regression.py`
- `tests/test_offtarget_backend.py`

**Docs**
- `data/GENOMIC_DATA_DOWNLOAD.md`
- `PLAN_B_hardening.md`

### Commands

```bash
git add scalpel/design/efficiency/ensemble.py scalpel/cli.py scalpel/api/server.py scalpel/offtarget/ tests/ data/GENOMIC_DATA_DOWNLOAD.md PLAN_B_hardening.md

git commit -m "Plan A: Core efficiency, CLI, API, off-target, and tests"

git push origin main
```

---

## Part 2: API Smoke Test

### Start Server
```bash
python -m scalpel.api.server
# Runs at http://127.0.0.1:8000
```

### Test Sequence
```
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT
```

### 1. `/api/design`
```bash
curl -X POST http://127.0.0.1:8000/api/design -H "Content-Type: application/json" -d "{\"sequence\": \"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT\", \"cas_type\": \"SpCas9\", \"n_guides\": 5}"
```
**Expect:** `200 OK`, JSON with `guides` array containing `spacer`, `pam`, `strand`, `position`, `efficiency_score`

### 2. `/api/design/stream`
```bash
curl -X POST http://127.0.0.1:8000/api/design/stream -H "Content-Type: application/json" -d "{\"sequence\": \"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT\", \"cas_type\": \"SpCas9\", \"n_guides\": 5}"
```
**Expect:** `200 OK`, SSE stream with progress and final guides

### 3. `/api/offtarget`
```bash
curl -X POST http://127.0.0.1:8000/api/offtarget -H "Content-Type: application/json" -d "{\"spacer\": \"ATGGATTTATCTGCTCTTCG\", \"genome\": \"GRCh38\", \"max_mismatches\": 3}"
```
**Expect:** `200 OK`, JSON with `off_targets`, `n_candidates_total`, `n_guides_returned`

### 4. 503 Test (No Index)
```bash
$env:SCALPEL_OFFTARGET_DB = "C:\nonexistent\path.duckdb"
python -m scalpel.api.server
# Then hit /api/offtarget
```
**Expect:** `503 Service Unavailable`

---

## Checklist

| Task | Status |
|------|--------|
| Stage Plan A files | ⬜ |
| Commit | ⬜ |
| `/api/design` works | ⬜ |
| `/api/design/stream` works | ⬜ |
| `/api/offtarget` works | ⬜ |
| 503 when no index | ⬜ |
