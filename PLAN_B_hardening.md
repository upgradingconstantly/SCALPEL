# PLAN B: Reliability Hardening Roadmap (Deferred)

This file captures the deferred hardening plan and is intentionally **not** implemented in Plan A.

## Goals
- Push runtime reliability toward 10/10.
- Push architecture quality toward 10/10.
- Expand test rigor beyond the Plan A target (~6/10).
- Keep optional ML support training-free by default.

## B1. Off-target Performance and Completeness
1. Add chunked candidate scanning and adaptive candidate caps for very large genome indexes.
2. Add backend metadata in responses:
   - `backend`
   - `db_path`
   - `candidate_count`
   - `search_ms`
3. Add explicit warning semantics for partial search conditions (cap-hit, timeout, truncated candidate set).

## B2. API/Dependency Hardening
1. Pin/validate FastAPI + Starlette + httpx compatibility for robust API integration tests.
2. Add stable API integration tests against ASGI transport.
3. Add regression tests for design/off-target/plan contracts.

## B3. Data Tooling Completion
1. Ship a supported gene database build path (or remove unsupported references entirely).
2. Standardize path conventions for genome/off-target assets across:
   - `scalpel/config.py`
   - `scalpel/offtarget/build_index.py`
   - docs

## B4. Optional Pretrained Model Path (No Local Training)
1. Add optional pretrained CRISPRon weight download flow.
2. Validate checksum/version at load time.
3. Keep rule-based fallback as first-class behavior when weights are unavailable.

## B5. Testing Expansion
1. Add larger off-target fixtures and mismatch edge-case suites.
2. Add modality-specific property tests for ranking invariants and bounds.
3. Add CI matrix for Python versions + dependency pins.

## B6. Acceptance Criteria
1. Off-target backend returns deterministic, benchmarked results with metadata.
2. API integration tests pass reliably with pinned dependencies.
3. Data setup docs and commands are internally consistent and executable.
4. Optional pretrained model path is documented, reproducible, and never mandatory.

## Status
- Deferred by request.
- Do not execute within Plan A sprint.
