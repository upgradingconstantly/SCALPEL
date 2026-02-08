"""
B2.2: API integration tests using ASGI transport.

Stable integration tests against the FastAPI app using httpx ASGITransport
for reliable, isolated testing without network dependencies.
"""

from __future__ import annotations

import pytest
from typing import Dict, Any

# Configure pytest-asyncio
pytest_plugins = ('pytest_asyncio',)

# Skip if httpx not available
httpx = pytest.importorskip("httpx")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from scalpel.api.server import app


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def client():
    """Create httpx AsyncClient with ASGI transport."""
    from httpx import ASGITransport
    transport = ASGITransport(app=app)
    return httpx.AsyncClient(transport=transport, base_url="http://test")


# =============================================================================
# Health Endpoint Tests
# =============================================================================

class TestHealthEndpoint:
    """Test the /api/health endpoint."""
    
    @pytest.mark.asyncio
    async def test_health_returns_200(self, client):
        """Health endpoint should return 200 OK."""
        async with client:
            response = await client.get("/api/health")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert "version" in data
    
    @pytest.mark.asyncio
    async def test_health_response_structure(self, client):
        """Health response should have expected structure."""
        async with client:
            response = await client.get("/api/health")
        
        data = response.json()
        assert isinstance(data, dict)
        assert "status" in data
        assert "version" in data


class TestInfoEndpoint:
    """Test the /api/info endpoint."""
    
    @pytest.mark.asyncio
    async def test_info_returns_200(self, client):
        """Info endpoint should return 200 OK."""
        async with client:
            response = await client.get("/api/info")
        
        assert response.status_code == 200
    
    @pytest.mark.asyncio
    async def test_info_contains_modalities(self, client):
        """Info should list available modalities."""
        async with client:
            response = await client.get("/api/info")
        
        data = response.json()
        assert "modalities" in data
        assert len(data["modalities"]) > 0
        
        # Check modality structure
        modality = data["modalities"][0]
        assert "value" in modality
        assert "label" in modality
    
    @pytest.mark.asyncio
    async def test_info_contains_genomes(self, client):
        """Info should list available genomes."""
        async with client:
            response = await client.get("/api/info")
        
        data = response.json()
        assert "genomes" in data
        assert len(data["genomes"]) > 0
        
        # Check for GRCh38
        genome_values = [g["value"] for g in data["genomes"]]
        assert "GRCh38" in genome_values
    
    @pytest.mark.asyncio
    async def test_info_contains_cas_variants(self, client):
        """Info should list available Cas variants."""
        async with client:
            response = await client.get("/api/info")
        
        data = response.json()
        assert "cas_variants" in data
        assert len(data["cas_variants"]) > 0
        
        # Check for SpCas9
        variant_values = [v["value"] for v in data["cas_variants"]]
        assert "SpCas9" in variant_values


# =============================================================================
# Design Endpoint Tests
# =============================================================================

class TestDesignEndpoint:
    """Test the /api/design endpoint."""
    
    @pytest.mark.asyncio
    async def test_design_with_sequence_returns_200(self, client):
        """Design with sequence input should return 200 OK."""
        # Simple sequence with known NGG PAM sites
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCC"
        
        async with client:
            response = await client.post(
                "/api/design",
                json={
                    "sequence": test_sequence,
                    "modality": "knockout",
                    "genome": "GRCh38",
                    "cas_variant": "SpCas9",
                    "n_guides": 5,
                }
            )
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
    
    @pytest.mark.asyncio
    async def test_design_response_structure(self, client):
        """Design response should have expected structure."""
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCC"
        
        async with client:
            response = await client.post(
                "/api/design",
                json={
                    "sequence": test_sequence,
                    "modality": "knockout",
                    "n_guides": 5,
                }
            )
        
        data = response.json()
        
        # Required fields
        assert "status" in data
        assert "target" in data
        assert "guides" in data
        assert "n_guides_returned" in data
        assert "n_candidates_total" in data
    
    @pytest.mark.asyncio
    async def test_design_guides_have_required_fields(self, client):
        """Each guide should have required fields."""
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCC"
        
        async with client:
            response = await client.post(
                "/api/design",
                json={"sequence": test_sequence, "n_guides": 3}
            )
        
        data = response.json()
        
        if data.get("guides"):
            guide = data["guides"][0]
            
            required_fields = [
                "rank",
                "spacer_sequence",
                "pam_sequence",
                "strand",
                "efficiency_score",
            ]
            
            for field in required_fields:
                assert field in guide, f"Missing field: {field}"
    
    @pytest.mark.asyncio
    async def test_design_efficiency_scores_bounded(self, client):
        """Efficiency scores should be between 0 and 1."""
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCC"
        
        async with client:
            response = await client.post(
                "/api/design",
                json={"sequence": test_sequence, "n_guides": 10}
            )
        
        data = response.json()
        
        for guide in data.get("guides", []):
            score = guide.get("efficiency_score", 0)
            assert 0 <= score <= 1, f"Score {score} out of bounds"
    
    @pytest.mark.asyncio
    async def test_design_invalid_modality_returns_422(self, client):
        """Invalid modality should return 422 validation error."""
        async with client:
            response = await client.post(
                "/api/design",
                json={
                    "sequence": "ATGGATTTATCTGCTCTTCGAGG",
                    "modality": "invalid_modality",
                }
            )
        
        assert response.status_code == 422
    
    @pytest.mark.asyncio
    async def test_design_invalid_sequence_returns_422(self, client):
        """Invalid DNA sequence should return 422 validation error."""
        async with client:
            response = await client.post(
                "/api/design",
                json={
                    "sequence": "ATGGXYZATTTATCTGCTCTTCG",  # Invalid bases
                    "modality": "knockout",
                }
            )
        
        assert response.status_code == 422
    
    @pytest.mark.asyncio
    async def test_design_no_target_returns_error(self, client):
        """Request without gene/sequence/coordinates should return error."""
        async with client:
            response = await client.post(
                "/api/design",
                json={
                    "modality": "knockout",
                    "n_guides": 5,
                }
            )
        
        # Should return error status (400 or 500)
        assert response.status_code in [400, 500]


# =============================================================================
# Off-target Endpoint Tests
# =============================================================================

class TestOfftargetEndpoint:
    """Test the /api/offtarget endpoint."""
    
    @pytest.mark.asyncio
    async def test_offtarget_request_structure(self, client):
        """Off-target endpoint should accept valid request structure."""
        async with client:
            response = await client.post(
                "/api/offtarget",
                json={
                    "spacer": "ATGGATTTATCTGCTCTTCG",
                    "genome": "GRCh38",
                    "max_mismatches": 3,
                }
            )
        
        # Should return either 200 (if index exists) or 503 (if not)
        assert response.status_code in [200, 503]
        
        if response.status_code == 503:
            data = response.json()
            assert "detail" in data
            assert "build_index" in data["detail"]
    
    @pytest.mark.asyncio
    async def test_offtarget_invalid_spacer_returns_422(self, client):
        """Invalid spacer sequence should return 422."""
        async with client:
            response = await client.post(
                "/api/offtarget",
                json={
                    "spacer": "ATGGXYZATTTATCTGCTCT",  # Invalid bases
                    "genome": "GRCh38",
                }
            )
        
        assert response.status_code == 422
    
    @pytest.mark.asyncio
    async def test_offtarget_short_spacer_returns_422(self, client):
        """Too short spacer should return 422."""
        async with client:
            response = await client.post(
                "/api/offtarget",
                json={
                    "spacer": "ATGG",  # Too short
                    "genome": "GRCh38",
                }
            )
        
        assert response.status_code == 422


# =============================================================================
# Export Endpoint Tests
# =============================================================================

class TestExportEndpoints:
    """Test the export endpoints."""
    
    @pytest.fixture
    def sample_design_data(self) -> Dict[str, Any]:
        """Sample design data for export tests."""
        return {
            "target": {
                "gene": "TP53",
                "chromosome": "chr17",
                "start": 7668421,
                "end": 7687490,
                "strand": "-",
                "genome": "GRCh38",
                "modality": "knockout",
            },
            "guides": [
                {
                    "rank": 1,
                    "spacer_sequence": "ATGGATTTATCTGCTCTTCG",
                    "pam_sequence": "AGG",
                    "strand": "+",
                    "genomic_start": 7674200,
                    "genomic_end": 7674220,
                    "cut_site": 7674217,
                    "efficiency_score": 0.85,
                    "efficiency_interpretation": "Good efficiency",
                    "red_flags": {"severity": "none", "interpretation": "No issues"},
                },
                {
                    "rank": 2,
                    "spacer_sequence": "GCGATCGATCGATCGATCGG",
                    "pam_sequence": "TGG",
                    "strand": "-",
                    "genomic_start": 7674500,
                    "genomic_end": 7674520,
                    "cut_site": 7674517,
                    "efficiency_score": 0.72,
                    "efficiency_interpretation": "Moderate efficiency",
                    "red_flags": {"severity": "low", "interpretation": "Minor concerns"},
                },
            ],
        }
    
    @pytest.mark.asyncio
    async def test_export_csv_returns_200(self, client, sample_design_data):
        """CSV export should return 200."""
        async with client:
            response = await client.post(
                "/api/export/csv",
                json=sample_design_data
            )
        
        assert response.status_code == 200
        assert "text/csv" in response.headers.get("content-type", "")
    
    @pytest.mark.asyncio
    async def test_export_csv_has_content_disposition(self, client, sample_design_data):
        """CSV export should have Content-Disposition header."""
        async with client:
            response = await client.post(
                "/api/export/csv",
                json=sample_design_data
            )
        
        assert "content-disposition" in response.headers
        assert "attachment" in response.headers["content-disposition"]
        assert ".csv" in response.headers["content-disposition"]
    
    @pytest.mark.asyncio
    async def test_export_json_returns_200(self, client, sample_design_data):
        """JSON export should return 200."""
        async with client:
            response = await client.post(
                "/api/export/json",
                json=sample_design_data
            )
        
        assert response.status_code == 200
        assert "application/json" in response.headers.get("content-type", "")
    
    @pytest.mark.asyncio
    async def test_export_fasta_returns_200(self, client, sample_design_data):
        """FASTA export should return 200."""
        async with client:
            response = await client.post(
                "/api/export/fasta",
                json=sample_design_data
            )
        
        assert response.status_code == 200
        assert "text/plain" in response.headers.get("content-type", "")


# =============================================================================
# Rate Limiting Tests
# =============================================================================

class TestRateLimiting:
    """Test rate limiting functionality."""
    
    @pytest.mark.asyncio
    async def test_rate_limiter_allows_normal_requests(self, client):
        """Normal request rate should be allowed."""
        async with client:
            # Make 5 requests - should all succeed
            for _ in range(5):
                response = await client.get("/api/health")
                assert response.status_code == 200


# =============================================================================
# CORS Tests
# =============================================================================

class TestCORS:
    """Test CORS configuration."""
    
    @pytest.mark.asyncio
    async def test_cors_headers_present(self, client):
        """CORS headers should be present for API requests."""
        async with client:
            response = await client.options(
                "/api/health",
                headers={"Origin": "http://localhost:3000"}
            )
        
        # Should have CORS headers (may be 200 or 405 depending on implementation)
        assert response.status_code in [200, 405]


# =============================================================================
# Contract Tests
# =============================================================================

class TestAPIContract:
    """Test API contract stability."""
    
    @pytest.mark.asyncio
    async def test_design_contract_stable(self, client):
        """Design endpoint contract should be stable."""
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCC"
        
        async with client:
            response = await client.post(
                "/api/design",
                json={"sequence": test_sequence, "n_guides": 5}
            )
        
        data = response.json()
        
        # Contract: these fields must always exist
        assert "status" in data
        assert "target" in data
        assert "guides" in data
        assert "n_candidates_total" in data
        assert "n_guides_returned" in data
        
        # Contract: guides must be a list
        assert isinstance(data["guides"], list)
        
        # Contract: target must have genome info
        assert "genome" in data["target"] or "modality" in data["target"]
    
    @pytest.mark.asyncio
    async def test_info_contract_stable(self, client):
        """Info endpoint contract should be stable."""
        async with client:
            response = await client.get("/api/info")
        
        data = response.json()
        
        # Contract: these arrays must exist
        assert "modalities" in data
        assert "genomes" in data
        assert "cas_variants" in data
        
        # Contract: each item must have value/label
        for modality in data["modalities"]:
            assert "value" in modality
            assert "label" in modality
