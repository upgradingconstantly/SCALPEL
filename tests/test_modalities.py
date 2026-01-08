"""
Tests for modality-specific designers.
"""

import pytest
from scalpel.models.enums import EditModality
from scalpel.models.data_classes import TargetSpecification, GeneInfo


class TestKnockoutDesigner:
    """Tests for CRISPR knockout design."""
    
    def test_skips_first_exon(self, sample_gene_info: GeneInfo):
        """First exon should be skipped (alt start sites)."""
        pass
    
    def test_skips_last_exon(self, sample_gene_info: GeneInfo):
        """Last exon should be skipped (NMD escape)."""
        pass
    
    def test_prefers_early_cds_position(self):
        """Guides earlier in CDS should score higher."""
        pass
    
    def test_prefers_constitutive_exons(self):
        """Constitutive exons should be preferred."""
        pass
    
    def test_domain_targeting_bonus(self):
        """Guides targeting functional domains should get bonus."""
        pass


class TestCRISPRiDesigner:
    """Tests for CRISPRi design."""
    
    def test_targets_near_tss(self):
        """Guides should target near TSS."""
        pass
    
    def test_optimal_tss_distance(self):
        """Optimal distance is TSS to TSS+50bp."""
        pass


class TestCRISPRaDesigner:
    """Tests for CRISPRa design."""
    
    def test_targets_promoter(self):
        """Guides should target promoter region."""
        pass
    
    def test_optimal_upstream_distance(self):
        """Optimal distance is 50-400bp upstream of TSS."""
        pass


class TestBaseEditorDesigner:
    """Tests for base editor design."""
    
    def test_cbe_targets_c(self):
        """CBE should target C for C→T conversion."""
        pass
    
    def test_abe_targets_a(self):
        """ABE should target A for A→G conversion."""
        pass
    
    def test_editing_window_position(self):
        """Edit should fall in editing window (pos 4-8 for CBE)."""
        pass
    
    def test_bystander_count(self):
        """Bystander edits should be counted and penalized."""
        pass


class TestPrimeEditorDesigner:
    """Tests for prime editor design."""
    
    def test_pbs_length_range(self):
        """PBS should be 8-15 nt."""
        pass
    
    def test_pbs_tm_calculation(self):
        """PBS Tm should be calculated correctly."""
        pass
    
    def test_rtt_encodes_edit(self):
        """RTT should encode the desired edit."""
        pass
    
    def test_nicking_guide_distance(self):
        """PE3 nicking guide should be 40-120bp away."""
        pass


class TestPluginRegistry:
    """Tests for modality plugin registration."""
    
    def test_all_modalities_registered(self):
        """All modalities should be registered."""
        # Import modalities to trigger registration decorators
        import scalpel.modalities  # noqa: F401
        from scalpel.core.plugin_registry import ModalityRegistry
        
        registered = ModalityRegistry.list_modalities()
        assert EditModality.KNOCKOUT in registered
        assert EditModality.INTERFERENCE in registered
        assert EditModality.ACTIVATION in registered
        assert EditModality.BASE_EDIT_CBE in registered
        assert EditModality.BASE_EDIT_ABE in registered
        assert EditModality.PRIME_EDIT in registered
    
    def test_get_plugin_instance(self):
        """Should be able to get plugin instances."""
        # Import modalities to trigger registration decorators
        import scalpel.modalities  # noqa: F401
        from scalpel.core.plugin_registry import ModalityRegistry
        
        ko_designer = ModalityRegistry.get_instance(EditModality.KNOCKOUT)
        assert ko_designer is not None
        assert ko_designer.display_name == "CRISPR Knockout"
