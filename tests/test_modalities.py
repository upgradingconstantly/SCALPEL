"""Tests for modality-specific designers."""

from __future__ import annotations

from scalpel.core.plugin_registry import ModalityRegistry
from scalpel.modalities.activation import ActivationDesigner
from scalpel.modalities.base_editor import ABEDesigner, CBEDesigner
from scalpel.modalities.interference import InterferenceDesigner
from scalpel.modalities.knockout import KnockoutDesigner
from scalpel.modalities.prime_editor import PrimeEditorDesigner
from scalpel.models.data_classes import (
    Exon,
    GeneInfo,
    ResolvedTarget,
    SpacerCandidate,
    TargetSpecification,
    Transcript,
)
from scalpel.models.enums import EditModality, Strand


class _TestInterferenceDesigner(InterferenceDesigner):
    def score_position(self, guide, target) -> float:  # pragma: no cover - shim for ABC
        return 0.5


class _TestActivationDesigner(ActivationDesigner):
    def score_position(self, guide, target) -> float:  # pragma: no cover - shim for ABC
        return 0.5


def _spacer(cut_site: int, strand: Strand = Strand.PLUS, seq: str = "ATCGATCGATCGATCGATCG") -> SpacerCandidate:
    return SpacerCandidate(
        spacer_sequence=seq,
        pam_sequence="AGG",
        strand=strand,
        genomic_start=cut_site - 20,
        genomic_end=cut_site,
        cut_site=cut_site,
        context_sequence=f"NNNNN{seq}AGGNNNNN",
    )


def _target(
    spec: TargetSpecification,
    strand: Strand = Strand.PLUS,
    tss: int = 1000,
    include_exons: bool = True,
) -> ResolvedTarget:
    if strand == Strand.PLUS:
        transcript_start = tss
        transcript_end = tss + 1200
    else:
        transcript_start = tss - 1200
        transcript_end = tss

    exons = []
    if include_exons:
        exons = [
            Exon(exon_number=1, start=200, end=260, is_constitutive=True),
            Exon(exon_number=2, start=320, end=420, is_constitutive=True),
            Exon(exon_number=3, start=500, end=560, is_constitutive=True),
        ]

    transcript = Transcript(
        transcript_id="ENST_TEST_001",
        gene_id="ENSG_TEST_001",
        chromosome="chr1",
        start=transcript_start,
        end=transcript_end,
        strand=strand,
        cds_start=250,
        cds_end=550,
        exons=exons,
    )

    gene_info = GeneInfo(
        gene_id="ENSG_TEST_001",
        symbol="TEST1",
        chromosome="chr1",
        start=200,
        end=600,
        strand=strand,
        transcripts=[transcript],
    )

    return ResolvedTarget(
        specification=spec,
        chromosome="chr1",
        start=200,
        end=600,
        strand=strand,
        sequence="A" * 1200,
        gene_info=gene_info,
        transcript=transcript,
    )


class TestKnockoutDesigner:
    """Tests for CRISPR knockout design."""

    def test_skips_first_exon(self, sample_target_spec: TargetSpecification) -> None:
        designer = KnockoutDesigner()
        target = _target(sample_target_spec, include_exons=True)
        spacers = [_spacer(230), _spacer(360), _spacer(530)]
        kept = designer._filter_by_exon(spacers, target, designer.get_default_parameters())
        assert all(s.cut_site != 230 for s in kept)

    def test_skips_last_exon(self, sample_target_spec: TargetSpecification) -> None:
        designer = KnockoutDesigner()
        target = _target(sample_target_spec, include_exons=True)
        spacers = [_spacer(230), _spacer(360), _spacer(530)]
        kept = designer._filter_by_exon(spacers, target, designer.get_default_parameters())
        assert all(s.cut_site != 530 for s in kept)

    def test_prefers_early_cds_position(self, sample_target_spec: TargetSpecification) -> None:
        designer = KnockoutDesigner()
        target = _target(sample_target_spec, include_exons=True)
        early = designer._score_position(_spacer(260), target)
        late = designer._score_position(_spacer(540), target)
        assert early > late

    def test_prefers_constitutive_exons(self) -> None:
        params = KnockoutDesigner().get_default_parameters()
        assert params["prefer_constitutive"] is True

    def test_domain_targeting_bonus(self, sample_target_spec: TargetSpecification) -> None:
        designer = KnockoutDesigner()
        target = _target(sample_target_spec, include_exons=True)
        target.gene_info.gene_id = "ENSG00000141510"  # TP53 demo domains exist
        boosted = designer._score_domain(_spacer(7674200), target)
        neutral = designer._score_domain(_spacer(100), target)
        assert boosted >= neutral


class TestCRISPRiDesigner:
    """Tests for CRISPRi design."""

    def test_targets_near_tss(self, sample_target_spec: TargetSpecification) -> None:
        designer = _TestInterferenceDesigner()
        target = _target(sample_target_spec, strand=Strand.PLUS, tss=1000, include_exons=False)
        tss = designer._get_tss(target)
        spacers = [_spacer(980), _spacer(1100), _spacer(1500)]
        kept = designer._filter_by_tss_window(spacers, tss, target)
        assert all(s in kept for s in spacers[:2])
        assert spacers[2] not in kept

    def test_optimal_tss_distance(self, sample_target_spec: TargetSpecification) -> None:
        designer = _TestInterferenceDesigner()
        target = _target(sample_target_spec, strand=Strand.PLUS, tss=1000, include_exons=False)
        tss = designer._get_tss(target)
        optimal = designer._score_tss_distance(_spacer(tss + 100), tss)
        distal = designer._score_tss_distance(_spacer(tss + 280), tss)
        assert optimal > distal


class TestCRISPRaDesigner:
    """Tests for CRISPRa design."""

    def test_targets_promoter(self, sample_target_spec: TargetSpecification) -> None:
        designer = _TestActivationDesigner()
        target = _target(sample_target_spec, strand=Strand.PLUS, tss=1000, include_exons=False)
        tss = designer._get_tss(target)
        spacers = [_spacer(700), _spacer(900), _spacer(1150)]
        kept = designer._filter_by_promoter(spacers, tss, target)
        assert spacers[0] in kept
        assert spacers[1] in kept
        assert spacers[2] not in kept

    def test_optimal_upstream_distance(self, sample_target_spec: TargetSpecification) -> None:
        designer = _TestActivationDesigner()
        target = _target(sample_target_spec, strand=Strand.PLUS, tss=1000, include_exons=False)
        optimal = designer._score_promoter_position(_spacer(850), 1000, target)
        edge = designer._score_promoter_position(_spacer(630), 1000, target)
        assert optimal > edge


class TestBaseEditorDesigner:
    """Tests for base editor design."""

    def test_cbe_targets_c(self) -> None:
        designer = CBEDesigner()
        with_c = _spacer(100, seq="AAACCAAACCAAACCAAACC")
        without_c = _spacer(100, seq="AAATAAATAAATAAATAAAT")
        kept = designer._filter_by_window([with_c, without_c], max_bystanders=5)
        assert with_c in kept
        assert without_c not in kept

    def test_abe_targets_a(self) -> None:
        designer = ABEDesigner()
        with_a = _spacer(100, seq="CCCAACCCAACCCAACCCAA")
        without_a = _spacer(100, seq="CCCGCCCGCCCGCCCGCCCG")
        kept = designer._filter_by_window([with_a, without_a], max_bystanders=5)
        assert with_a in kept
        assert without_a not in kept

    def test_editing_window_position(self) -> None:
        designer = CBEDesigner()
        strong = designer._score_window_position("AAACCAAACCAAACCAAACC")
        weak = designer._score_window_position("AAAAAAAAAAAAAAAAAAAA")
        assert strong >= weak

    def test_bystander_count(self) -> None:
        designer = CBEDesigner()
        count, positions = designer.count_bystanders("AAACCAAACCAAACCAAACC")
        assert count == len(positions)
        assert all(4 <= p <= 8 for p in positions)


class TestPrimeEditorDesigner:
    """Tests for prime editor design."""

    def test_pbs_length_range(self) -> None:
        designer = PrimeEditorDesigner()
        variants = designer.design_pbs_variants(
            spacer="ATCGATCGATCGATCGATCG",
            target_sequence="A" * 100,
            nick_position=10,
        )
        lengths = {v["length"] for v in variants}
        assert min(lengths) == designer.pbs_length_range[0]
        assert max(lengths) == designer.pbs_length_range[1]

    def test_pbs_tm_calculation(self) -> None:
        tm = PrimeEditorDesigner()._calculate_tm("ATCGATCGAT")
        assert tm > 0

    def test_rtt_encodes_edit(self) -> None:
        params = PrimeEditorDesigner().get_default_parameters()
        assert params["rtt_length_min"] >= 1
        assert params["rtt_length_max"] >= params["rtt_length_min"]

    def test_nicking_guide_distance(self) -> None:
        params = PrimeEditorDesigner().get_default_parameters()
        assert params["pe3_distance_min"] == 40
        assert params["pe3_distance_max"] == 120


class TestPluginRegistry:
    """Tests for modality plugin registration."""

    def test_all_modalities_registered(self) -> None:
        # Import modalities to trigger registration decorators.
        import scalpel.modalities  # noqa: F401

        registered = ModalityRegistry.list_modalities()
        assert EditModality.KNOCKOUT in registered
        assert EditModality.INTERFERENCE in registered
        assert EditModality.ACTIVATION in registered
        assert EditModality.BASE_EDIT_CBE in registered
        assert EditModality.BASE_EDIT_ABE in registered
        assert EditModality.PRIME_EDIT in registered

    def test_get_plugin_instance(self) -> None:
        import scalpel.modalities  # noqa: F401

        ko_designer = ModalityRegistry.get_instance(EditModality.KNOCKOUT)
        assert ko_designer is not None
        assert ko_designer.display_name == "CRISPR Knockout"
