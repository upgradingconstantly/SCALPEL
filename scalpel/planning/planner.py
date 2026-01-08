"""
Experiment Planning Generator.

Creates comprehensive experiment planning documents with controls,
validation assays, and troubleshooting guidance.
"""

from __future__ import annotations

from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
from datetime import datetime
import hashlib
import json

from scalpel.models.enums import EditModality, ValidationTier
from scalpel.models.data_classes import (
    ExperimentPlan,
    ControlRecommendation,
    ValidationAssay,
    FailureMode,
    DesignedGuide,
    GeneInfo,
)


# =============================================================================
# Control Recommendations Database
# =============================================================================

# Validated positive controls for different modalities
POSITIVE_CONTROLS = {
    EditModality.KNOCKOUT: [
        {
            "name": "AAVS1 Safe Harbor",
            "gene": "AAVS1",
            "guide": "GGGGCCACTAGGGACAGGAT",
            "description": "Safe harbor locus with high cutting efficiency",
            "expected_result": "80-95% indel frequency by Sanger/NGS",
        },
        {
            "name": "EMX1 Standard",
            "gene": "EMX1",
            "guide": "GAGTCCGAGCAGAAGAAGAA",
            "description": "Widely validated Cas9 target in EMX1",
            "expected_result": "60-80% indel frequency",
        },
    ],
    EditModality.INTERFERENCE: [
        {
            "name": "PCNA CRISPRi",
            "gene": "PCNA",
            "guide": "GCTGTCACCTTGCGCTCGCT",
            "description": "Essential gene with robust knockdown",
            "expected_result": ">80% knockdown by qPCR",
        },
    ],
    EditModality.ACTIVATION: [
        {
            "name": "TTN CRISPRa",
            "gene": "TTN",
            "guide": "GACGGGTGCGATTTCTTGGT",
            "description": "Large gene with low basal expression",
            "expected_result": ">10-fold activation by qPCR",
        },
    ],
    EditModality.BASE_EDIT_CBE: [
        {
            "name": "HEK3 Site",
            "gene": "HEK3",
            "guide": "GGCCCAGACTGAGCACGTGA",
            "description": "Standard CBE validation site",
            "expected_result": ">50% C>T conversion",
        },
    ],
    EditModality.BASE_EDIT_ABE: [
        {
            "name": "HEK4 Site",
            "gene": "HEK4",
            "guide": "GGCACTGCGGCTGGAGGTGG",
            "description": "Standard ABE validation site",
            "expected_result": ">50% A>G conversion",
        },
    ],
    EditModality.PRIME_EDIT: [
        {
            "name": "HEK3 PE Site",
            "gene": "HEK3",
            "guide": "GGCCCAGACTGAGCACGTGA",
            "description": "Standard prime edit site for PE3",
            "expected_result": ">20% correct edit",
        },
    ],
}

# Non-targeting controls
NEGATIVE_CONTROLS = [
    {
        "name": "Non-targeting Scramble 1",
        "guide": "GCGAGGTATTCGGCTCCGCG",
        "description": "Scrambled sequence with no genomic match",
        "expected_result": "No cutting/editing activity",
        "interpretation": "Baseline for toxicity and off-target effects",
    },
    {
        "name": "Non-targeting Scramble 2",
        "guide": "GTATTACTGATATTGGTGGG",
        "description": "Second non-targeting control for validation",
        "expected_result": "No cutting/editing activity",
        "interpretation": "Confirms specificity of experimental guides",
    },
]


# =============================================================================
# Validation Assay Templates
# =============================================================================

VALIDATION_ASSAYS = {
    # Tier 1: Genotyping
    ValidationTier.TIER1_GENOTYPING: [
        {
            "name": "Sanger Sequencing",
            "description": "Direct sequencing of target locus to confirm editing",
            "timing": "48-72h post-transfection",
            "success_criteria": "Visible mixed peaks at cut site (TIDE analysis)",
            "estimated_cost": "$10-20 per sample",
        },
        {
            "name": "T7 Endonuclease I Assay",
            "description": "Mismatch cleavage assay for quick screening",
            "timing": "48h post-transfection",
            "success_criteria": "Visible cleavage bands on gel",
            "estimated_cost": "$5-10 per sample",
        },
        {
            "name": "NGS Amplicon Sequencing",
            "description": "Deep sequencing for precise editing quantification",
            "timing": "48-72h post-transfection",
            "success_criteria": ">30% indel frequency",
            "estimated_cost": "$50-100 per sample",
        },
    ],
    # Tier 2: Molecular Validation
    ValidationTier.TIER2_MOLECULAR: [
        {
            "name": "qRT-PCR",
            "description": "Quantify mRNA knockdown/activation",
            "timing": "48-72h for transient, 7d for stable",
            "success_criteria": ">70% knockdown or >2-fold activation",
            "estimated_cost": "$20-50 per sample",
        },
        {
            "name": "Western Blot",
            "description": "Protein level validation",
            "timing": "72h-7d depending on protein half-life",
            "success_criteria": "Visible reduction/absence of target band",
            "estimated_cost": "$50-100 per sample",
        },
    ],
    # Tier 3: Functional Validation
    ValidationTier.TIER3_FUNCTIONAL: [
        {
            "name": "Cell Viability Assay",
            "description": "MTT/CellTiter-Glo for growth effects",
            "timing": "3-7 days post-editing",
            "success_criteria": "Statistically significant change from control",
            "estimated_cost": "$10-30 per sample",
        },
        {
            "name": "Phenotypic Assay",
            "description": "Modality-specific functional readout",
            "timing": "Depends on phenotype",
            "success_criteria": "Gene-specific functional change",
            "estimated_cost": "Variable",
        },
    ],
}


# =============================================================================
# Failure Modes
# =============================================================================

COMMON_FAILURE_MODES = [
    {
        "failure": "No editing detected",
        "probability": "Medium",
        "symptoms": ["No indels in sequencing", "Wildtype bands only"],
        "potential_causes": [
            "Low transfection efficiency",
            "Inactive guide",
            "Cas9/gRNA not expressed",
        ],
        "troubleshooting": [
            "Check transfection with GFP control",
            "Verify guide expression by qPCR",
            "Try different delivery method",
            "Use validated positive control guide",
        ],
    },
    {
        "failure": "High toxicity",
        "probability": "Low",
        "symptoms": ["Massive cell death", "Poor recovery"],
        "potential_causes": [
            "Off-target editing",
            "Essential gene knockout",
            "Delivery toxicity",
        ],
        "troubleshooting": [
            "Reduce Cas9/gRNA dose",
            "Use non-integrating delivery",
            "Check off-target profile",
            "Consider CRISPRi instead of knockout",
        ],
    },
    {
        "failure": "No phenotype despite editing",
        "probability": "Medium",
        "symptoms": ["High editing rate", "No functional change"],
        "potential_causes": [
            "In-frame indels",
            "Functional redundancy",
            "NMD escape",
            "Alternative splicing",
        ],
        "troubleshooting": [
            "Target earlier in CDS",
            "Design multiple guides",
            "Check protein by Western",
            "Consider functional domain targeting",
        ],
    },
]


class ExperimentPlanner:
    """
    Generates comprehensive experiment planning documents.
    
    Creates recommendations for:
    - Positive and negative controls
    - Tiered validation assays
    - Troubleshooting guides
    - Sequencing strategies
    """
    
    def __init__(self, modality: EditModality):
        self.modality = modality
    
    def generate_plan(
        self,
        guides: List[DesignedGuide],
        gene_info: Optional[GeneInfo] = None,
        experiment_context: Optional[Dict[str, Any]] = None,
    ) -> ExperimentPlan:
        """
        Generate a complete experiment plan.
        
        Args:
            guides: Designed guides to include in plan
            gene_info: Information about target gene
            experiment_context: Additional context (cell type, delivery, etc.)
        
        Returns:
            ExperimentPlan with all recommendations
        """
        target_gene = gene_info.symbol if gene_info else "Unknown"
        context = experiment_context or {}
        
        # Get controls
        positive_controls = self._get_positive_controls()
        negative_controls = self._get_negative_controls()
        
        # Get validation assays
        validation_assays = self._get_validation_assays(context)
        
        # Get failure modes
        failure_modes = self._get_failure_modes()
        
        # Generate provenance hash
        input_hash = self._generate_hash(guides, gene_info)
        
        return ExperimentPlan(
            target_gene=target_gene,
            edit_modality=self.modality,
            selected_guides=guides,
            positive_controls=positive_controls,
            negative_controls=negative_controls,
            validation_assays=validation_assays,
            failure_modes=failure_modes,
            input_hash=input_hash,
        )
    
    def _get_positive_controls(self) -> List[ControlRecommendation]:
        """Get modality-appropriate positive controls."""
        controls = []
        
        modality_controls = POSITIVE_CONTROLS.get(self.modality, [])
        for ctrl in modality_controls[:2]:  # Limit to 2
            controls.append(ControlRecommendation(
                control_type="positive",
                name=ctrl["name"],
                description=ctrl["description"],
                gene_target=ctrl.get("gene"),
                guide_sequence=ctrl.get("guide"),
                expected_result=ctrl["expected_result"],
                interpretation=f"Validates {self.modality.value} activity in your system",
            ))
        
        return controls
    
    def _get_negative_controls(self) -> List[ControlRecommendation]:
        """Get non-targeting negative controls."""
        controls = []
        
        for ctrl in NEGATIVE_CONTROLS[:2]:  # Limit to 2
            controls.append(ControlRecommendation(
                control_type="negative",
                name=ctrl["name"],
                description=ctrl["description"],
                guide_sequence=ctrl["guide"],
                expected_result=ctrl["expected_result"],
                interpretation=ctrl["interpretation"],
            ))
        
        return controls
    
    def _get_validation_assays(
        self,
        context: Dict[str, Any],
    ) -> List[ValidationAssay]:
        """Get tiered validation assays."""
        assays = []
        
        # Add Tier 1 (always required)
        for assay in VALIDATION_ASSAYS[ValidationTier.TIER1_GENOTYPING][:2]:
            assays.append(ValidationAssay(
                tier=1,
                name=assay["name"],
                description=assay["description"],
                timing=assay["timing"],
                success_criteria=assay["success_criteria"],
                estimated_cost=assay.get("estimated_cost"),
            ))
        
        # Add Tier 2
        for assay in VALIDATION_ASSAYS[ValidationTier.TIER2_MOLECULAR][:2]:
            assays.append(ValidationAssay(
                tier=2,
                name=assay["name"],
                description=assay["description"],
                timing=assay["timing"],
                success_criteria=assay["success_criteria"],
                estimated_cost=assay.get("estimated_cost"),
            ))
        
        # Add Tier 3 (for full validation)
        for assay in VALIDATION_ASSAYS[ValidationTier.TIER3_FUNCTIONAL][:1]:
            assays.append(ValidationAssay(
                tier=3,
                name=assay["name"],
                description=assay["description"],
                timing=assay["timing"],
                success_criteria=assay["success_criteria"],
                estimated_cost=assay.get("estimated_cost"),
            ))
        
        return assays
    
    def _get_failure_modes(self) -> List[FailureMode]:
        """Get common failure modes and troubleshooting."""
        modes = []
        
        for fm in COMMON_FAILURE_MODES:
            modes.append(FailureMode(
                failure=fm["failure"],
                probability=fm["probability"],
                symptoms=fm["symptoms"],
                potential_causes=fm["potential_causes"],
                troubleshooting=fm["troubleshooting"],
            ))
        
        return modes
    
    def _generate_hash(
        self,
        guides: List[DesignedGuide],
        gene_info: Optional[GeneInfo],
    ) -> str:
        """Generate provenance hash for reproducibility."""
        data = {
            "modality": self.modality.value,
            "gene": gene_info.symbol if gene_info else None,
            "guides": [g.spacer.spacer_sequence for g in guides[:5]],
            "timestamp": datetime.utcnow().isoformat(),
        }
        content = json.dumps(data, sort_keys=True)
        return hashlib.sha256(content.encode()).hexdigest()[:16]
    
    def to_markdown(self, plan: ExperimentPlan) -> str:
        """Convert plan to Markdown format."""
        lines = [
            f"# Experiment Plan: {plan.target_gene} {plan.edit_modality.value.title()}",
            "",
            f"**Generated:** {plan.created_at.strftime('%Y-%m-%d %H:%M')}",
            f"**Provenance:** `{plan.input_hash}`",
            "",
            "---",
            "",
            "## Selected Guides",
            "",
        ]
        
        for i, guide in enumerate(plan.selected_guides[:5], 1):
            lines.append(f"### Guide {i}: `{guide.spacer.spacer_sequence}`")
            lines.append(f"- PAM: {guide.spacer.pam_sequence}")
            lines.append(f"- Efficiency: {guide.efficiency_score.overall_score:.2f}")
            lines.append(f"- Position: chr{guide.spacer.genomic_start}")
            lines.append("")
        
        lines.extend([
            "---",
            "",
            "## Controls",
            "",
            "### Positive Controls",
            "",
        ])
        
        for ctrl in plan.positive_controls:
            lines.append(f"**{ctrl.name}** ({ctrl.gene_target})")
            lines.append(f"- Guide: `{ctrl.guide_sequence}`")
            lines.append(f"- Expected: {ctrl.expected_result}")
            lines.append("")
        
        lines.extend([
            "### Negative Controls",
            "",
        ])
        
        for ctrl in plan.negative_controls:
            lines.append(f"**{ctrl.name}**")
            lines.append(f"- Guide: `{ctrl.guide_sequence}`")
            lines.append(f"- Purpose: {ctrl.interpretation}")
            lines.append("")
        
        lines.extend([
            "---",
            "",
            "## Validation Assays",
            "",
        ])
        
        for tier in [1, 2, 3]:
            tier_assays = [a for a in plan.validation_assays if a.tier == tier]
            if tier_assays:
                tier_names = {1: "Genotyping", 2: "Molecular", 3: "Functional"}
                lines.append(f"### Tier {tier}: {tier_names.get(tier, 'Other')}")
                lines.append("")
                for assay in tier_assays:
                    lines.append(f"**{assay.name}**")
                    lines.append(f"- Timing: {assay.timing}")
                    lines.append(f"- Success: {assay.success_criteria}")
                    if assay.estimated_cost:
                        lines.append(f"- Cost: {assay.estimated_cost}")
                    lines.append("")
        
        lines.extend([
            "---",
            "",
            "## Troubleshooting",
            "",
        ])
        
        for fm in plan.failure_modes:
            lines.append(f"### ⚠️ {fm.failure}")
            lines.append(f"**Probability:** {fm.probability}")
            lines.append("")
            lines.append("**Symptoms:**")
            for s in fm.symptoms:
                lines.append(f"- {s}")
            lines.append("")
            lines.append("**Troubleshooting:**")
            for t in fm.troubleshooting[:3]:
                lines.append(f"1. {t}")
            lines.append("")
        
        return "\n".join(lines)


def generate_experiment_plan(
    guides: List[DesignedGuide],
    modality: EditModality,
    gene_info: Optional[GeneInfo] = None,
) -> ExperimentPlan:
    """Convenience function to generate experiment plan."""
    planner = ExperimentPlanner(modality)
    return planner.generate_plan(guides, gene_info)
