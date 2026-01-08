"""
SCALPEL CLI - Scriptable CRISPR design workflows.

Usage:
    scalpel design --gene TP53 --modality knockout --genome GRCh38
    scalpel offtarget --spacer ATCGATCGATCGATCGATCG
    scalpel plan --input results.json --output plan.md
    
Pipe-friendly:
    scalpel design --gene BRCA1 | scalpel offtarget | scalpel plan
"""

import json
import sys
from pathlib import Path
from typing import Optional, List
from enum import Enum

import typer
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn

from scalpel import __version__
from scalpel.models.enums import Genome, EditModality, CasVariantType
from scalpel.config import get_config

# Initialize Typer app and Rich console
app = typer.Typer(
    name="scalpel",
    help="Computational CRISPR design platform",
    no_args_is_help=True,
    rich_markup_mode="rich",
)
console = Console()


# =============================================================================
# Version callback
# =============================================================================

def version_callback(value: bool) -> None:
    if value:
        console.print(f"[bold blue]SCALPEL[/bold blue] version {__version__}")
        raise typer.Exit()


# =============================================================================
# Main app options
# =============================================================================

@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        help="Show version and exit",
        callback=version_callback,
        is_eager=True,
    ),
) -> None:
    """
    SCALPEL - Computational CRISPR Design Platform
    
    Design gRNAs, analyze off-targets, and generate experiment plans.
    """
    pass


# =============================================================================
# Design command
# =============================================================================

@app.command()
def design(
    gene: Optional[str] = typer.Option(
        None,
        "--gene", "-g",
        help="Gene symbol to target (e.g., TP53, BRCA1)",
    ),
    sequence: Optional[str] = typer.Option(
        None,
        "--sequence", "-s",
        help="Target DNA sequence directly",
    ),
    coordinates: Optional[str] = typer.Option(
        None,
        "--coordinates", "-c",
        help="Genomic coordinates (chr:start-end)",
    ),
    modality: EditModality = typer.Option(
        EditModality.KNOCKOUT,
        "--modality", "-m",
        help="Editing modality",
    ),
    genome: str = typer.Option(
        "GRCh38",
        "--genome",
        help="Reference genome (GRCh38, GRCh37, GRCm39, etc.)",
    ),
    cas_variant: CasVariantType = typer.Option(
        CasVariantType.SPCAS9,
        "--cas",
        help="Cas protein variant",
    ),
    n_guides: int = typer.Option(
        10,
        "--n-guides", "-n",
        help="Number of guides to return",
    ),
    stream: bool = typer.Option(
        False,
        "--stream",
        help="Stream results progressively (JSONL format, one guide per line)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Output file (JSON). If not specified, prints to stdout.",
    ),
    format: str = typer.Option(
        "json",
        "--format", "-f",
        help="Output format: json, tsv, table",
    ),
) -> None:
    """
    Design gRNAs for a target gene or sequence.
    
    Examples:
        scalpel design --gene TP53 --modality knockout
        scalpel design --sequence ATCG... --cas SpCas9-NG
        scalpel design --coordinates chr17:7668421-7687490
    """
    # Validate input
    if not any([gene, sequence, coordinates]):
        console.print("[red]Error:[/red] Must specify --gene, --sequence, or --coordinates")
        raise typer.Exit(1)
    
    # Parse genome
    try:
        genome_enum = Genome.from_string(genome)
    except ValueError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        # Resolve target
        task = progress.add_task("Resolving target...", total=None)
        
        # Import target resolver
        from scalpel.genome.target_resolver import TargetResolver, TargetNotFoundError
        from scalpel.models.data_classes import TargetSpecification
        from scalpel.models.enums import TargetType
        
        try:
            # Build target specification
            if gene:
                target_type = TargetType.GENE_SYMBOL
                target_value = gene
            elif sequence:
                target_type = TargetType.SEQUENCE
                target_value = sequence
            else:
                target_type = TargetType.GENOMIC_COORDINATES
                target_value = coordinates
            
            spec = TargetSpecification(
                target_type=target_type,
                target_value=target_value,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant,
            )
            
            # Resolve target
            resolver = TargetResolver(genome_enum)
            resolved = resolver.resolve(spec)
            
            progress.update(task, description="Extracting spacers...")
            
            # Extract spacers using the Stage 3 engine
            from scalpel.design import SpacerExtractor
            
            extractor = SpacerExtractor(cas_variant)
            
            # Get sequence for extraction
            # For demo genes, we use a synthetic sequence with PAM sites
            if resolved.sequence and "N" not in resolved.sequence[:100]:
                target_sequence = resolved.sequence
            else:
                # Generate demo sequence with PAM sites for testing
                target_sequence = _generate_demo_sequence(resolved.end - resolved.start)
            
            spacers = extractor.extract_spacers(
                target_sequence,
                chromosome=resolved.chromosome,
                start_position=resolved.start,
            )
            
            progress.update(task, description=f"Found {len(spacers)} candidates, scoring...")
            
            # Score spacers for efficiency (Stage 4)
            from scalpel.design.efficiency import EnsembleScorer
            
            scorer = EnsembleScorer()
            scored_guides = scorer.score_batch(spacers, modality.value, n_top=n_guides)
            
            progress.update(task, description=f"Ranked {len(scored_guides)} guides!")
            
            # Convert to output format
            guide_output = []
            for i, sg in enumerate(scored_guides, 1):
                guide_output.append({
                    "rank": i,
                    "spacer_sequence": sg.spacer.spacer_sequence,
                    "pam_sequence": sg.spacer.pam_sequence,
                    "strand": sg.spacer.strand.value,
                    "genomic_start": sg.spacer.genomic_start,
                    "genomic_end": sg.spacer.genomic_end,
                    "cut_site": sg.spacer.cut_site,
                    "efficiency_score": round(sg.efficiency.overall_score, 3),
                    "efficiency_interpretation": sg.efficiency.interpretation,
                    "off_target_count": None,  # Will be filled by offtarget command
                })
            
            # Add red flag detection (Stage 7)
            from scalpel.core.red_flags import detect_red_flags, summarize_red_flags
            
            for guide in guide_output:
                # Find the corresponding spacer
                spacer_seq = guide["spacer_sequence"]
                matching_spacer = next(
                    (sg.spacer for sg in scored_guides if sg.spacer.spacer_sequence == spacer_seq),
                    None
                )
                if matching_spacer:
                    flags = detect_red_flags(matching_spacer, gene_info=resolved.gene_info)
                    flag_summary = summarize_red_flags(flags)
                    guide["red_flags"] = flag_summary
            
            progress.update(task, description="Design complete!")
            
            # Build results
            results = {
                "status": "success",
                "target": {
                    "gene": resolved.gene_info.symbol if resolved.gene_info else None,
                    "gene_id": resolved.gene_info.gene_id if resolved.gene_info else None,
                    "chromosome": resolved.chromosome,
                    "start": resolved.start,
                    "end": resolved.end,
                    "strand": resolved.strand.value,
                    "genome": genome_enum.value,
                    "modality": modality.value,
                    "cas_variant": cas_variant.value,
                    "sequence_length": len(target_sequence),
                },
                "n_guides_requested": n_guides,
                "n_guides_found": len(spacers),
                "guides": guide_output,
            }
            
            # Add transcript info if available
            if resolved.transcript:
                results["target"]["transcript_id"] = resolved.transcript.transcript_id
                results["target"]["is_mane_select"] = resolved.transcript.is_mane_select
                if resolved.transcript.exons:
                    results["target"]["n_exons"] = len(resolved.transcript.exons)
        
        except TargetNotFoundError as e:
            results = {
                "status": "error",
                "error": str(e),
                "target": {
                    "gene": gene,
                    "sequence": sequence,
                    "coordinates": coordinates,
                },
                "guides": [],
            }
        except Exception as e:
            results = {
                "status": "error",
                "error": f"Unexpected error: {e}",
                "guides": [],
            }
    
    # Output
    if stream:
        # Streaming mode: output JSONL (one JSON object per line)
        # First output target info
        print(json.dumps({"type": "init", "target": results.get("target", {})}), flush=True)
        # Then output guides one by one
        for guide in results.get("guides", []):
            print(json.dumps({"type": "guide", **guide}), flush=True)
        # Finally output done
        print(json.dumps({"type": "done", "n_guides_found": results.get("n_guides_found", 0)}), flush=True)
    else:
        _output_results(results, output, format)


# =============================================================================
# Off-target command
# =============================================================================

@app.command()
def offtarget(
    spacer: Optional[str] = typer.Option(
        None,
        "--spacer", "-s",
        help="20bp spacer sequence to analyze",
    ),
    input_file: Optional[Path] = typer.Option(
        None,
        "--input", "-i",
        help="Input JSON from design command (reads from stdin if -)",
    ),
    genome: str = typer.Option(
        "GRCh38",
        "--genome",
        help="Reference genome",
    ),
    max_mismatches: int = typer.Option(
        4,
        "--max-mismatches",
        help="Maximum mismatches to consider",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Output file (JSON)",
    ),
    format: str = typer.Option(
        "json",
        "--format", "-f",
        help="Output format: json, tsv, table",
    ),
) -> None:
    """
    Analyze off-target sites for a spacer sequence.
    
    Examples:
        scalpel offtarget --spacer ATCGATCGATCGATCGATCG
        scalpel design --gene TP53 | scalpel offtarget
    """
    # Handle piped input
    spacers_to_analyze: List[str] = []
    
    if spacer:
        spacers_to_analyze.append(spacer)
    elif input_file:
        if str(input_file) == "-":
            # Read from stdin
            data = json.load(sys.stdin)
        else:
            with open(input_file) as f:
                data = json.load(f)
        
        # Extract spacers from design results
        if "guides" in data:
            for guide in data["guides"]:
                if "spacer_sequence" in guide:
                    spacers_to_analyze.append(guide["spacer_sequence"])
    elif not sys.stdin.isatty():
        # Stdin is piped
        data = json.load(sys.stdin)
        if "guides" in data:
            for guide in data["guides"]:
                if "spacer_sequence" in guide:
                    spacers_to_analyze.append(guide["spacer_sequence"])
    
    if not spacers_to_analyze:
        console.print("[red]Error:[/red] No spacer sequences provided")
        raise typer.Exit(1)
    
    # Parse genome
    try:
        genome_enum = Genome.from_string(genome)
    except ValueError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)
    
    # Run off-target analysis
    from scalpel.offtarget import OffTargetSearcher, RiskCalculator
    
    searcher = OffTargetSearcher(genome_enum, max_mismatches=max_mismatches)
    risk_calc = RiskCalculator(genome_enum)
    
    analyses = []
    for spacer_seq in spacers_to_analyze:
        # Search for off-targets
        analysis = searcher.search(spacer_seq)
        
        # Generate summary
        summary = risk_calc.generate_summary(analysis)
        
        # Convert to output format
        site_output = []
        for site in analysis.sites[:20]:  # Limit to top 20
            site_output.append({
                "chromosome": site.chromosome,
                "position": site.position,
                "strand": site.strand.value,
                "sequence": site.sequence,
                "pam": site.pam,
                "mismatches": site.mismatch_count,
                "cfd_score": round(site.cutting_probability, 4),
                "risk_score": round(site.risk_score, 3),
            })
        
        analyses.append({
            "spacer": spacer_seq,
            "total_sites": summary["total_sites"],
            "high_risk_count": summary["high_risk_count"],
            "specificity_score": summary["specificity_score"],
            "interpretation": summary["interpretation"],
            "sites_by_mismatch": summary["sites_by_mismatch"],
            "sites": site_output,
        })
    
    results = {
        "status": "success",
        "genome": genome_enum.value,
        "max_mismatches": max_mismatches,
        "analyses": analyses,
    }
    
    _output_results(results, output, format)


# =============================================================================
# Plan command
# =============================================================================

@app.command()
def plan(
    input_file: Optional[Path] = typer.Option(
        None,
        "--input", "-i",
        help="Input JSON from design/offtarget commands (reads from stdin if -)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Output file (Markdown)",
    ),
    format: str = typer.Option(
        "markdown",
        "--format", "-f",
        help="Output format: markdown, html, json",
    ),
) -> None:
    """
    Generate an experiment plan from design results.
    
    Examples:
        scalpel plan --input design_results.json --output plan.md
        scalpel design --gene TP53 | scalpel offtarget | scalpel plan
    """
    # Load input
    if input_file:
        if str(input_file) == "-":
            data = json.load(sys.stdin)
        else:
            with open(input_file) as f:
                data = json.load(f)
    elif not sys.stdin.isatty():
        data = json.load(sys.stdin)
    else:
        console.print("[red]Error:[/red] No input provided")
        raise typer.Exit(1)
    
    # Parse input and generate plan
    from scalpel.planning import ExperimentPlanner
    from scalpel.models.data_classes import SpacerCandidate, EfficiencyScore, DesignedGuide, GeneInfo
    from scalpel.models.enums import EditModality, Strand
    
    # Get modality from input
    modality_str = data.get("target", {}).get("modality", "knockout")
    try:
        modality = EditModality(modality_str)
    except ValueError:
        modality = EditModality.KNOCKOUT
    
    # Create planner
    planner = ExperimentPlanner(modality)
    
    # Convert guides from JSON to DesignedGuide objects
    guides = []
    for guide_data in data.get("guides", [])[:5]:
        spacer = SpacerCandidate(
            spacer_sequence=guide_data.get("spacer_sequence", ""),
            pam_sequence=guide_data.get("pam_sequence", ""),
            strand=Strand(guide_data.get("strand", "+")),
            genomic_start=guide_data.get("genomic_start", 0),
            genomic_end=guide_data.get("genomic_end", 0),
            cut_site=guide_data.get("cut_site", 0),
            context_sequence="",
        )
        efficiency = EfficiencyScore(
            overall_score=guide_data.get("efficiency_score", 0.5),
            interpretation=guide_data.get("efficiency_interpretation", ""),
        )
        guides.append(DesignedGuide(
            spacer=spacer,
            efficiency_score=efficiency,
            composite_score=guide_data.get("efficiency_score", 0.5),
        ))
    
    # Get gene info
    gene_info = None
    target = data.get("target", {})
    if target.get("gene"):
        gene_info = GeneInfo(
            gene_id=target.get("gene_id", ""),
            symbol=target.get("gene", ""),
            chromosome=target.get("chromosome", ""),
            start=target.get("start", 0),
            end=target.get("end", 0),
            strand=Strand(target.get("strand", "+")),
        )
    
    # Generate plan
    plan = planner.generate_plan(guides, gene_info)
    
    # Output
    if format.lower() == "json":
        plan_output = {
            "target_gene": plan.target_gene,
            "modality": plan.edit_modality.value,
            "n_guides": len(plan.selected_guides),
            "positive_controls": [
                {"name": c.name, "guide": c.guide_sequence, "expected": c.expected_result}
                for c in plan.positive_controls
            ],
            "negative_controls": [
                {"name": c.name, "guide": c.guide_sequence}
                for c in plan.negative_controls
            ],
            "validation_assays": [
                {"tier": a.tier, "name": a.name, "timing": a.timing}
                for a in plan.validation_assays
            ],
            "failure_modes": [
                {"failure": f.failure, "probability": f.probability}
                for f in plan.failure_modes
            ],
            "provenance": plan.input_hash,
        }
        import json as json_mod
        plan_content = json_mod.dumps(plan_output, indent=2)
    else:
        plan_content = planner.to_markdown(plan)
    
    if output:
        output.write_text(plan_content)
        console.print(f"[green]Plan written to {output}[/green]")
    else:
        console.print(plan_content)


# =============================================================================
# Info command
# =============================================================================

@app.command()
def info() -> None:
    """Show configuration and available resources."""
    config = get_config()
    
    console.print(Panel.fit(
        f"[bold blue]SCALPEL[/bold blue] v{__version__}\n"
        f"Data directory: {config.data_dir}\n"
        f"Genomes directory: {config.genomes_dir}\n"
        f"Models directory: {config.models_dir}",
        title="Configuration",
    ))
    
    # Available modalities
    table = Table(title="Available Modalities")
    table.add_column("Modality", style="cyan")
    table.add_column("Description")
    
    for modality in EditModality:
        table.add_row(modality.value, modality.name.replace("_", " ").title())
    
    console.print(table)
    
    # Available Cas variants
    table = Table(title="Supported Cas Variants")
    table.add_column("Variant", style="green")
    table.add_column("PAM")
    
    cas_info = {
        CasVariantType.SPCAS9: "NGG",
        CasVariantType.SPCAS9_NG: "NG",
        CasVariantType.SPCAS9_VQR: "NGA",
        CasVariantType.SACAS9: "NNGRRT",
        CasVariantType.CAS12A: "TTTV",
        CasVariantType.CAS12A_RR: "TYCV",
    }
    
    for variant, pam in cas_info.items():
        table.add_row(variant.value, pam)
    
    console.print(table)


# =============================================================================
# Utility functions
# =============================================================================

def _generate_demo_sequence(length: int) -> str:
    """
    Generate a demo DNA sequence with PAM sites for testing.
    
    Creates a realistic sequence with:
    - ~50% GC content
    - Multiple NGG PAM sites
    - No long homopolymers
    """
    import random
    random.seed(42)  # Reproducible for demos
    
    # Build sequence with periodic PAM sites
    sequence = []
    pam_interval = 50  # Insert PAM approximately every 50bp
    
    for i in range(min(length, 10000)):  # Cap at 10kb
        if i > 20 and i % pam_interval == 0:
            # Insert an NGG PAM site
            sequence.extend(["A", "G", "G"])
            i += 2
        else:
            # Random base with ~50% GC
            base = random.choice(["A", "C", "G", "T"])
            sequence.append(base)
    
    return "".join(sequence[:length])


def _output_results(results: dict, output: Optional[Path], format: str) -> None:
    """Output results in specified format."""
    if format == "json":
        json_str = json.dumps(results, indent=2, default=str)
        if output:
            output.write_text(json_str)
            console.print(f"[green]Results written to {output}[/green]")
        else:
            # Print to stdout for piping
            print(json_str)
    
    elif format == "tsv":
        # Convert to TSV (simplified)
        lines = []
        if "guides" in results:
            headers = ["rank", "spacer", "pam", "score"]
            lines.append("\t".join(headers))
            for i, guide in enumerate(results["guides"], 1):
                row = [
                    str(i),
                    guide.get("spacer_sequence", ""),
                    guide.get("pam_sequence", ""),
                    str(guide.get("composite_score", "")),
                ]
                lines.append("\t".join(row))
        
        tsv_str = "\n".join(lines)
        if output:
            output.write_text(tsv_str)
        else:
            print(tsv_str)
    
    elif format == "table":
        # Rich table output
        if "guides" in results and results["guides"]:
            table = Table(title="Designed Guides")
            table.add_column("Rank", style="dim")
            table.add_column("Spacer", style="cyan")
            table.add_column("PAM", style="green")
            table.add_column("Score", style="yellow")
            
            for i, guide in enumerate(results["guides"], 1):
                table.add_row(
                    str(i),
                    guide.get("spacer_sequence", ""),
                    guide.get("pam_sequence", ""),
                    f"{guide.get('composite_score', 0):.3f}",
                )
            
            console.print(table)
        else:
            console.print(json.dumps(results, indent=2, default=str))


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    app()
