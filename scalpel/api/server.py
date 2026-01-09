"""
SCALPEL Web API

FastAPI-based web interface for the CRISPR design platform.
"""

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List, AsyncGenerator
from pathlib import Path
import uvicorn
import asyncio
import json

app = FastAPI(
    title="SCALPEL",
    description="Computational CRISPR Design Platform",
    version="0.1.0",
)

# CORS for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# =============================================================================
# Request/Response Models
# =============================================================================

class DesignRequest(BaseModel):
    gene: Optional[str] = None
    coordinates: Optional[str] = None
    sequence: Optional[str] = None
    modality: str = "knockout"
    genome: str = "GRCh38"
    cas_variant: str = "SpCas9"
    n_guides: int = 10


class OffTargetRequest(BaseModel):
    spacer: str
    genome: str = "GRCh38"
    max_mismatches: int = 4


class PlanRequest(BaseModel):
    design_data: dict


# =============================================================================
# API Endpoints
# =============================================================================

@app.get("/", response_class=HTMLResponse)
async def root():
    """Serve the main web interface."""
    index_file = STATIC_DIR / "index.html"
    if index_file.exists():
        return FileResponse(str(index_file), media_type="text/html")
    return get_html_ui()


@app.get("/api/health")
async def health():
    """Health check endpoint."""
    return {"status": "healthy", "version": "0.1.0"}


@app.post("/api/design")
async def design_guides(request: DesignRequest):
    """Design gRNAs for a target."""
    try:
        from scalpel.genome import TargetResolver, GeneDatabase, get_demo_gene
        from scalpel.genome.target_resolver import TargetNotFoundError
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        from scalpel.core.red_flags import detect_red_flags, summarize_red_flags
        from scalpel.models.enums import Genome, EditModality, CasVariantType, TargetType
        from scalpel.models.data_classes import TargetSpecification
        
        # Parse enums
        genome_enum = Genome.from_string(request.genome)
        modality = EditModality(request.modality)
        cas_variant = CasVariantType(request.cas_variant)
        
        # Resolve target
        resolver = TargetResolver(genome_enum)
        
        if request.gene:
            spec = TargetSpecification(
                target_type=TargetType.GENE_SYMBOL,
                target_value=request.gene,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        elif request.coordinates:
            spec = TargetSpecification(
                target_type=TargetType.GENOMIC_COORDINATES,
                target_value=request.coordinates,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        elif request.sequence:
            spec = TargetSpecification(
                target_type=TargetType.SEQUENCE,
                target_value=request.sequence,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        else:
            raise HTTPException(status_code=400, detail="Must provide gene, coordinates, or sequence")
        
        # Get or generate sequence
        target_sequence = resolved.sequence
        # Generate demo sequence if sequence is empty, too short, or all N's (placeholder)
        if not target_sequence or len(target_sequence) < 50 or set(target_sequence) == {'N'}:
            # Generate demo sequence
            import random
            random.seed(hash(resolved.gene_info.symbol if resolved.gene_info else "TEST"))
            bases = "ACGT"
            target_sequence = ""
            for i in range(resolved.end - resolved.start):
                if i % 50 == 0 and i > 0:
                    pam = random.choice(["AGG", "TGG", "CGG", "GGG"])
                    target_sequence += pam
                else:
                    target_sequence += random.choice(bases)
        
        # Extract spacers
        extractor = SpacerExtractor(cas_variant)
        spacers = extractor.extract_spacers(
            target_sequence,
            chromosome=resolved.chromosome,
            start_position=resolved.start,
        )
        
        # Score
        scorer = EnsembleScorer()
        scored_guides = scorer.score_batch(spacers, modality.value, n_top=None)  # Return ALL candidates
        
        # Build output
        guides = []
        for i, sg in enumerate(scored_guides, 1):
            flags = detect_red_flags(sg.spacer, gene_info=resolved.gene_info)
            flag_summary = summarize_red_flags(flags)
            
            guides.append({
                "rank": i,
                "spacer_sequence": sg.spacer.spacer_sequence,
                "pam_sequence": sg.spacer.pam_sequence,
                "strand": sg.spacer.strand.value,
                "genomic_start": sg.spacer.genomic_start,
                "genomic_end": sg.spacer.genomic_end,
                "cut_site": sg.spacer.cut_site,
                "efficiency_score": round(sg.efficiency.overall_score, 3),
                "efficiency_interpretation": sg.efficiency.interpretation,
                "red_flags": flag_summary,
            })
        
        return {
            "status": "success",
            "target": {
                "gene": resolved.gene_info.symbol if resolved.gene_info else None,
                "chromosome": resolved.chromosome,
                "start": resolved.start,
                "end": resolved.end,
                "genome": genome_enum.value,
                "modality": modality.value,
            },
            "n_guides_found": len(spacers),
            "guides": guides,
        }
        
    except TargetNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


async def design_stream_generator(request: DesignRequest) -> AsyncGenerator[str, None]:
    """Async generator that yields guides as Server-Sent Events."""
    try:
        from scalpel.genome import TargetResolver, GeneDatabase, get_demo_gene
        from scalpel.genome.target_resolver import TargetNotFoundError
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        from scalpel.core.red_flags import detect_red_flags, summarize_red_flags
        from scalpel.models.enums import Genome, EditModality, CasVariantType, TargetType
        from scalpel.models.data_classes import TargetSpecification
        
        # Parse enums
        genome_enum = Genome.from_string(request.genome)
        modality = EditModality(request.modality)
        cas_variant = CasVariantType(request.cas_variant)
        
        # Resolve target
        resolver = TargetResolver(genome_enum)
        
        if request.gene:
            spec = TargetSpecification(
                target_type=TargetType.GENE_SYMBOL,
                target_value=request.gene,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        elif request.coordinates:
            spec = TargetSpecification(
                target_type=TargetType.GENOMIC_COORDINATES,
                target_value=request.coordinates,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        elif request.sequence:
            spec = TargetSpecification(
                target_type=TargetType.SEQUENCE,
                target_value=request.sequence,
                genome=genome_enum,
                modality=modality,
                cas_variant=cas_variant
            )
            resolved = resolver.resolve(spec)
        else:
            yield f"event: error\ndata: {json.dumps({'error': 'Must provide gene, coordinates, or sequence'})}\n\n"
            return
        
        # Send init event with target info
        init_data = {
            "type": "init",
            "target": {
                "gene": resolved.gene_info.symbol if resolved.gene_info else None,
                "chromosome": resolved.chromosome,
                "start": resolved.start,
                "end": resolved.end,
                "genome": genome_enum.value,
                "modality": modality.value,
            }
        }
        yield f"event: init\ndata: {json.dumps(init_data)}\n\n"
        await asyncio.sleep(0.05)  # Small delay for client to render
        
        # Get or generate sequence
        target_sequence = resolved.sequence
        # Generate demo sequence if sequence is empty, too short, or all N's (placeholder)
        if not target_sequence or len(target_sequence) < 50 or set(target_sequence) == {'N'}:
            import random
            random.seed(hash(resolved.gene_info.symbol if resolved.gene_info else "TEST"))
            bases = "ACGT"
            target_sequence = ""
            for i in range(resolved.end - resolved.start):
                if i % 50 == 0 and i > 0:
                    pam = random.choice(["AGG", "TGG", "CGG", "GGG"])
                    target_sequence += pam
                else:
                    target_sequence += random.choice(bases)
        
        # Extract spacers
        extractor = SpacerExtractor(cas_variant)
        spacers = extractor.extract_spacers(
            target_sequence,
            chromosome=resolved.chromosome,
            start_position=resolved.start,
        )
        
        # Send progress event
        yield f"event: progress\ndata: {json.dumps({'stage': 'scoring', 'total_spacers': len(spacers)})}\n\n"
        await asyncio.sleep(0.05)
        
        # Score and stream guides one by one
        scorer = EnsembleScorer()
        scored_guides = scorer.score_batch(spacers, modality.value, n_top=None)  # Return ALL candidates
        
        for i, sg in enumerate(scored_guides, 1):
            flags = detect_red_flags(sg.spacer, gene_info=resolved.gene_info)
            flag_summary = summarize_red_flags(flags)
            
            guide_data = {
                "rank": i,
                "spacer_sequence": sg.spacer.spacer_sequence,
                "pam_sequence": sg.spacer.pam_sequence,
                "strand": sg.spacer.strand.value,
                "genomic_start": sg.spacer.genomic_start,
                "genomic_end": sg.spacer.genomic_end,
                "cut_site": sg.spacer.cut_site,
                "efficiency_score": round(sg.efficiency.overall_score, 3),
                "efficiency_interpretation": sg.efficiency.interpretation,
                "red_flags": flag_summary,
            }
            yield f"event: guide\ndata: {json.dumps(guide_data)}\n\n"
            await asyncio.sleep(0.1)  # Delay to show progressive loading
        
        # Done event
        done_data = {"n_guides_returned": len(scored_guides), "n_guides_found": len(spacers)}
        yield f"event: done\ndata: {json.dumps(done_data)}\n\n"
        
    except TargetNotFoundError as e:
        yield f"event: error\ndata: {json.dumps({'error': str(e)})}\n\n"
    except Exception as e:
        yield f"event: error\ndata: {json.dumps({'error': str(e)})}\n\n"


@app.post("/api/design/stream")
async def design_guides_stream(request: DesignRequest):
    """Stream gRNA design results using Server-Sent Events.
    
    Events:
        - init: Target information
        - progress: Processing stage updates
        - guide: Individual guide (streamed one at a time)
        - done: Completion summary
        - error: Error information
    """
    return StreamingResponse(
        design_stream_generator(request),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",  # Disable nginx buffering
        }
    )


@app.post("/api/offtarget")
async def analyze_offtargets(request: OffTargetRequest):
    """Analyze off-target sites for a spacer."""
    try:
        from scalpel.offtarget import OffTargetSearcher, RiskCalculator
        from scalpel.models.enums import Genome
        
        genome_enum = Genome.from_string(request.genome)
        searcher = OffTargetSearcher(genome_enum, max_mismatches=request.max_mismatches)
        risk_calc = RiskCalculator(genome_enum)
        
        analysis = searcher.search(request.spacer)
        summary = risk_calc.generate_summary(analysis)
        
        sites = []
        for site in analysis.sites[:20]:
            sites.append({
                "chromosome": site.chromosome,
                "position": site.position,
                "sequence": site.sequence,
                "pam": site.pam,
                "mismatches": site.mismatch_count,
                "cfd_score": round(site.cutting_probability, 4),
            })
        
        return {
            "status": "success",
            "spacer": request.spacer,
            "total_sites": summary["total_sites"],
            "specificity_score": summary["specificity_score"],
            "interpretation": summary["interpretation"],
            "sites_by_mismatch": summary["sites_by_mismatch"],
            "sites": sites,
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/plan")
async def generate_plan(request: PlanRequest):
    """Generate experiment plan from design data."""
    try:
        from scalpel.planning import ExperimentPlanner
        from scalpel.models.enums import EditModality, Strand
        from scalpel.models.data_classes import SpacerCandidate, EfficiencyScore, DesignedGuide, GeneInfo
        
        data = request.design_data
        modality_str = data.get("target", {}).get("modality", "knockout")
        modality = EditModality(modality_str)
        
        planner = ExperimentPlanner(modality)
        
        # Convert guides
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
        
        # Gene info
        gene_info = None
        target = data.get("target", {})
        if target.get("gene"):
            gene_info = GeneInfo(
                gene_id="",
                symbol=target.get("gene", ""),
                chromosome=target.get("chromosome", ""),
                start=target.get("start", 0),
                end=target.get("end", 0),
                strand=Strand(target.get("strand", "+")),
            )
        
        plan = planner.generate_plan(guides, gene_info)
        
        return {
            "status": "success",
            "target_gene": plan.target_gene,
            "modality": plan.edit_modality.value,
            "positive_controls": [
                {"name": c.name, "guide": c.guide_sequence, "expected": c.expected_result}
                for c in plan.positive_controls
            ],
            "negative_controls": [
                {"name": c.name, "guide": c.guide_sequence}
                for c in plan.negative_controls
            ],
            "validation_assays": [
                {"tier": a.tier, "name": a.name, "timing": a.timing, "success": a.success_criteria}
                for a in plan.validation_assays
            ],
            "failure_modes": [
                {"failure": f.failure, "probability": f.probability, "troubleshooting": f.troubleshooting[:3]}
                for f in plan.failure_modes
            ],
            "provenance": plan.input_hash,
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def get_html_ui():
    """Return the HTML/CSS/JS for the web interface."""
    return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCALPEL - Computational CRISPR Design Platform</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.9.0/p5.min.js"></script>
    <style>
        :root {
            /* Primary - Scientific Blue */
            --sci-primary: #1a365d;
            --sci-primary-light: #2c5282;
            --sci-accent: #3182ce;

            /* Backgrounds - Clean Academic White */
            --sci-bg: #f8fafc;
            --sci-bg-alt: #f1f5f9;
            --sci-panel: #ffffff;
            --sci-sidebar: #1e3a5f;

            /* Status Colors */
            --sci-success: #059669;
            --sci-good: #0891b2;
            --sci-warning: #d97706;
            --sci-danger: #dc2626;

            /* DNA Nucleotide Colors */
            --dna-adenine: #22c55e;
            --dna-thymine: #ef4444;
            --dna-guanine: #facc15;
            --dna-cytosine: #3b82f6;
            --dna-pam: #a855f7;

            /* Text */
            --sci-text: #0f172a;
            --sci-text-secondary: #475569;
            --sci-text-muted: #94a3b8;

            /* Borders & Shadows */
            --sci-border: #e2e8f0;
            --sci-shadow: rgba(0,0,0,0.1);
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--sci-bg);
            color: var(--sci-text);
            min-height: 100vh;
        }

        .container {
            max-width: 1600px;
            margin: 0 auto;
            padding: 0;
        }

        header {
            background: var(--sci-sidebar);
            color: white;
            padding: 1.5rem 2rem;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .header-left {
            display: flex;
            align-items: center;
            gap: 1.5rem;
        }

        #header-canvas {
            border-radius: 8px;
        }

        .logo-text {
            display: flex;
            flex-direction: column;
        }

        .logo {
            font-size: 2rem;
            font-weight: 700;
            letter-spacing: 2px;
            color: white;
        }

        .tagline {
            color: rgba(255,255,255,0.7);
            font-size: 0.9rem;
            font-weight: 400;
        }

        .main-grid {
            display: grid;
            grid-template-columns: 320px 1fr;
            min-height: calc(100vh - 80px);
        }

        aside {
            background: var(--sci-panel);
            border-right: 1px solid var(--sci-border);
            padding: 1.5rem;
        }

        .panel {
            background: var(--sci-panel);
            border-radius: 8px;
            border: 1px solid var(--sci-border);
            padding: 1.5rem;
            box-shadow: 0 1px 3px var(--sci-shadow);
        }

        .panel-title {
            font-size: 0.75rem;
            font-weight: 600;
            color: var(--sci-text-muted);
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 1.25rem;
            padding-bottom: 0.75rem;
            border-bottom: 1px solid var(--sci-border);
        }

        .form-group {
            margin-bottom: 1rem;
        }

        label {
            display: block;
            font-size: 0.8rem;
            font-weight: 500;
            color: var(--sci-text-secondary);
            margin-bottom: 0.4rem;
        }

        input, select {
            width: 100%;
            padding: 0.6rem 0.75rem;
            background: var(--sci-bg);
            border: 1px solid var(--sci-border);
            border-radius: 6px;
            color: var(--sci-text);
            font-size: 0.9rem;
            transition: all 0.2s;
        }

        input:focus, select:focus {
            outline: none;
            border-color: var(--sci-accent);
            box-shadow: 0 0 0 3px rgba(49, 130, 206, 0.15);
        }

        input::placeholder {
            color: var(--sci-text-muted);
        }

        .btn-row {
            display: flex;
            gap: 0.5rem;
            margin-top: 1rem;
        }

        .btn {
            flex: 1;
            padding: 0.7rem 1rem;
            border: none;
            border-radius: 6px;
            font-size: 0.85rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s;
        }

        .btn-primary {
            background: var(--sci-accent);
            color: white;
        }

        .btn-primary:hover {
            background: var(--sci-primary-light);
        }

        .btn-primary:disabled {
            opacity: 0.5;
            cursor: not-allowed;
        }

        .btn-secondary {
            background: var(--sci-bg-alt);
            color: var(--sci-text-secondary);
            border: 1px solid var(--sci-border);
        }

        .btn-secondary:hover {
            background: var(--sci-border);
        }

        main {
            padding: 1.5rem 2rem;
            background: var(--sci-bg);
        }

        .tabs {
            display: flex;
            gap: 0;
            margin-bottom: 1.5rem;
            border-bottom: 2px solid var(--sci-border);
        }

        .tab {
            padding: 0.75rem 1.5rem;
            background: transparent;
            border: none;
            color: var(--sci-text-muted);
            cursor: pointer;
            transition: all 0.2s;
            font-weight: 500;
            position: relative;
        }

        .tab.active {
            color: var(--sci-accent);
        }

        .tab.active::after {
            content: '';
            position: absolute;
            bottom: -2px;
            left: 0;
            right: 0;
            height: 2px;
            background: var(--sci-accent);
        }

        .summary-stats {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1rem;
            margin-bottom: 1.5rem;
        }

        .stat-card {
            background: var(--sci-panel);
            padding: 1rem 1.25rem;
            border-radius: 8px;
            border: 1px solid var(--sci-border);
        }

        .stat-value {
            font-size: 1.5rem;
            font-weight: 700;
            color: var(--sci-primary);
        }

        .stat-label {
            font-size: 0.75rem;
            color: var(--sci-text-muted);
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-top: 0.25rem;
        }

        .visualization-row {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1rem;
            margin-bottom: 1.5rem;
        }

        .viz-panel {
            background: var(--sci-panel);
            border-radius: 8px;
            border: 1px solid var(--sci-border);
            padding: 1rem;
        }

        .viz-title {
            font-size: 0.75rem;
            font-weight: 600;
            color: var(--sci-text-muted);
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-bottom: 0.75rem;
        }

        .results-area {
            display: flex;
            flex-direction: column;
            gap: 1rem;
        }

        .guide-card {
            background: var(--sci-panel);
            border-radius: 8px;
            padding: 1rem 1.25rem;
            border: 1px solid var(--sci-border);
            transition: all 0.2s;
        }

        .guide-card:hover {
            border-color: var(--sci-accent);
            box-shadow: 0 2px 8px var(--sci-shadow);
        }

        .guide-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.75rem;
        }

        .guide-rank {
            background: var(--sci-primary);
            color: white;
            width: 28px;
            height: 28px;
            border-radius: 6px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: 600;
            font-size: 0.85rem;
        }

        .guide-sequence {
            font-family: 'SF Mono', 'Monaco', 'Consolas', monospace;
            font-size: 0.95rem;
            letter-spacing: 0.5px;
        }

        .nucleotide-A { color: var(--dna-adenine); }
        .nucleotide-T { color: var(--dna-thymine); }
        .nucleotide-G { color: var(--dna-guanine); }
        .nucleotide-C { color: var(--dna-cytosine); }

        .score-badge {
            padding: 0.3rem 0.6rem;
            border-radius: 4px;
            font-size: 0.75rem;
            font-weight: 600;
        }

        .score-excellent {
            background: rgba(5, 150, 105, 0.1);
            color: var(--sci-success);
        }

        .score-good {
            background: rgba(8, 145, 178, 0.1);
            color: var(--sci-good);
        }

        .score-moderate {
            background: rgba(217, 119, 6, 0.1);
            color: var(--sci-warning);
        }

        .guide-details {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1rem;
            padding-top: 0.75rem;
            border-top: 1px solid var(--sci-border);
        }

        .detail-item {
            text-align: center;
        }

        .detail-label {
            font-size: 0.65rem;
            color: var(--sci-text-muted);
            text-transform: uppercase;
            margin-bottom: 0.2rem;
        }

        .detail-value {
            font-size: 0.85rem;
            font-weight: 500;
            color: var(--sci-text);
        }

        .loading {
            display: none;
            text-align: center;
            padding: 3rem;
        }

        .loading.active {
            display: block;
        }

        .spinner {
            width: 40px;
            height: 40px;
            border: 3px solid var(--sci-border);
            border-top-color: var(--sci-accent);
            border-radius: 50%;
            animation: spin 0.8s linear infinite;
            margin: 0 auto 1rem;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(-5px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .status-message {
            padding: 0.75rem 1rem;
            border-radius: 6px;
            margin-bottom: 1rem;
            display: none;
            font-size: 0.85rem;
        }

        .status-message.error {
            display: block;
            background: rgba(220, 38, 38, 0.1);
            border: 1px solid rgba(220, 38, 38, 0.2);
            color: var(--sci-danger);
        }

        .status-message.success {
            display: block;
            background: rgba(5, 150, 105, 0.1);
            border: 1px solid rgba(5, 150, 105, 0.2);
            color: var(--sci-success);
        }

        .empty-state {
            text-align: center;
            padding: 3rem 2rem;
            color: var(--sci-text-muted);
        }

        .empty-state svg {
            width: 60px;
            height: 60px;
            margin-bottom: 1rem;
            opacity: 0.4;
        }

        canvas {
            display: block;
        }

        .red-flag {
            margin-top: 0.75rem;
            padding: 0.6rem;
            background: rgba(217, 119, 6, 0.08);
            border-radius: 6px;
            font-size: 0.8rem;
            color: var(--sci-warning);
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="header-left">
                <div id="header-canvas"></div>
                <div class="logo-text">
                    <div class="logo">SCALPEL</div>
                    <div class="tagline">Computational CRISPR Design Platform</div>
                </div>
            </div>
        </header>

        <div class="main-grid">
            <aside>
                <div class="panel-title">Design Parameters</div>

                <div class="form-group">
                    <label>Gene Symbol</label>
                    <input type="text" id="gene" placeholder="e.g., TP53, BRCA1, MYC">
                </div>

                <div class="form-group">
                    <label>Modality</label>
                    <select id="modality">
                        <option value="knockout">Knockout (NHEJ)</option>
                        <option value="interference">CRISPRi (Interference)</option>
                        <option value="activation">CRISPRa (Activation)</option>
                        <option value="base_edit_cbe">Base Edit (CBE)</option>
                        <option value="base_edit_abe">Base Edit (ABE)</option>
                        <option value="prime_edit">Prime Edit</option>
                    </select>
                </div>

                <div class="form-group">
                    <label>Genome</label>
                    <select id="genome">
                        <option value="GRCh38">Human (GRCh38)</option>
                        <option value="GRCh37">Human (GRCh37)</option>
                        <option value="GRCm39">Mouse (GRCm39)</option>
                    </select>
                </div>

                <div class="form-group">
                    <label>Cas Variant</label>
                    <select id="cas_variant">
                        <option value="SpCas9">SpCas9 (NGG)</option>
                        <option value="SpCas9-NG">SpCas9-NG (NG)</option>
                        <option value="SaCas9">SaCas9 (NNGRRT)</option>
                        <option value="Cas12a">Cas12a (TTTV)</option>
                    </select>
                </div>

                <div class="form-group">
                    <label>Number of Guides</label>
                    <input type="number" id="n_guides" value="10" min="1" max="50">
                </div>

                <div class="btn-row">
                    <button class="btn btn-secondary" onclick="loadSample()">Sample</button>
                    <button class="btn btn-primary" id="designBtn" onclick="designGuides()">Design</button>
                </div>
            </aside>

            <main>
                <div class="tabs">
                    <button class="tab active" onclick="switchTab('guides')">Guides</button>
                    <button class="tab" onclick="switchTab('offtarget')">Off-Targets</button>
                    <button class="tab" onclick="switchTab('plan')">Experiment Plan</button>
                </div>

                <div id="status" class="status-message"></div>

                <div id="loading" class="loading">
                    <div class="spinner"></div>
                    <p style="color: var(--sci-text-secondary);">Analyzing target sequence...</p>
                </div>

                <div id="summary" class="summary-stats" style="display: none;">
                    <div class="stat-card">
                        <div class="stat-value" id="stat-total">0</div>
                        <div class="stat-label">Candidates</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-efficiency">0%</div>
                        <div class="stat-label">Avg Efficiency</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-gene">-</div>
                        <div class="stat-label">Target</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-modality">-</div>
                        <div class="stat-label">Modality</div>
                    </div>
                </div>

                <div id="visualizations" class="visualization-row" style="display: none;">
                    <div class="viz-panel">
                        <div class="viz-title">Efficiency Distribution</div>
                        <div id="efficiency-chart"></div>
                    </div>
                    <div class="viz-panel">
                        <div class="viz-title">Genomic Position Map</div>
                        <div id="genome-map"></div>
                    </div>
                </div>

                <div id="sequence-viz" class="viz-panel" style="display: none; margin-bottom: 1.5rem;">
                    <div class="viz-title">Target Sequence</div>
                    <div id="sequence-canvas"></div>
                </div>

                <div id="results-guides" class="results-area">
                    <div class="empty-state">
                        <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="1.5" d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
                        </svg>
                        <p>Enter a gene symbol and click "Design" to analyze gRNA candidates</p>
                        <p style="margin-top: 0.5rem; font-size: 0.85rem;">Or click "Sample" to load TP53 demo data</p>
                    </div>
                </div>

                <div id="results-offtarget" class="results-area" style="display: none;"></div>
                <div id="results-plan" class="results-area" style="display: none;"></div>
            </main>
        </div>
    </div>
    
    <script>
        let currentDesignData = null;
        let headerSketch, efficiencySketch, genomeSketch, sequenceSketch;

        // ========== p5.js DNA Helix Animation ==========
        const headerSketchFn = (p) => {
            let angle = 0;
            const nucleotides = [];

            p.setup = () => {
                p.createCanvas(60, 60);
                p.noStroke();
                for (let i = 0; i < 8; i++) {
                    nucleotides.push({
                        y: i * 8,
                        base: ['A', 'T', 'G', 'C'][Math.floor(Math.random() * 4)]
                    });
                }
            };

            p.draw = () => {
                p.background(30, 58, 95);
                p.translate(30, 5);

                for (let i = 0; i < nucleotides.length; i++) {
                    const y = nucleotides[i].y;
                    const x1 = Math.sin(angle + i * 0.5) * 18;
                    const x2 = Math.sin(angle + i * 0.5 + Math.PI) * 18;

                    // Backbone
                    p.fill(200, 200, 200, 150);
                    p.ellipse(x1, y, 4, 4);
                    p.ellipse(x2, y, 4, 4);

                    // Base pair connection
                    p.stroke(100, 150, 200, 100);
                    p.strokeWeight(1);
                    p.line(x1, y, x2, y);
                    p.noStroke();

                    // Nucleotides with colors
                    const base = nucleotides[i].base;
                    if (base === 'A') p.fill(34, 197, 94);
                    else if (base === 'T') p.fill(239, 68, 68);
                    else if (base === 'G') p.fill(250, 204, 21);
                    else p.fill(59, 130, 246);

                    p.ellipse(x1, y, 6, 6);
                    p.ellipse(x2, y, 6, 6);
                }
                angle += 0.02;
            };
        };

        // ========== p5.js Efficiency Chart ==========
        const efficiencySketchFn = (p) => {
            let guides = [];
            let targetGuides = [];

            p.setup = () => {
                const container = document.getElementById('efficiency-chart');
                p.createCanvas(container.offsetWidth - 20, 150);
            };

            p.draw = () => {
                p.background(255);
                if (guides.length === 0) {
                    p.fill(148, 163, 184);
                    p.textAlign(p.CENTER);
                    p.textSize(12);
                    p.text('Efficiency data will appear here', p.width/2, p.height/2);
                    return;
                }

                // Animate towards target
                for (let i = 0; i < targetGuides.length; i++) {
                    if (i >= guides.length) guides.push(0);
                    guides[i] += (targetGuides[i] - guides[i]) * 0.1;
                }

                const barWidth = Math.min(40, (p.width - 60) / guides.length - 4);
                const maxHeight = p.height - 40;

                // Axes
                p.stroke(226, 232, 240);
                p.line(40, 10, 40, p.height - 20);
                p.line(40, p.height - 20, p.width - 10, p.height - 20);

                // Y-axis labels
                p.fill(148, 163, 184);
                p.textSize(10);
                p.textAlign(p.RIGHT);
                p.noStroke();
                p.text('100%', 35, 15);
                p.text('50%', 35, maxHeight/2 + 10);
                p.text('0%', 35, p.height - 18);

                // Bars
                for (let i = 0; i < guides.length; i++) {
                    const x = 50 + i * (barWidth + 4);
                    const h = guides[i] * maxHeight;
                    const y = p.height - 20 - h;

                    // Color based on score
                    if (guides[i] >= 0.8) p.fill(5, 150, 105);
                    else if (guides[i] >= 0.6) p.fill(8, 145, 178);
                    else p.fill(217, 119, 6);

                    p.noStroke();
                    p.rect(x, y, barWidth, h, 2);

                    // Label
                    p.fill(71, 85, 105);
                    p.textAlign(p.CENTER);
                    p.textSize(9);
                    p.text(i + 1, x + barWidth/2, p.height - 8);
                }
            };

            p.updateData = (newGuides) => {
                targetGuides = newGuides.map(g => g.efficiency_score);
            };

            p.windowResized = () => {
                const container = document.getElementById('efficiency-chart');
                if (container) p.resizeCanvas(container.offsetWidth - 20, 150);
            };
        };

        // ========== p5.js Genome Position Map ==========
        const genomeSketchFn = (p) => {
            let guides = [];
            let genomicRange = { start: 0, end: 1000000 };

            p.setup = () => {
                const container = document.getElementById('genome-map');
                p.createCanvas(container.offsetWidth - 20, 100);
            };

            p.draw = () => {
                p.background(255);

                if (guides.length === 0) {
                    p.fill(148, 163, 184);
                    p.textAlign(p.CENTER);
                    p.textSize(12);
                    p.text('Genomic positions will appear here', p.width/2, p.height/2);
                    return;
                }

                // Chromosome bar
                p.fill(241, 245, 249);
                p.stroke(226, 232, 240);
                p.rect(20, 40, p.width - 40, 20, 10);

                // Position labels
                p.fill(148, 163, 184);
                p.noStroke();
                p.textSize(9);
                p.textAlign(p.LEFT);
                p.text(genomicRange.start.toLocaleString(), 20, 75);
                p.textAlign(p.RIGHT);
                p.text(genomicRange.end.toLocaleString(), p.width - 20, 75);

                // Guide markers
                const range = genomicRange.end - genomicRange.start;
                guides.forEach((g, i) => {
                    const x = 20 + ((g.cut_site - genomicRange.start) / range) * (p.width - 40);

                    // Marker line
                    p.stroke(49, 130, 206);
                    p.strokeWeight(2);
                    p.line(x, 30, x, 70);

                    // Marker dot
                    p.noStroke();
                    if (g.strand === '+') p.fill(34, 197, 94);
                    else p.fill(239, 68, 68);
                    p.ellipse(x, 35, 8, 8);

                    // Rank label
                    p.fill(26, 54, 93);
                    p.textAlign(p.CENTER);
                    p.textSize(8);
                    p.text(g.rank, x, 25);
                });

                // Legend
                p.textSize(9);
                p.fill(34, 197, 94);
                p.ellipse(p.width - 80, 12, 6, 6);
                p.fill(71, 85, 105);
                p.textAlign(p.LEFT);
                p.text('+ strand', p.width - 72, 15);

                p.fill(239, 68, 68);
                p.ellipse(p.width - 80, 24, 6, 6);
                p.fill(71, 85, 105);
                p.text('- strand', p.width - 72, 27);
            };

            p.updateData = (newGuides, target) => {
                guides = newGuides;
                if (target) {
                    genomicRange = { start: target.start, end: target.end };
                }
            };

            p.windowResized = () => {
                const container = document.getElementById('genome-map');
                if (container) p.resizeCanvas(container.offsetWidth - 20, 100);
            };
        };

        // ========== p5.js Sequence Visualizer ==========
        const sequenceSketchFn = (p) => {
            let sequence = '';
            let pamStart = -1;

            p.setup = () => {
                const container = document.getElementById('sequence-canvas');
                p.createCanvas(container.offsetWidth - 20, 60);
                p.textFont('monospace');
            };

            p.draw = () => {
                p.background(255);

                if (!sequence) return;

                const charWidth = 14;
                const startX = 10;
                const y = 35;

                for (let i = 0; i < sequence.length; i++) {
                    const char = sequence[i];
                    const x = startX + i * charWidth;

                    // Background for PAM
                    if (pamStart >= 0 && i >= pamStart && i < pamStart + 3) {
                        p.fill(168, 85, 247, 40);
                        p.noStroke();
                        p.rect(x - 2, y - 14, charWidth, 20, 2);
                    }

                    // Nucleotide color
                    if (char === 'A') p.fill(34, 197, 94);
                    else if (char === 'T') p.fill(239, 68, 68);
                    else if (char === 'G') p.fill(250, 204, 21);
                    else if (char === 'C') p.fill(59, 130, 246);
                    else p.fill(100);

                    p.textSize(14);
                    p.textAlign(p.CENTER);
                    p.text(char, x + charWidth/2, y);

                    // Position markers
                    if (i % 5 === 0) {
                        p.fill(148, 163, 184);
                        p.textSize(8);
                        p.text(i + 1, x + charWidth/2, y + 18);
                    }
                }

                // Legend
                p.textSize(9);
                p.textAlign(p.LEFT);
                p.fill(168, 85, 247);
                p.text('PAM', p.width - 40, 15);
            };

            p.updateSequence = (seq, pam) => {
                sequence = seq;
                pamStart = seq.length - 3;  // PAM at end typically
            };

            p.windowResized = () => {
                const container = document.getElementById('sequence-canvas');
                if (container) p.resizeCanvas(container.offsetWidth - 20, 60);
            };
        };

        // Initialize p5.js sketches
        document.addEventListener('DOMContentLoaded', () => {
            headerSketch = new p5(headerSketchFn, 'header-canvas');
        });

        // ========== Load Sample Data ==========
        function loadSample() {
            document.getElementById('gene').value = 'TP53';
            document.getElementById('modality').value = 'knockout';
            document.getElementById('genome').value = 'GRCh38';
            document.getElementById('cas_variant').value = 'SpCas9';
            document.getElementById('n_guides').value = '10';
            showStatus('Sample loaded: TP53 knockout. Click "Design" to analyze.', 'success');
        }

        // ========== Tab Switching ==========
        function switchTab(tab) {
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.querySelector(`.tab:nth-child(${tab === 'guides' ? 1 : tab === 'offtarget' ? 2 : 3})`).classList.add('active');

            document.getElementById('results-guides').style.display = tab === 'guides' ? 'flex' : 'none';
            document.getElementById('results-offtarget').style.display = tab === 'offtarget' ? 'flex' : 'none';
            document.getElementById('results-plan').style.display = tab === 'plan' ? 'flex' : 'none';

            // Show/hide visualizations based on tab
            const vizRow = document.getElementById('visualizations');
            const seqViz = document.getElementById('sequence-viz');
            if (tab === 'guides' && currentDesignData && currentDesignData.guides.length > 0) {
                vizRow.style.display = 'grid';
                seqViz.style.display = 'block';
            }

            if (tab === 'offtarget' && currentDesignData && currentDesignData.guides.length > 0) {
                analyzeOffTarget(currentDesignData.guides[0].spacer_sequence);
            }
            if (tab === 'plan' && currentDesignData) {
                generatePlan();
            }
        }

        // ========== Design Guides ==========
        async function designGuides() {
            const gene = document.getElementById('gene').value.trim();
            if (!gene) {
                showStatus('Please enter a gene symbol', 'error');
                return;
            }

            document.getElementById('loading').classList.add('active');
            document.getElementById('designBtn').disabled = true;
            document.getElementById('summary').style.display = 'none';
            document.getElementById('visualizations').style.display = 'none';
            document.getElementById('sequence-viz').style.display = 'none';
            document.getElementById('results-guides').innerHTML = '';
            hideStatus();

            currentDesignData = { guides: [], target: {} };

            try {
                const response = await fetch('/api/design/stream', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        gene: gene,
                        modality: document.getElementById('modality').value,
                        genome: document.getElementById('genome').value,
                        cas_variant: document.getElementById('cas_variant').value,
                        n_guides: parseInt(document.getElementById('n_guides').value)
                    })
                });

                if (!response.ok) {
                    const errorData = await response.json();
                    throw new Error(errorData.detail || 'Design failed');
                }

                const reader = response.body.getReader();
                const decoder = new TextDecoder();
                let buffer = '';

                while (true) {
                    const { done, value } = await reader.read();
                    if (done) break;

                    buffer += decoder.decode(value, { stream: true });
                    const lines = buffer.split('\\n');
                    buffer = lines.pop();

                    for (const line of lines) {
                        if (line.startsWith('data: ')) {
                            const data = JSON.parse(line.substring(6));
                            handleStreamEvent(data);
                        }
                    }
                }

                // Initialize visualizations after data loaded
                initializeVisualizations();
                showStatus(`Found ${currentDesignData.n_guides_found || 0} candidates, showing top ${currentDesignData.guides.length}`, 'success');

            } catch (error) {
                showStatus(error.message, 'error');
            } finally {
                document.getElementById('loading').classList.remove('active');
                document.getElementById('designBtn').disabled = false;
            }
        }

        function initializeVisualizations() {
            if (currentDesignData.guides.length === 0) return;

            document.getElementById('visualizations').style.display = 'grid';
            document.getElementById('sequence-viz').style.display = 'block';

            // Initialize efficiency chart
            if (!efficiencySketch) {
                efficiencySketch = new p5(efficiencySketchFn, 'efficiency-chart');
            }
            setTimeout(() => {
                if (efficiencySketch && efficiencySketch.updateData) {
                    efficiencySketch.updateData(currentDesignData.guides);
                }
            }, 100);

            // Initialize genome map
            if (!genomeSketch) {
                genomeSketch = new p5(genomeSketchFn, 'genome-map');
            }
            setTimeout(() => {
                if (genomeSketch && genomeSketch.updateData) {
                    genomeSketch.updateData(currentDesignData.guides, currentDesignData.target);
                }
            }, 100);

            // Initialize sequence visualizer
            if (!sequenceSketch) {
                sequenceSketch = new p5(sequenceSketchFn, 'sequence-canvas');
            }
            setTimeout(() => {
                if (sequenceSketch && sequenceSketch.updateSequence && currentDesignData.guides[0]) {
                    const seq = currentDesignData.guides[0].spacer_sequence + currentDesignData.guides[0].pam_sequence;
                    sequenceSketch.updateSequence(seq);
                }
            }, 100);
        }

        function handleStreamEvent(data) {
            const container = document.getElementById('results-guides');

            if (data.type === 'init' && data.target) {
                currentDesignData.target = data.target;
                document.getElementById('summary').style.display = 'grid';
                document.getElementById('stat-gene').textContent = data.target.gene || '-';
                document.getElementById('stat-modality').textContent = data.target.modality;
                document.getElementById('stat-total').textContent = '...';
                document.getElementById('stat-efficiency').textContent = '...';
            } else if (data.stage === 'scoring') {
                document.getElementById('stat-total').textContent = data.total_spacers;
            } else if (data.spacer_sequence) {
                currentDesignData.guides.push(data);
                container.innerHTML += renderGuideCard(data);

                const avgEff = currentDesignData.guides.reduce((sum, g) => sum + g.efficiency_score, 0) / currentDesignData.guides.length;
                document.getElementById('stat-efficiency').textContent = (avgEff * 100).toFixed(0) + '%';

                // Update visualizations in real-time
                if (efficiencySketch && efficiencySketch.updateData) {
                    efficiencySketch.updateData(currentDesignData.guides);
                }
            } else if (data.n_guides_found !== undefined) {
                currentDesignData.n_guides_found = data.n_guides_found;
            } else if (data.error) {
                showStatus(data.error, 'error');
            }
        }

        function colorizeSequence(seq) {
            return seq.split('').map(c => {
                const cls = 'nucleotide-' + c;
                return `<span class="${cls}">${c}</span>`;
            }).join('');
        }

        function renderGuideCard(guide) {
            return `
                <div class="guide-card" style="animation: fadeIn 0.3s ease-out;">
                    <div class="guide-header">
                        <div style="display: flex; align-items: center; gap: 0.75rem;">
                            <div class="guide-rank">${guide.rank}</div>
                            <div class="guide-sequence">${colorizeSequence(guide.spacer_sequence)}</div>
                        </div>
                        <span class="score-badge ${getScoreClass(guide.efficiency_score)}">
                            ${(guide.efficiency_score * 100).toFixed(0)}%
                        </span>
                    </div>
                    <div class="guide-details">
                        <div class="detail-item">
                            <div class="detail-label">PAM</div>
                            <div class="detail-value">${guide.pam_sequence}</div>
                        </div>
                        <div class="detail-item">
                            <div class="detail-label">Strand</div>
                            <div class="detail-value">${guide.strand}</div>
                        </div>
                        <div class="detail-item">
                            <div class="detail-label">Cut Site</div>
                            <div class="detail-value">${guide.cut_site.toLocaleString()}</div>
                        </div>
                        <div class="detail-item">
                            <div class="detail-label">Position</div>
                            <div class="detail-value">${guide.genomic_start.toLocaleString()}</div>
                        </div>
                    </div>
                    ${guide.red_flags && guide.red_flags.total_flags > 0 ? `
                        <div class="red-flag">Warning: ${guide.red_flags.interpretation}</div>
                    ` : ''}
                </div>
            `;
        }

        async function analyzeOffTarget(spacer) {
            const container = document.getElementById('results-offtarget');
            container.innerHTML = '<div class="loading active"><div class="spinner"></div><p style="color: var(--sci-text-secondary);">Analyzing off-targets...</p></div>';

            try {
                const response = await fetch('/api/offtarget', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ spacer, genome: document.getElementById('genome').value })
                });

                const data = await response.json();

                container.innerHTML = `
                    <div class="panel">
                        <div class="panel-title">Off-Target Analysis</div>
                        <div style="font-family: monospace; color: var(--sci-accent); margin-bottom: 1rem; font-size: 0.9rem;">${spacer}</div>
                        <div class="summary-stats" style="margin-bottom: 1rem; grid-template-columns: repeat(2, 1fr);">
                            <div class="stat-card">
                                <div class="stat-value">${data.total_sites}</div>
                                <div class="stat-label">Total Sites</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value">${data.specificity_score}%</div>
                                <div class="stat-label">Specificity</div>
                            </div>
                        </div>
                        <p style="color: var(--sci-text-secondary); margin-bottom: 1rem; font-size: 0.9rem;">${data.interpretation}</p>
                        <table style="width: 100%; border-collapse: collapse; font-size: 0.85rem;">
                            <thead>
                                <tr style="border-bottom: 2px solid var(--sci-border);">
                                    <th style="text-align: left; padding: 0.5rem; color: var(--sci-text-muted); font-weight: 600;">Chr</th>
                                    <th style="text-align: left; padding: 0.5rem; color: var(--sci-text-muted); font-weight: 600;">Sequence</th>
                                    <th style="text-align: center; padding: 0.5rem; color: var(--sci-text-muted); font-weight: 600;">MM</th>
                                    <th style="text-align: right; padding: 0.5rem; color: var(--sci-text-muted); font-weight: 600;">CFD</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${data.sites.slice(0, 10).map(s => `
                                    <tr style="border-bottom: 1px solid var(--sci-border);">
                                        <td style="padding: 0.5rem; color: var(--sci-text);">${s.chromosome}</td>
                                        <td style="padding: 0.5rem; font-family: monospace;">${colorizeSequence(s.sequence)}</td>
                                        <td style="text-align: center; padding: 0.5rem; color: var(--sci-text);">${s.mismatches}</td>
                                        <td style="text-align: right; padding: 0.5rem; color: var(--sci-text);">${(s.cfd_score * 100).toFixed(1)}%</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                `;
            } catch (error) {
                container.innerHTML = `<div class="status-message error">${error.message}</div>`;
            }
        }

        async function generatePlan() {
            const container = document.getElementById('results-plan');
            container.innerHTML = '<div class="loading active"><div class="spinner"></div><p style="color: var(--sci-text-secondary);">Generating experiment plan...</p></div>';

            try {
                const response = await fetch('/api/plan', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ design_data: currentDesignData })
                });

                const data = await response.json();

                container.innerHTML = `
                    <div class="panel">
                        <div class="panel-title">Experiment Plan: ${data.target_gene}</div>
                        <p style="color: var(--sci-text-muted); margin-bottom: 1.5rem; font-size: 0.8rem;">Provenance Hash: ${data.provenance}</p>

                        <h4 style="color: var(--sci-success); margin: 0 0 0.75rem; font-size: 0.85rem; font-weight: 600;">POSITIVE CONTROLS</h4>
                        ${data.positive_controls.map(c => `
                            <div style="background: var(--sci-bg-alt); padding: 0.75rem; border-radius: 6px; margin-bottom: 0.5rem;">
                                <strong style="color: var(--sci-text);">${c.name}</strong>
                                <div style="font-family: monospace; color: var(--sci-accent); font-size: 0.85rem; margin: 0.25rem 0;">${c.guide}</div>
                                <div style="font-size: 0.8rem; color: var(--sci-text-muted);">${c.expected}</div>
                            </div>
                        `).join('')}

                        <h4 style="color: var(--sci-good); margin: 1.5rem 0 0.75rem; font-size: 0.85rem; font-weight: 600;">VALIDATION ASSAYS</h4>
                        ${data.validation_assays.map(a => `
                            <div style="background: var(--sci-bg-alt); padding: 0.75rem; border-radius: 6px; margin-bottom: 0.5rem;">
                                <strong style="color: var(--sci-text);">Tier ${a.tier}: ${a.name}</strong>
                                <div style="font-size: 0.8rem; color: var(--sci-text-muted);">Timing: ${a.timing}</div>
                            </div>
                        `).join('')}

                        <h4 style="color: var(--sci-warning); margin: 1.5rem 0 0.75rem; font-size: 0.85rem; font-weight: 600;">TROUBLESHOOTING</h4>
                        ${data.failure_modes.map(f => `
                            <div style="background: rgba(217, 119, 6, 0.08); padding: 0.75rem; border-radius: 6px; margin-bottom: 0.5rem;">
                                <strong style="color: var(--sci-text);">${f.failure}</strong>
                                <span style="color: var(--sci-text-muted); font-size: 0.8rem;"> (${f.probability})</span>
                                <ul style="margin: 0.5rem 0 0 1.25rem; font-size: 0.8rem; color: var(--sci-text-secondary);">
                                    ${f.troubleshooting.map(t => `<li>${t}</li>`).join('')}
                                </ul>
                            </div>
                        `).join('')}
                    </div>
                `;
            } catch (error) {
                container.innerHTML = `<div class="status-message error">${error.message}</div>`;
            }
        }

        function getScoreClass(score) {
            if (score >= 0.8) return 'score-excellent';
            if (score >= 0.6) return 'score-good';
            return 'score-moderate';
        }

        function showStatus(message, type) {
            const status = document.getElementById('status');
            status.textContent = message;
            status.className = 'status-message ' + type;
        }

        function hideStatus() {
            document.getElementById('status').className = 'status-message';
        }
    </script>
</body>
</html>"""


if __name__ == "__main__":
    print("Starting SCALPEL Web Server...")
    print("Open http://localhost:8000 in your browser")
    uvicorn.run(app, host="0.0.0.0", port=8000)
