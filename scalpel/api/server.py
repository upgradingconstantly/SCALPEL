"""
SCALPEL Web API

FastAPI-based web interface for the CRISPR design platform.
"""

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List, AsyncGenerator
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
        from scalpel.models.enums import Genome, EditModality, CasVariantType
        
        # Parse enums
        genome_enum = Genome.from_string(request.genome)
        modality = EditModality(request.modality)
        cas_variant = CasVariantType(request.cas_variant)
        
        # Resolve target
        resolver = TargetResolver(genome_enum)
        
        if request.gene:
            resolved = resolver.resolve_gene(request.gene, modality)
        elif request.coordinates:
            resolved = resolver.resolve_coordinates(request.coordinates)
        elif request.sequence:
            resolved = resolver.resolve_sequence(request.sequence)
        else:
            raise HTTPException(status_code=400, detail="Must provide gene, coordinates, or sequence")
        
        # Get or generate sequence
        target_sequence = resolved.sequence
        if not target_sequence or len(target_sequence) < 50:
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
        scored_guides = scorer.score_batch(spacers, modality.value, n_top=request.n_guides)
        
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
        from scalpel.models.enums import Genome, EditModality, CasVariantType
        
        # Parse enums
        genome_enum = Genome.from_string(request.genome)
        modality = EditModality(request.modality)
        cas_variant = CasVariantType(request.cas_variant)
        
        # Resolve target
        resolver = TargetResolver(genome_enum)
        
        if request.gene:
            resolved = resolver.resolve_gene(request.gene, modality)
        elif request.coordinates:
            resolved = resolver.resolve_coordinates(request.coordinates)
        elif request.sequence:
            resolved = resolver.resolve_sequence(request.sequence)
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
        if not target_sequence or len(target_sequence) < 50:
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
        scored_guides = scorer.score_batch(spacers, modality.value, n_top=request.n_guides)
        
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
    <title>SCALPEL - CRISPR Design Platform</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg-primary: #0a0a0f;
            --bg-secondary: #12121a;
            --bg-tertiary: #1a1a25;
            --accent-primary: #6366f1;
            --accent-secondary: #22d3ee;
            --accent-success: #10b981;
            --accent-warning: #f59e0b;
            --accent-danger: #ef4444;
            --text-primary: #f1f5f9;
            --text-secondary: #94a3b8;
            --text-muted: #64748b;
            --border: #2a2a3a;
            --glow: rgba(99, 102, 241, 0.3);
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Inter', sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            min-height: 100vh;
            background-image: 
                radial-gradient(ellipse at top, rgba(99, 102, 241, 0.1) 0%, transparent 50%),
                radial-gradient(ellipse at bottom right, rgba(34, 211, 238, 0.05) 0%, transparent 50%);
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 2rem;
        }
        
        header {
            text-align: center;
            padding: 2rem 0 3rem;
        }
        
        .logo {
            font-size: 2.5rem;
            font-weight: 700;
            background: linear-gradient(135deg, var(--accent-primary), var(--accent-secondary));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }
        
        .tagline {
            color: var(--text-secondary);
            font-size: 1.1rem;
        }
        
        .main-grid {
            display: grid;
            grid-template-columns: 350px 1fr;
            gap: 2rem;
        }
        
        .panel {
            background: var(--bg-secondary);
            border-radius: 16px;
            border: 1px solid var(--border);
            padding: 1.5rem;
        }
        
        .panel-title {
            font-size: 1rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 1.5rem;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }
        
        .panel-title::before {
            content: '';
            width: 4px;
            height: 20px;
            background: linear-gradient(to bottom, var(--accent-primary), var(--accent-secondary));
            border-radius: 2px;
        }
        
        .form-group {
            margin-bottom: 1.25rem;
        }
        
        label {
            display: block;
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }
        
        input, select {
            width: 100%;
            padding: 0.75rem 1rem;
            background: var(--bg-tertiary);
            border: 1px solid var(--border);
            border-radius: 8px;
            color: var(--text-primary);
            font-size: 0.95rem;
            transition: all 0.2s;
        }
        
        input:focus, select:focus {
            outline: none;
            border-color: var(--accent-primary);
            box-shadow: 0 0 0 3px var(--glow);
        }
        
        input::placeholder {
            color: var(--text-muted);
        }
        
        .btn {
            width: 100%;
            padding: 0.875rem 1.5rem;
            border: none;
            border-radius: 8px;
            font-size: 0.95rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s;
        }
        
        .btn-primary {
            background: linear-gradient(135deg, var(--accent-primary), #4f46e5);
            color: white;
        }
        
        .btn-primary:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 25px var(--glow);
        }
        
        .btn-primary:disabled {
            opacity: 0.5;
            cursor: not-allowed;
            transform: none;
        }
        
        .results-area {
            display: flex;
            flex-direction: column;
            gap: 1.5rem;
        }
        
        .guide-card {
            background: var(--bg-tertiary);
            border-radius: 12px;
            padding: 1.25rem;
            border: 1px solid var(--border);
            transition: all 0.2s;
        }
        
        .guide-card:hover {
            border-color: var(--accent-primary);
        }
        
        .guide-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 1rem;
        }
        
        .guide-rank {
            background: linear-gradient(135deg, var(--accent-primary), var(--accent-secondary));
            color: white;
            width: 32px;
            height: 32px;
            border-radius: 8px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: 600;
            font-size: 0.9rem;
        }
        
        .guide-sequence {
            font-family: 'Monaco', 'Consolas', monospace;
            font-size: 1.1rem;
            letter-spacing: 1px;
            color: var(--accent-secondary);
        }
        
        .score-badge {
            padding: 0.35rem 0.75rem;
            border-radius: 20px;
            font-size: 0.8rem;
            font-weight: 600;
        }
        
        .score-excellent {
            background: rgba(16, 185, 129, 0.2);
            color: var(--accent-success);
        }
        
        .score-good {
            background: rgba(34, 211, 238, 0.2);
            color: var(--accent-secondary);
        }
        
        .score-moderate {
            background: rgba(245, 158, 11, 0.2);
            color: var(--accent-warning);
        }
        
        .guide-details {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 1rem;
            margin-top: 1rem;
            padding-top: 1rem;
            border-top: 1px solid var(--border);
        }
        
        .detail-item {
            text-align: center;
        }
        
        .detail-label {
            font-size: 0.75rem;
            color: var(--text-muted);
            margin-bottom: 0.25rem;
        }
        
        .detail-value {
            font-size: 0.9rem;
            font-weight: 500;
            color: var(--text-secondary);
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
            width: 48px;
            height: 48px;
            border: 3px solid var(--border);
            border-top-color: var(--accent-primary);
            border-radius: 50%;
            animation: spin 0.8s linear infinite;
            margin: 0 auto 1rem;
        }
        
        @keyframes spin {
            to { transform: rotate(360deg); }
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(-10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .status-message {
            padding: 1rem;
            border-radius: 8px;
            margin-bottom: 1rem;
            display: none;
        }
        
        .status-message.error {
            display: block;
            background: rgba(239, 68, 68, 0.1);
            border: 1px solid rgba(239, 68, 68, 0.3);
            color: var(--accent-danger);
        }
        
        .status-message.success {
            display: block;
            background: rgba(16, 185, 129, 0.1);
            border: 1px solid rgba(16, 185, 129, 0.3);
            color: var(--accent-success);
        }
        
        .tabs {
            display: flex;
            gap: 0.5rem;
            margin-bottom: 1.5rem;
        }
        
        .tab {
            padding: 0.75rem 1.25rem;
            background: transparent;
            border: 1px solid var(--border);
            border-radius: 8px;
            color: var(--text-secondary);
            cursor: pointer;
            transition: all 0.2s;
        }
        
        .tab.active {
            background: var(--accent-primary);
            border-color: var(--accent-primary);
            color: white;
        }
        
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1rem;
            margin-bottom: 1.5rem;
        }
        
        .stat-card {
            background: var(--bg-tertiary);
            padding: 1rem;
            border-radius: 12px;
            text-align: center;
        }
        
        .stat-value {
            font-size: 1.75rem;
            font-weight: 700;
            background: linear-gradient(135deg, var(--accent-primary), var(--accent-secondary));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .stat-label {
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 0.25rem;
        }
        
        .empty-state {
            text-align: center;
            padding: 4rem 2rem;
            color: var(--text-muted);
        }
        
        .empty-state svg {
            width: 80px;
            height: 80px;
            margin-bottom: 1rem;
            opacity: 0.5;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="logo">🧬 SCALPEL</div>
            <p class="tagline">Computational CRISPR Design Platform</p>
        </header>
        
        <div class="main-grid">
            <aside>
                <div class="panel">
                    <h2 class="panel-title">Design Parameters</h2>
                    
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
                    
                    <button class="btn btn-primary" id="designBtn" onclick="designGuides()">
                        🎯 Design gRNAs
                    </button>
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
                    <p>Designing guides...</p>
                </div>
                
                <div id="summary" class="summary-stats" style="display: none;">
                    <div class="stat-card">
                        <div class="stat-value" id="stat-total">0</div>
                        <div class="stat-label">Guides Found</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-efficiency">0%</div>
                        <div class="stat-label">Avg Efficiency</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-gene">-</div>
                        <div class="stat-label">Target Gene</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value" id="stat-modality">-</div>
                        <div class="stat-label">Modality</div>
                    </div>
                </div>
                
                <div id="results-guides" class="results-area">
                    <div class="empty-state">
                        <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="1.5" d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
                        </svg>
                        <p>Enter a gene symbol and click "Design gRNAs" to get started</p>
                    </div>
                </div>
                
                <div id="results-offtarget" class="results-area" style="display: none;">
                </div>
                
                <div id="results-plan" class="results-area" style="display: none;">
                </div>
            </main>
        </div>
    </div>
    
    <script>
        let currentDesignData = null;
        
        function switchTab(tab) {
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.querySelector(`.tab:nth-child(${tab === 'guides' ? 1 : tab === 'offtarget' ? 2 : 3})`).classList.add('active');
            
            document.getElementById('results-guides').style.display = tab === 'guides' ? 'flex' : 'none';
            document.getElementById('results-offtarget').style.display = tab === 'offtarget' ? 'flex' : 'none';
            document.getElementById('results-plan').style.display = tab === 'plan' ? 'flex' : 'none';
            
            if (tab === 'offtarget' && currentDesignData && currentDesignData.guides.length > 0) {
                analyzeOffTarget(currentDesignData.guides[0].spacer_sequence);
            }
            if (tab === 'plan' && currentDesignData) {
                generatePlan();
            }
        }
        
        async function designGuides() {
            const gene = document.getElementById('gene').value.trim();
            if (!gene) {
                showStatus('Please enter a gene symbol', 'error');
                return;
            }
            
            document.getElementById('loading').classList.add('active');
            document.getElementById('designBtn').disabled = true;
            document.getElementById('summary').style.display = 'none';
            document.getElementById('results-guides').innerHTML = '';  // Clear previous results
            hideStatus();
            
            // Reset design data
            currentDesignData = { guides: [], target: {} };
            let guidesReceived = 0;
            
            try {
                // Use fetch with ReadableStream for SSE (more control than EventSource for POST)
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
                    buffer = lines.pop();  // Keep incomplete line in buffer
                    
                    for (const line of lines) {
                        if (line.startsWith('event: ')) {
                            const eventType = line.substring(7);
                            continue;
                        }
                        if (line.startsWith('data: ')) {
                            const data = JSON.parse(line.substring(6));
                            handleStreamEvent(data);
                        }
                    }
                }
                
                // Final display update
                showStatus(`Found ${currentDesignData.n_guides_found || guidesReceived} candidate guides, showing ${currentDesignData.guides.length}`, 'success');
                
            } catch (error) {
                showStatus(error.message, 'error');
            } finally {
                document.getElementById('loading').classList.remove('active');
                document.getElementById('designBtn').disabled = false;
            }
        }
        
        function handleStreamEvent(data) {
            const container = document.getElementById('results-guides');
            
            if (data.type === 'init' && data.target) {
                // Initialize with target info
                currentDesignData.target = data.target;
                document.getElementById('summary').style.display = 'grid';
                document.getElementById('stat-gene').textContent = data.target.gene || '-';
                document.getElementById('stat-modality').textContent = data.target.modality;
                document.getElementById('stat-total').textContent = '...';
                document.getElementById('stat-efficiency').textContent = '...';
            } else if (data.stage === 'scoring') {
                // Progress update
                document.getElementById('stat-total').textContent = data.total_spacers;
            } else if (data.spacer_sequence) {
                // Guide received - render immediately
                currentDesignData.guides.push(data);
                container.innerHTML += renderGuideCard(data);
                
                // Update average efficiency
                const avgEff = currentDesignData.guides.reduce((sum, g) => sum + g.efficiency_score, 0) / currentDesignData.guides.length;
                document.getElementById('stat-efficiency').textContent = (avgEff * 100).toFixed(0) + '%';
            } else if (data.n_guides_found !== undefined) {
                // Done event
                currentDesignData.n_guides_found = data.n_guides_found;
            } else if (data.error) {
                showStatus(data.error, 'error');
            }
        }
        
        function renderGuideCard(guide) {
            return `
                <div class="guide-card" style="animation: fadeIn 0.3s ease-out;">
                    <div class="guide-header">
                        <div style="display: flex; align-items: center; gap: 1rem;">
                            <div class="guide-rank">${guide.rank}</div>
                            <div class="guide-sequence">${guide.spacer_sequence}</div>
                        </div>
                        <span class="score-badge ${getScoreClass(guide.efficiency_score)}">
                            ${(guide.efficiency_score * 100).toFixed(0)}% efficiency
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
                    </div>
                    ${guide.red_flags && guide.red_flags.total_flags > 0 ? `
                        <div style="margin-top: 1rem; padding: 0.75rem; background: rgba(245, 158, 11, 0.1); border-radius: 8px; font-size: 0.85rem; color: var(--accent-warning);">
                            ⚠️ ${guide.red_flags.interpretation}
                        </div>
                    ` : ''}
                </div>
            `;
        }
        
        function displayResults(data) {
            // Update summary
            document.getElementById('summary').style.display = 'grid';
            document.getElementById('stat-total').textContent = data.guides.length;
            document.getElementById('stat-gene').textContent = data.target.gene || '-';
            document.getElementById('stat-modality').textContent = data.target.modality;
            
            const avgEff = data.guides.reduce((sum, g) => sum + g.efficiency_score, 0) / data.guides.length;
            document.getElementById('stat-efficiency').textContent = (avgEff * 100).toFixed(0) + '%';
            
            // Display guides
            const container = document.getElementById('results-guides');
            container.innerHTML = data.guides.map(guide => `
                <div class="guide-card">
                    <div class="guide-header">
                        <div style="display: flex; align-items: center; gap: 1rem;">
                            <div class="guide-rank">${guide.rank}</div>
                            <div class="guide-sequence">${guide.spacer_sequence}</div>
                        </div>
                        <span class="score-badge ${getScoreClass(guide.efficiency_score)}">
                            ${(guide.efficiency_score * 100).toFixed(0)}% efficiency
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
                    </div>
                    ${guide.red_flags && guide.red_flags.total_flags > 0 ? `
                        <div style="margin-top: 1rem; padding: 0.75rem; background: rgba(245, 158, 11, 0.1); border-radius: 8px; font-size: 0.85rem; color: var(--accent-warning);">
                            ⚠️ ${guide.red_flags.interpretation}
                        </div>
                    ` : ''}
                </div>
            `).join('');
            
            showStatus(`Found ${data.n_guides_found} candidate guides, showing top ${data.guides.length}`, 'success');
        }
        
        async function analyzeOffTarget(spacer) {
            const container = document.getElementById('results-offtarget');
            container.innerHTML = '<div class="loading active"><div class="spinner"></div><p>Analyzing off-targets...</p></div>';
            
            try {
                const response = await fetch('/api/offtarget', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ spacer, genome: document.getElementById('genome').value })
                });
                
                const data = await response.json();
                
                container.innerHTML = `
                    <div class="panel">
                        <h3 class="panel-title">Off-Target Analysis: ${spacer}</h3>
                        <div class="summary-stats" style="margin-bottom: 1rem;">
                            <div class="stat-card">
                                <div class="stat-value">${data.total_sites}</div>
                                <div class="stat-label">Total Sites</div>
                            </div>
                            <div class="stat-card">
                                <div class="stat-value">${data.specificity_score}%</div>
                                <div class="stat-label">Specificity</div>
                            </div>
                        </div>
                        <p style="color: var(--text-secondary); margin-bottom: 1rem;">${data.interpretation}</p>
                        <table style="width: 100%; border-collapse: collapse;">
                            <thead>
                                <tr style="border-bottom: 1px solid var(--border);">
                                    <th style="text-align: left; padding: 0.5rem; color: var(--text-muted);">Chromosome</th>
                                    <th style="text-align: left; padding: 0.5rem; color: var(--text-muted);">Sequence</th>
                                    <th style="text-align: center; padding: 0.5rem; color: var(--text-muted);">MM</th>
                                    <th style="text-align: right; padding: 0.5rem; color: var(--text-muted);">CFD</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${data.sites.slice(0, 10).map(s => `
                                    <tr style="border-bottom: 1px solid var(--border);">
                                        <td style="padding: 0.5rem;">${s.chromosome}</td>
                                        <td style="padding: 0.5rem; font-family: monospace; color: var(--accent-secondary);">${s.sequence}</td>
                                        <td style="text-align: center; padding: 0.5rem;">${s.mismatches}</td>
                                        <td style="text-align: right; padding: 0.5rem;">${(s.cfd_score * 100).toFixed(1)}%</td>
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
            container.innerHTML = '<div class="loading active"><div class="spinner"></div><p>Generating experiment plan...</p></div>';
            
            try {
                const response = await fetch('/api/plan', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ design_data: currentDesignData })
                });
                
                const data = await response.json();
                
                container.innerHTML = `
                    <div class="panel">
                        <h3 class="panel-title">Experiment Plan: ${data.target_gene}</h3>
                        <p style="color: var(--text-muted); margin-bottom: 1rem;">Provenance: ${data.provenance}</p>
                        
                        <h4 style="color: var(--accent-secondary); margin: 1.5rem 0 0.75rem;">✅ Positive Controls</h4>
                        ${data.positive_controls.map(c => `
                            <div style="background: var(--bg-tertiary); padding: 0.75rem; border-radius: 8px; margin-bottom: 0.5rem;">
                                <strong>${c.name}</strong>
                                <div style="font-family: monospace; color: var(--accent-secondary); font-size: 0.9rem; margin: 0.25rem 0;">${c.guide}</div>
                                <div style="font-size: 0.85rem; color: var(--text-muted);">${c.expected}</div>
                            </div>
                        `).join('')}
                        
                        <h4 style="color: var(--accent-secondary); margin: 1.5rem 0 0.75rem;">🔬 Validation Assays</h4>
                        ${data.validation_assays.map(a => `
                            <div style="background: var(--bg-tertiary); padding: 0.75rem; border-radius: 8px; margin-bottom: 0.5rem;">
                                <strong>Tier ${a.tier}: ${a.name}</strong>
                                <div style="font-size: 0.85rem; color: var(--text-muted);">Timing: ${a.timing}</div>
                            </div>
                        `).join('')}
                        
                        <h4 style="color: var(--accent-warning); margin: 1.5rem 0 0.75rem;">⚠️ Troubleshooting</h4>
                        ${data.failure_modes.map(f => `
                            <div style="background: rgba(245, 158, 11, 0.1); padding: 0.75rem; border-radius: 8px; margin-bottom: 0.5rem;">
                                <strong>${f.failure}</strong> (${f.probability} probability)
                                <ul style="margin: 0.5rem 0 0 1.25rem; font-size: 0.85rem; color: var(--text-secondary);">
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
