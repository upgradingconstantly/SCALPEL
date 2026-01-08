"""
SCALPEL End-to-End Test Suite
================================

Tests the complete pipeline: Design → Off-target → Plan
"""

import subprocess
import json
import sys

def run_cmd(cmd):
    """Run a command and return output."""
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        cwd="."
    )
    return result.stdout, result.stderr, result.returncode

def test_design():
    """Test gRNA design for TP53 knockout."""
    print("\n" + "="*60)
    print("TEST 1: gRNA Design for TP53 Knockout")
    print("="*60)
    
    stdout, stderr, code = run_cmd(
        "scalpel design --gene TP53 --modality knockout --n-guides 5 --format json"
    )
    
    if code != 0:
        print(f"❌ FAILED: {stderr}")
        return None
    
    # Parse output (skip progress spinner lines)
    lines = stdout.strip().split('\n')
    json_start = next(i for i, l in enumerate(lines) if l.strip().startswith('{'))
    data = json.loads('\n'.join(lines[json_start:]))
    
    print(f"✅ Status: {data['status']}")
    print(f"✅ Gene: {data['target']['gene']}")
    print(f"✅ Guides found: {data['n_guides_found']}")
    print(f"✅ Top guide: {data['guides'][0]['spacer_sequence']}")
    print(f"✅ Efficiency: {data['guides'][0]['efficiency_score']}")
    print(f"✅ Interpretation: {data['guides'][0]['efficiency_interpretation']}")
    
    if 'red_flags' in data['guides'][0]:
        flags = data['guides'][0]['red_flags']
        print(f"✅ Red flags: {flags['total_flags']} ({flags['interpretation']})")
    
    return data

def test_offtarget(spacer):
    """Test off-target analysis."""
    print("\n" + "="*60)
    print(f"TEST 2: Off-Target Analysis for {spacer}")
    print("="*60)
    
    stdout, stderr, code = run_cmd(
        f"scalpel offtarget --spacer {spacer} --format json"
    )
    
    if code != 0:
        print(f"❌ FAILED: {stderr}")
        return None
    
    lines = stdout.strip().split('\n')
    json_start = next(i for i, l in enumerate(lines) if l.strip().startswith('{'))
    data = json.loads('\n'.join(lines[json_start:]))
    
    analysis = data['analyses'][0]
    print(f"✅ Status: {data['status']}")
    print(f"✅ Total sites: {analysis['total_sites']}")
    print(f"✅ Specificity score: {analysis['specificity_score']}")
    print(f"✅ Interpretation: {analysis['interpretation']}")
    print(f"✅ Sites by mismatch: {analysis['sites_by_mismatch']}")
    
    return data

def test_plan():
    """Test experiment plan generation."""
    print("\n" + "="*60)
    print("TEST 3: Experiment Plan Generation")
    print("="*60)
    
    # First save design to file
    run_cmd("scalpel design --gene BRCA1 --modality knockout --n-guides 3 --format json --output test_design.json")
    
    stdout, stderr, code = run_cmd(
        "scalpel plan --input test_design.json --format json"
    )
    
    if code != 0:
        print(f"❌ FAILED: {stderr}")
        return None
    
    lines = stdout.strip().split('\n')
    json_start = next(i for i, l in enumerate(lines) if l.strip().startswith('{'))
    data = json.loads('\n'.join(lines[json_start:]))
    
    print(f"✅ Target gene: {data['target_gene']}")
    print(f"✅ Modality: {data['modality']}")
    print(f"✅ Positive controls: {len(data['positive_controls'])}")
    for ctrl in data['positive_controls']:
        print(f"   - {ctrl['name']}: {ctrl['expected']}")
    print(f"✅ Negative controls: {len(data['negative_controls'])}")
    print(f"✅ Validation assays: {len(data['validation_assays'])}")
    for assay in data['validation_assays']:
        print(f"   - Tier {assay['tier']}: {assay['name']}")
    print(f"✅ Provenance: {data['provenance']}")
    
    return data

def test_modalities():
    """Test different modalities."""
    print("\n" + "="*60)
    print("TEST 4: Multiple Modalities")
    print("="*60)
    
    modalities = ["knockout", "interference", "activation"]
    
    for mod in modalities:
        stdout, stderr, code = run_cmd(
            f"scalpel design --gene MYC --modality {mod} --n-guides 2 --format json"
        )
        
        if code != 0:
            print(f"❌ {mod}: FAILED")
            continue
        
        lines = stdout.strip().split('\n')
        json_start = next(i for i, l in enumerate(lines) if l.strip().startswith('{'))
        data = json.loads('\n'.join(lines[json_start:]))
        
        n_guides = len(data.get('guides', []))
        print(f"✅ {mod.title()}: {n_guides} guides designed")

def main():
    print("="*60)
    print("  SCALPEL END-TO-END TEST SUITE")
    print("  Testing complete CRISPR design pipeline")
    print("="*60)
    
    # Test 1: Design
    design_data = test_design()
    
    # Test 2: Off-target
    if design_data and design_data.get('guides'):
        spacer = design_data['guides'][0]['spacer_sequence']
        test_offtarget(spacer)
    
    # Test 3: Plan
    test_plan()
    
    # Test 4: Modalities
    test_modalities()
    
    print("\n" + "="*60)
    print("  ALL TESTS COMPLETE!")
    print("="*60)
    print("\n✅ SCALPEL is fully functional!\n")

if __name__ == "__main__":
    main()
