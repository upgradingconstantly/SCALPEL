// SCALPEL Professional UI - Phase 2 JavaScript
// ==============================================

// ===== STATE =====
let currentDesignData = { guides: [], target: {}, allGuides: [] };
let headerSketch, efficiencySketch, genomeSketch, sequenceSketch;
let currentPage = 1;
let perPage = 50;
let infiniteScrollEnabled = false;
let isLoadingMore = false;
let selectedGuides = new Set();
let searchQuery = '';
let strandFilter = 'all';
let efficiencyMin = 0;
let efficiencyMax = 100;
let recentSearches = [];
const RECENT_SEARCHES_MAX = 10;
const INFINITE_SCROLL_BATCH = 75;

// ===== INIT =====
document.addEventListener('DOMContentLoaded', () => {
    headerSketch = new p5(headerSketchFn, 'header-canvas');
    loadRecentSearches();
    loadColorPreset();
    setupKeyboardShortcuts();
    setupInfiniteScroll();
    setupGeneInputDropdown();
    loadAvailableOptions();  // Load dynamic dropdowns
});

// ===== DYNAMIC OPTIONS =====
async function loadAvailableOptions() {
    try {
        const response = await fetch('/api/info');
        const data = await response.json();

        // Only show built genomes (7 currently available)
        const builtGenomes = ['GRCh38', 'GRCh37', 'GRCm39', 'BDGP6', 'WBcel235', 'Sscrofa11.1', 'TAIR10'];
        const availableGenomes = data.genomes.filter(g => builtGenomes.includes(g.value));

        // Populate genome dropdown
        const genomeSelect = document.getElementById('genome');
        genomeSelect.innerHTML = availableGenomes.map(g =>
            `<option value="${g.value}">${g.label}</option>`
        ).join('');

        // Populate modality dropdown
        const modalitySelect = document.getElementById('modality');
        modalitySelect.innerHTML = data.modalities.map(m =>
            `<option value="${m.value}">${m.label}</option>`
        ).join('');

        // Populate Cas variant dropdown
        const casSelect = document.getElementById('cas_variant');
        casSelect.innerHTML = data.cas_variants.map(c =>
            `<option value="${c.value}">${c.value} (${c.pam})</option>`
        ).join('');
    } catch (err) {
        console.error('Failed to load options:', err);
    }
}

// ===== RECENT SEARCHES =====
function loadRecentSearches() {
    const saved = localStorage.getItem('scalpel-recent-searches');
    if (saved) recentSearches = JSON.parse(saved);
}

function saveRecentSearch(gene) {
    gene = gene.toUpperCase().trim();
    if (!gene) return;
    recentSearches = recentSearches.filter(g => g !== gene);
    recentSearches.unshift(gene);
    if (recentSearches.length > RECENT_SEARCHES_MAX) {
        recentSearches = recentSearches.slice(0, RECENT_SEARCHES_MAX);
    }
    localStorage.setItem('scalpel-recent-searches', JSON.stringify(recentSearches));
}

function setupGeneInputDropdown() {
    const input = document.getElementById('gene');
    const dropdown = document.getElementById('recentDropdown');

    input.addEventListener('focus', () => {
        if (recentSearches.length > 0) {
            renderRecentDropdown();
            dropdown.classList.add('visible');
        }
    });

    input.addEventListener('blur', () => {
        setTimeout(() => dropdown.classList.remove('visible'), 150);
    });
}

function renderRecentDropdown() {
    const dropdown = document.getElementById('recentDropdown');
    if (recentSearches.length === 0) {
        dropdown.innerHTML = '';
        return;
    }
    dropdown.innerHTML = `
        <div class="recent-dropdown-header">Recently Searched</div>
        ${recentSearches.map(g => `<div class="recent-item" onclick="selectRecentSearch('${g}')">${g}</div>`).join('')}
    `;
}

function selectRecentSearch(gene) {
    document.getElementById('gene').value = gene;
    document.getElementById('recentDropdown').classList.remove('visible');
}

// ===== KEYBOARD SHORTCUTS =====
function setupKeyboardShortcuts() {
    document.addEventListener('keydown', (e) => {
        // Escape - close modals
        if (e.key === 'Escape') {
            closeAllModals();
        }

        // Enter in gene field - start design
        if (e.key === 'Enter' && e.target.id === 'gene') {
            e.preventDefault();
            designGuides();
        }

        // Arrow keys for pagination (when not in input)
        if (!['INPUT', 'TEXTAREA', 'SELECT'].includes(e.target.tagName)) {
            if (e.key === 'ArrowLeft') {
                prevPage();
            } else if (e.key === 'ArrowRight') {
                nextPage();
            }
        }

        // Ctrl+C with selected guides
        if (e.ctrlKey && e.key === 'c' && selectedGuides.size > 0) {
            e.preventDefault();
            bulkCopySelected();
        }
    });
}

function closeAllModals() {
    document.querySelectorAll('.modal-overlay').forEach(m => m.classList.remove('active'));
}

// ===== SETTINGS =====
function openSettings() {
    document.getElementById('settingsModal').classList.add('active');
}

function closeSettings() {
    document.getElementById('settingsModal').classList.remove('active');
}

function setColorPreset(preset) {
    const presets = {
        default: { A: '#22c55e', T: '#ef4444', G: '#eab308', C: '#3b82f6' },
        ucsc: { A: '#00aa00', T: '#ff0000', G: '#bbbb00', C: '#0000ff' },
        colorblind: { A: '#0072B2', T: '#D55E00', G: '#F0E442', C: '#009E73' }
    };
    const colors = presets[preset] || presets.default;

    document.documentElement.style.setProperty('--dna-adenine', colors.A);
    document.documentElement.style.setProperty('--dna-thymine', colors.T);
    document.documentElement.style.setProperty('--dna-guanine', colors.G);
    document.documentElement.style.setProperty('--dna-cytosine', colors.C);

    document.querySelectorAll('.color-preset').forEach(b => b.classList.remove('active'));
    document.querySelector(`[data-preset="${preset}"]`)?.classList.add('active');
    localStorage.setItem('scalpel-colors', preset);

    if (currentDesignData.guides.length > 0) renderGuides();
}

function loadColorPreset() {
    const saved = localStorage.getItem('scalpel-colors');
    if (saved) setColorPreset(saved);
}

// ===== P5.js SKETCHES =====
const headerSketchFn = (p) => {
    let angle = 0;
    const nucleotides = [];

    p.setup = () => {
        p.createCanvas(50, 50);
        p.noStroke();
        for (let i = 0; i < 7; i++) {
            nucleotides.push({ y: i * 7.5, base: ['A', 'T', 'G', 'C'][Math.floor(Math.random() * 4)] });
        }
    };

    p.draw = () => {
        p.background(30, 58, 95);
        p.translate(25, 3);
        for (let i = 0; i < nucleotides.length; i++) {
            const y = nucleotides[i].y;
            const x1 = Math.sin(angle + i * 0.5) * 14;
            const x2 = Math.sin(angle + i * 0.5 + Math.PI) * 14;
            p.fill(200, 200, 200, 150);
            p.ellipse(x1, y, 3, 3);
            p.ellipse(x2, y, 3, 3);
            p.stroke(100, 150, 200, 100);
            p.strokeWeight(0.5);
            p.line(x1, y, x2, y);
            p.noStroke();
            const base = nucleotides[i].base;
            if (base === 'A') p.fill(34, 197, 94);
            else if (base === 'T') p.fill(239, 68, 68);
            else if (base === 'G') p.fill(234, 179, 8);
            else p.fill(59, 130, 246);
            p.ellipse(x1, y, 5, 5);
            p.ellipse(x2, y, 5, 5);
        }
        angle += 0.02;
    };
};

const efficiencySketchFn = (p) => {
    let bins = [];

    p.setup = () => {
        const c = document.getElementById('efficiency-chart');
        p.createCanvas(c.offsetWidth - 16, 180);
    };

    p.draw = () => {
        p.background(255);
        if (bins.length === 0) {
            p.fill(148, 163, 184);
            p.textAlign(p.CENTER);
            p.textSize(12);
            p.noStroke();
            p.text('Efficiency distribution will appear here', p.width / 2, p.height / 2);
            return;
        }

        const maxCount = Math.max(...bins, 1);
        const barW = (p.width - 60) / bins.length;
        const maxH = p.height - 50;

        // Axes
        p.stroke(226, 232, 240);
        p.line(40, 15, 40, p.height - 30);
        p.line(40, p.height - 30, p.width - 15, p.height - 30);

        // Bars
        p.noStroke();
        for (let i = 0; i < bins.length; i++) {
            const h = (bins[i] / maxCount) * maxH;
            const x = 45 + i * barW;
            const pct = (i + 0.5) * 10;
            if (pct >= 80) p.fill(5, 150, 105);
            else if (pct >= 60) p.fill(8, 145, 178);
            else p.fill(217, 119, 6);
            p.rect(x, p.height - 30 - h, barW - 3, h, 2);
        }

        // Labels
        p.fill(100, 116, 139);
        p.textSize(10);
        p.textAlign(p.CENTER);
        p.noStroke();
        for (let i = 0; i <= 10; i += 2) {
            const x = 45 + i * barW;
            p.text(i * 10 + '%', x, p.height - 12);
        }
    };

    p.updateData = (guides) => {
        bins = new Array(10).fill(0);
        guides.forEach(g => {
            const idx = Math.min(9, Math.floor((g.efficiency_score || 0) * 10));
            bins[idx]++;
        });
    };

    p.windowResized = () => {
        const c = document.getElementById('efficiency-chart');
        if (c) p.resizeCanvas(c.offsetWidth - 16, 180);
    };
};

const genomeSketchFn = (p) => {
    let guides = [];
    let range = { start: 0, end: 1000000 };
    let isFullscreen = false;

    p.setup = () => {
        const c = document.getElementById('genome-map');
        p.createCanvas(c.offsetWidth - 16, 140);
    };

    p.draw = () => {
        p.background(255);
        if (guides.length === 0) {
            p.fill(148, 163, 184);
            p.textAlign(p.CENTER);
            p.textSize(12);
            p.noStroke();
            p.text('Guide positions will appear here', p.width / 2, p.height / 2);
            return;
        }

        const yMid = isFullscreen ? 200 : 70;
        const barH = isFullscreen ? 28 : 20;

        // Chromosome bar
        p.fill(241, 245, 249);
        p.stroke(200, 210, 220);
        p.strokeWeight(1);
        p.rect(40, yMid - barH / 2, p.width - 80, barH, barH / 2);

        // Position labels - BELOW the bar with padding
        p.fill(100, 116, 139);
        p.noStroke();
        p.textSize(10);
        p.textAlign(p.LEFT);
        p.text(range.start.toLocaleString(), 40, yMid + barH / 2 + 22);
        p.textAlign(p.RIGHT);
        p.text(range.end.toLocaleString(), p.width - 40, yMid + barH / 2 + 22);

        // Compute positions and stagger to avoid overlap
        const rng = range.end - range.start;
        const positions = guides.map((g, i) => ({
            x: 40 + ((g.cut_site - range.start) / rng) * (p.width - 80),
            g,
            index: i
        }));

        // Sort by x position
        positions.sort((a, b) => a.x - b.x);

        // Stagger dots and labels to avoid overlap
        let lastDotX = -50;
        let lastLabelX = -100;
        let alternateUp = false;

        positions.forEach((pos, idx) => {
            let { x, g } = pos;

            // Multi-level staggering for labels when too close
            let dotY = yMid - barH - 10;
            let labelRow = 0; // 0 = top, 1 = above bar, 2 = below bar

            if (x - lastDotX < 15) {
                alternateUp = !alternateUp;
                dotY = alternateUp ? yMid - barH - 22 : yMid - barH - 6;
            } else {
                alternateUp = false;
            }

            // Use bottom row for labels if too close
            if (x - lastLabelX < 30) {
                labelRow = (labelRow + 1) % 3;
            }
            lastDotX = x;

            // Draw guide line - TOUCH THE BOTTOM of the bar completely
            p.stroke(49, 130, 206);
            p.strokeWeight(2);
            p.line(x, dotY + 4, x, yMid + barH / 2); // Line to exact bottom of bar

            // Draw strand dot
            p.noStroke();
            p.fill(g.strand === '+' ? [34, 197, 94] : [239, 68, 68]);
            p.ellipse(x, dotY, 10, 10);

            // Draw rank label - multi-row stagger
            let labelY;
            if (labelRow === 0) {
                labelY = dotY - 12;
            } else if (labelRow === 1) {
                labelY = dotY - 24;
            } else {
                labelY = yMid + barH / 2 + 35; // Below bar
            }

            if (x - lastLabelX < 25 && labelRow === 0) {
                labelY = dotY - 24;
            }
            lastLabelX = x;

            p.fill(26, 54, 93);
            p.textAlign(p.CENTER);
            p.textSize(9);
            p.text(g.rank, x, labelY);
        });

        // Legend
        p.textSize(10);
        p.fill(34, 197, 94);
        p.noStroke();
        p.ellipse(p.width - 75, 15, 8);
        p.fill(71, 85, 105);
        p.textAlign(p.LEFT);
        p.text('+', p.width - 67, 19);
        p.fill(239, 68, 68);
        p.ellipse(p.width - 45, 15, 8);
        p.fill(71, 85, 105);
        p.text('-', p.width - 37, 19);
    };

    p.updateData = (g, t) => {
        guides = g.slice(0, isFullscreen ? 50 : 25);
        if (t && t.start && t.end) range = { start: t.start, end: t.end };
    };

    p.setFullscreen = (fs) => {
        isFullscreen = fs;
        const containerId = isFullscreen ? 'genome-map-fullscreen' : 'genome-map';
        const c = document.getElementById(containerId);
        if (c) p.resizeCanvas(c.offsetWidth - 20, isFullscreen ? 420 : 140);
    };

    p.windowResized = () => {
        if (!isFullscreen) {
            const c = document.getElementById('genome-map');
            if (c) p.resizeCanvas(c.offsetWidth - 16, 140);
        }
    };
};

const sequenceSketchFn = (p) => {
    let sequence = '';
    let pamStart = -1;

    p.setup = () => {
        const c = document.getElementById('sequence-canvas');
        p.createCanvas(c.offsetWidth - 16, 65);
        p.textFont('monospace');
    };

    p.draw = () => {
        p.background(250, 250, 248); // Warm neutral background for DNA colors
        if (!sequence) return;

        const charW = 15;
        const y = 35;

        for (let i = 0; i < sequence.length; i++) {
            const ch = sequence[i];
            const x = 10 + i * charW;

            // Highlight PAM
            if (pamStart >= 0 && i >= pamStart && i < pamStart + 3) {
                p.fill(168, 85, 247, 40);
                p.noStroke();
                p.rect(x - 1, y - 14, charW, 22, 2);
            }

            // Get colors from CSS
            const style = getComputedStyle(document.documentElement);
            if (ch === 'A') p.fill(style.getPropertyValue('--dna-adenine') || '#22c55e');
            else if (ch === 'T') p.fill(style.getPropertyValue('--dna-thymine') || '#ef4444');
            else if (ch === 'G') p.fill(style.getPropertyValue('--dna-guanine') || '#eab308');
            else if (ch === 'C') p.fill(style.getPropertyValue('--dna-cytosine') || '#3b82f6');
            else p.fill(100);

            p.textSize(15);
            p.textAlign(p.CENTER);
            p.noStroke();
            p.text(ch, x + charW / 2, y);

            if (i % 5 === 0) {
                p.fill(100, 116, 139);
                p.textSize(8);
                p.text(i + 1, x + charW / 2, y + 18);
            }
        }

        // PAM label - larger
        p.textSize(12);
        p.textAlign(p.LEFT);
        p.fill(168, 85, 247);
        p.noStroke();
        p.text('PAM', p.width - 35, 15);
    };

    p.updateSequence = (seq) => {
        sequence = seq;
        pamStart = seq.length - 3;
    };

    p.windowResized = () => {
        const c = document.getElementById('sequence-canvas');
        if (c) p.resizeCanvas(c.offsetWidth - 16, 65);
    };
};

// ===== TAB SWITCHING =====
function switchTab(tab) {
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    event.target.classList.add('active');

    document.getElementById('results-guides').style.display = tab === 'guides' ? 'block' : 'none';
    document.getElementById('results-offtarget').style.display = tab === 'offtarget' ? 'block' : 'none';
    document.getElementById('results-plan').style.display = tab === 'plan' ? 'block' : 'none';
    document.getElementById('results-visualization').style.display = tab === 'visualization' ? 'block' : 'none';

    // Show/hide visualizations based on tab
    if (tab === 'guides' && currentDesignData.guides?.length > 0) {
        document.getElementById('visualizations').style.display = 'flex';
        document.getElementById('sequence-viz').style.display = 'block';
    } else {
        document.getElementById('visualizations').style.display = 'none';
        document.getElementById('sequence-viz').style.display = 'none';
    }

    if (tab === 'offtarget' && currentDesignData.guides?.length > 0) {
        analyzeOffTarget(currentDesignData.guides[0].spacer_sequence);
    }
    if (tab === 'plan' && currentDesignData.guides?.length > 0) {
        generatePlan();
    }
    if (tab === 'visualization' && currentDesignData.guides?.length > 0) {
        initDNAHelixVisualization();
    }
}

// ===== 3D DNA HELIX VISUALIZATION (Three.js) =====
let dnaScene, dnaCamera, dnaRenderer, dnaControls;
let dnaHelixGroup;
let helixRotationSpeed = 0.003;
let dnaAnimationId;
let dnaTooltip;

// Create text sprite for nucleotide labels
function createTextSprite(text, color = '#ffffff') {
    const canvas = document.createElement('canvas');
    canvas.width = 64;
    canvas.height = 64;
    const ctx = canvas.getContext('2d');

    ctx.fillStyle = color;
    ctx.font = 'bold 48px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, 32, 32);

    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;

    const material = new THREE.SpriteMaterial({
        map: texture,
        transparent: true,
        depthTest: false
    });
    const sprite = new THREE.Sprite(material);
    sprite.scale.set(0.6, 0.6, 1);
    return sprite;
}

function initDNAHelixVisualization() {
    const container = document.getElementById('dna-helix-viz');
    if (!container) return;

    // Clear previous
    container.innerHTML = '';
    if (dnaAnimationId) cancelAnimationFrame(dnaAnimationId);

    // Create tooltip element
    if (!dnaTooltip) {
        dnaTooltip = document.createElement('div');
        dnaTooltip.id = 'dna-tooltip';
        dnaTooltip.style.cssText = `
            position: absolute;
            background: rgba(26, 54, 93, 0.95);
            color: white;
            padding: 10px 14px;
            border-radius: 6px;
            font-size: 13px;
            font-family: sans-serif;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.2s;
            z-index: 1000;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
            line-height: 1.5;
        `;
        document.body.appendChild(dnaTooltip);
    }

    // Get sequence from current guide
    const sequence = currentDesignData.guides[0]?.spacer_sequence || 'GCAAGAGGTACAGCAAGGCC';
    const complementMap = { A: 'T', T: 'A', G: 'C', C: 'G' };
    const baseColors = {
        A: 0x22c55e, T: 0xef4444, G: 0xeab308, C: 0x3b82f6
    };

    // Scene setup
    dnaScene = new THREE.Scene();
    dnaScene.background = new THREE.Color(0xffffff);

    // Camera
    const width = container.offsetWidth;
    const height = container.offsetHeight || 400;
    dnaCamera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    dnaCamera.position.set(0, 5, 25);

    // Renderer
    dnaRenderer = new THREE.WebGLRenderer({ antialias: true });
    dnaRenderer.setSize(width, height);
    dnaRenderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    dnaRenderer.shadowMap.enabled = true;
    container.appendChild(dnaRenderer.domElement);

    // Orbit Controls for interactivity
    dnaControls = new THREE.OrbitControls(dnaCamera, dnaRenderer.domElement);
    dnaControls.enableDamping = true;
    dnaControls.dampingFactor = 0.05;
    dnaControls.minDistance = 10;
    dnaControls.maxDistance = 50;

    // Lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
    dnaScene.add(ambientLight);

    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(10, 20, 10);
    directionalLight.castShadow = true;
    dnaScene.add(directionalLight);

    const backLight = new THREE.DirectionalLight(0xffffff, 0.3);
    backLight.position.set(-10, -10, -10);
    dnaScene.add(backLight);

    // DNA Helix Group
    dnaHelixGroup = new THREE.Group();

    // Create DNA helix
    const helixRadius = 2;
    const helixPitch = 0.8;
    const basePairs = sequence.length;

    for (let i = 0; i < basePairs; i++) {
        const base = sequence[i];
        const complement = complementMap[base] || 'N';
        const angle = i * 0.6;
        const y = (i - basePairs / 2) * helixPitch;

        // Strand 1 position (5' strand)
        const x1 = Math.cos(angle) * helixRadius;
        const z1 = Math.sin(angle) * helixRadius;

        // Strand 2 position (3' strand - opposite)
        const x2 = Math.cos(angle + Math.PI) * helixRadius;
        const z2 = Math.sin(angle + Math.PI) * helixRadius;

        // Nucleotide sphere - Strand 1
        const baseGeom = new THREE.SphereGeometry(0.5, 32, 32);
        const baseMat = new THREE.MeshPhongMaterial({
            color: baseColors[base] || 0x999999,
            shininess: 80
        });
        const baseMesh = new THREE.Mesh(baseGeom, baseMat);
        baseMesh.position.set(x1, y, z1);
        baseMesh.userData = { base, index: i, strand: "5' strand" };
        baseMesh.castShadow = true;
        dnaHelixGroup.add(baseMesh);

        // Add text label sprite - Strand 1
        const baseLabel = createTextSprite(base);
        baseLabel.position.set(x1, y, z1);
        dnaHelixGroup.add(baseLabel);

        // Nucleotide sphere - Strand 2
        const compGeom = new THREE.SphereGeometry(0.5, 32, 32);
        const compMat = new THREE.MeshPhongMaterial({
            color: baseColors[complement] || 0x999999,
            shininess: 80
        });
        const compMesh = new THREE.Mesh(compGeom, compMat);
        compMesh.position.set(x2, y, z2);
        compMesh.userData = { base: complement, index: i, strand: "3' strand" };
        compMesh.castShadow = true;
        dnaHelixGroup.add(compMesh);

        // Add text label sprite - Strand 2
        const compLabel = createTextSprite(complement);
        compLabel.position.set(x2, y, z2);
        dnaHelixGroup.add(compLabel);

        // Hydrogen bond (connecting line)
        const bondGeom = new THREE.CylinderGeometry(0.05, 0.05, helixRadius * 2, 8);
        const bondMat = new THREE.MeshPhongMaterial({ color: 0xcccccc, transparent: true, opacity: 0.5 });
        const bondMesh = new THREE.Mesh(bondGeom, bondMat);
        bondMesh.position.set(0, y, 0);
        bondMesh.rotation.z = Math.PI / 2;
        bondMesh.rotation.y = angle;
        dnaHelixGroup.add(bondMesh);
    }

    // Backbone ribbons
    const backboneCurve1 = [];
    const backboneCurve2 = [];
    for (let i = 0; i <= basePairs * 10; i++) {
        const t = i / 10;
        const angle = t * 0.6;
        const y = (t - basePairs / 2) * helixPitch;
        backboneCurve1.push(new THREE.Vector3(
            Math.cos(angle) * helixRadius,
            y,
            Math.sin(angle) * helixRadius
        ));
        backboneCurve2.push(new THREE.Vector3(
            Math.cos(angle + Math.PI) * helixRadius,
            y,
            Math.sin(angle + Math.PI) * helixRadius
        ));
    }

    const tubeGeom1 = new THREE.TubeGeometry(
        new THREE.CatmullRomCurve3(backboneCurve1), 200, 0.15, 8, false
    );
    const tubeMat1 = new THREE.MeshPhongMaterial({ color: 0x2563eb, shininess: 50 });
    dnaHelixGroup.add(new THREE.Mesh(tubeGeom1, tubeMat1));

    const tubeGeom2 = new THREE.TubeGeometry(
        new THREE.CatmullRomCurve3(backboneCurve2), 200, 0.15, 8, false
    );
    const tubeMat2 = new THREE.MeshPhongMaterial({ color: 0xdc2626, shininess: 50 });
    dnaHelixGroup.add(new THREE.Mesh(tubeGeom2, tubeMat2));

    dnaScene.add(dnaHelixGroup);

    // Click interaction - show tooltip (non-blocking)
    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();

    container.addEventListener('click', (event) => {
        const rect = container.getBoundingClientRect();
        mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
        mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

        raycaster.setFromCamera(mouse, dnaCamera);
        const intersects = raycaster.intersectObjects(dnaHelixGroup.children);

        if (intersects.length > 0) {
            const obj = intersects[0].object;
            if (obj.userData.base) {
                const baseNames = { A: 'Adenine', T: 'Thymine', G: 'Guanine', C: 'Cytosine' };
                dnaTooltip.innerHTML = `
                    <strong>${obj.userData.base}</strong> - ${baseNames[obj.userData.base]}<br>
                    Position: ${obj.userData.index + 1}<br>
                    Strand: ${obj.userData.strand}
                `;
                dnaTooltip.style.left = (event.clientX + 15) + 'px';
                dnaTooltip.style.top = (event.clientY - 10) + 'px';
                dnaTooltip.style.opacity = '1';
            }
        } else {
            // Clicked in container but not on nucleotide - hide tooltip
            dnaTooltip.style.opacity = '0';
        }
        event.stopPropagation();
    });

    // Hide tooltip immediately when clicking anywhere else
    document.addEventListener('click', () => {
        dnaTooltip.style.opacity = '0';
    });

    // Animation loop
    function animate() {
        dnaAnimationId = requestAnimationFrame(animate);
        dnaHelixGroup.rotation.y += helixRotationSpeed;
        dnaControls.update();
        dnaRenderer.render(dnaScene, dnaCamera);
    }
    animate();

    // Handle resize
    window.addEventListener('resize', () => {
        const newWidth = container.offsetWidth;
        const newHeight = container.offsetHeight || 400;
        dnaCamera.aspect = newWidth / newHeight;
        dnaCamera.updateProjectionMatrix();
        dnaRenderer.setSize(newWidth, newHeight);
    });
}

function updateHelixSpeed(val) {
    helixRotationSpeed = val / 1000;
}

function resetHelixView() {
    if (dnaCamera && dnaControls) {
        dnaCamera.position.set(0, 5, 25);
        dnaControls.reset();
    }
}

// ===== DESIGN GUIDES =====
async function designGuides() {
    const gene = document.getElementById('gene').value.trim();
    if (!gene) return;

    saveRecentSearch(gene);
    document.getElementById('recentDropdown').classList.remove('visible');

    document.getElementById('loading').classList.add('active');
    document.getElementById('designBtn').disabled = true;
    document.getElementById('summary').style.display = 'none';
    document.getElementById('visualizations').style.display = 'none';
    document.getElementById('sequence-viz').style.display = 'none';
    document.getElementById('guides-container').innerHTML = '';
    document.getElementById('status').style.display = 'none';

    // Show progress
    updateProgress(0, 'Finding candidates...');

    currentDesignData = { guides: [], allGuides: [], target: {} };
    currentPage = 1;
    selectedGuides.clear();

    try {
        const response = await fetch('/api/design', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                gene,
                modality: document.getElementById('modality').value,
                genome: document.getElementById('genome').value,
                cas_variant: document.getElementById('cas_variant').value,
                n_guides: 999999
            })
        });

        if (!response.ok) {
            const err = await response.json();
            // Handle Pydantic validation errors (array) vs regular errors (string)
            let errorMsg = 'Design failed';
            if (err.detail) {
                if (Array.isArray(err.detail)) {
                    // Pydantic validation error - extract messages
                    errorMsg = err.detail.map(e => e.msg || e.message || JSON.stringify(e)).join('; ');
                } else if (typeof err.detail === 'string') {
                    errorMsg = err.detail;
                } else {
                    errorMsg = JSON.stringify(err.detail);
                }
            }
            throw new Error(errorMsg);
        }

        const data = await response.json();

        currentDesignData.target = data.target || {};
        currentDesignData.guides = data.guides || [];
        currentDesignData.allGuides = data.guides || [];
        currentDesignData.n_guides_found = data.guides?.length || 0;

        updateProgress(100, `Found ${currentDesignData.guides.length} candidates`);
        setTimeout(() => finalizeResults(), 300);

    } catch (err) {
        showError(err.message);
    } finally {
        document.getElementById('loading').classList.remove('active');
        document.getElementById('designBtn').disabled = false;
    }
}

function updateProgress(percent, text) {
    const progressFill = document.getElementById('progressFill');
    const progressText = document.getElementById('progressText');
    if (progressFill) progressFill.style.width = percent + '%';
    if (progressText) progressText.textContent = text;
}

function finalizeResults() {
    if (!currentDesignData.guides || currentDesignData.guides.length === 0) {
        showError('No guides found for this target');
        return;
    }

    document.getElementById('summary').style.display = 'grid';
    document.getElementById('stat-gene').textContent = currentDesignData.target?.gene || '-';
    document.getElementById('stat-modality').textContent = currentDesignData.target?.modality || '-';

    const avgEff = currentDesignData.guides.reduce((s, g) => s + (g.efficiency_score || 0), 0) / currentDesignData.guides.length;
    document.getElementById('stat-efficiency').textContent = (avgEff * 100).toFixed(0) + '%';
    document.getElementById('stat-total').textContent = currentDesignData.guides.length.toLocaleString();

    const status = document.getElementById('status');
    status.innerHTML = `<span class="status-dot"></span>Found ${currentDesignData.guides.length.toLocaleString()} candidates`;
    status.style.display = 'inline-flex';

    initVisualizations();
    renderGuides();
}

function initVisualizations() {
    document.getElementById('visualizations').style.display = 'flex';
    document.getElementById('sequence-viz').style.display = 'block';

    if (!efficiencySketch) {
        efficiencySketch = new p5(efficiencySketchFn, 'efficiency-chart');
    }
    setTimeout(() => efficiencySketch?.updateData?.(currentDesignData.guides), 50);

    if (!genomeSketch) {
        genomeSketch = new p5(genomeSketchFn, 'genome-map');
    }
    setTimeout(() => genomeSketch?.updateData?.(currentDesignData.guides, currentDesignData.target), 50);

    if (!sequenceSketch) {
        sequenceSketch = new p5(sequenceSketchFn, 'sequence-canvas');
    }
    setTimeout(() => {
        if (currentDesignData.guides[0]) {
            const seq = currentDesignData.guides[0].spacer_sequence + currentDesignData.guides[0].pam_sequence;
            sequenceSketch?.updateSequence?.(seq);
        }
    }, 50);
}

// ===== RENDER GUIDES =====
function renderGuides() {
    let filtered = [...currentDesignData.guides];

    // Apply filters
    if (searchQuery) {
        filtered = filtered.filter(g => g.spacer_sequence.includes(searchQuery.toUpperCase()));
    }
    if (strandFilter !== 'all') {
        filtered = filtered.filter(g => g.strand === strandFilter);
    }
    if (efficiencyMin > 0 || efficiencyMax < 100) {
        filtered = filtered.filter(g => {
            const eff = g.efficiency_score * 100;
            return eff >= efficiencyMin && eff <= efficiencyMax;
        });
    }

    // Paginate
    const start = (currentPage - 1) * perPage;
    const pageGuides = filtered.slice(start, start + perPage);
    const totalPages = Math.ceil(filtered.length / perPage) || 1;

    const container = document.getElementById('guides-container');
    container.innerHTML = pageGuides.map(g => renderGuideCard(g)).join('');

    // Update pagination
    document.getElementById('pageInfo').textContent = `Page ${currentPage} of ${totalPages}`;
    document.getElementById('prevBtn').disabled = currentPage <= 1;
    document.getElementById('nextBtn').disabled = currentPage >= totalPages;

    updateCompareBtn();
}

function renderGuideCard(g) {
    const checked = selectedGuides.has(g.rank) ? 'checked' : '';
    const gcContent = calculateGC(g.spacer_sequence);

    return `
        <div class="guide-card" data-rank="${g.rank}" title="GC: ${gcContent}%">
            <div class="guide-header">
                <div class="guide-left">
                    <input type="checkbox" class="guide-checkbox" ${checked} onchange="toggleGuide(${g.rank})">
                    <div class="guide-rank">${g.rank}</div>
                    <div class="guide-sequence">${colorize(g.spacer_sequence)}</div>
                </div>
                <div class="guide-actions">
                    <button class="copy-btn" onclick="copySeq('${g.spacer_sequence}')">Copy</button>
                    <span class="score-badge ${scoreClass(g.efficiency_score)}">${(g.efficiency_score * 100).toFixed(0)}%</span>
                </div>
            </div>
            <div class="guide-details">
                <div class="detail-item"><div class="detail-label">PAM</div><div class="detail-value">${g.pam_sequence}</div></div>
                <div class="detail-item"><div class="detail-label">Strand</div><div class="detail-value">${g.strand}</div></div>
                <div class="detail-item"><div class="detail-label">Cut Site</div><div class="detail-value">${g.cut_site.toLocaleString()}</div></div>
                <div class="detail-item"><div class="detail-label">Position</div><div class="detail-value">${g.genomic_start.toLocaleString()}</div></div>
            </div>
        </div>
    `;
}

function colorize(seq) {
    return seq.split('').map(c => `<span class="nucleotide-${c}">${c}</span>`).join('');
}

function scoreClass(s) {
    if (s >= 0.8) return 'score-excellent';
    if (s >= 0.6) return 'score-good';
    return 'score-moderate';
}

function calculateGC(seq) {
    const gc = (seq.match(/[GC]/gi) || []).length;
    return ((gc / seq.length) * 100).toFixed(0);
}

// ===== GUIDE ACTIONS =====
function toggleGuide(rank) {
    if (selectedGuides.has(rank)) selectedGuides.delete(rank);
    else selectedGuides.add(rank);
    updateCompareBtn();
}

function updateCompareBtn() {
    const btn = document.getElementById('compareBtn');
    btn.textContent = `Compare (${selectedGuides.size})`;
    btn.disabled = selectedGuides.size < 2;

    const bulkBtn = document.getElementById('bulkCopyBtn');
    if (bulkBtn) {
        bulkBtn.style.display = selectedGuides.size > 0 ? 'inline-block' : 'none';
    }
}

function copySeq(seq) {
    navigator.clipboard.writeText(seq);
}

function bulkCopySelected() {
    const seqs = currentDesignData.guides
        .filter(g => selectedGuides.has(g.rank))
        .map(g => g.spacer_sequence);
    navigator.clipboard.writeText(seqs.join('\n'));
}

// ===== PAGINATION =====
function prevPage() {
    if (currentPage > 1) {
        currentPage--;
        renderGuides();
        // Don't auto-scroll - stay where user is
    }
}

function nextPage() {
    currentPage++;
    renderGuides();
    // Don't auto-scroll - stay where user is
}

function setPerPage(val) {
    perPage = parseInt(val);
    currentPage = 1;
    renderGuides();
}

function scrollToResults() {
    document.getElementById('guides-container').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// ===== INFINITE SCROLL =====
function setupInfiniteScroll() {
    // Listen on window scroll instead of main container
    window.addEventListener('scroll', () => {
        if (!infiniteScrollEnabled || isLoadingMore) return;

        const scrollY = window.scrollY;
        const windowHeight = window.innerHeight;
        const docHeight = document.documentElement.scrollHeight;

        if (scrollY + windowHeight >= docHeight - 300) {
            loadMoreGuides();
        }
    });
}

function toggleInfiniteScroll() {
    infiniteScrollEnabled = document.getElementById('infiniteScrollToggle').checked;
    document.querySelector('.pagination-controls').style.opacity = infiniteScrollEnabled ? '0.5' : '1';
}

function loadMoreGuides() {
    if (isLoadingMore) return;

    const filtered = getFilteredGuides();
    const currentlyShown = currentPage * perPage;

    if (currentlyShown >= filtered.length) return;

    isLoadingMore = true;

    // Load next batch
    const nextBatch = filtered.slice(currentlyShown, currentlyShown + INFINITE_SCROLL_BATCH);
    const container = document.getElementById('guides-container');

    nextBatch.forEach(g => {
        container.insertAdjacentHTML('beforeend', renderGuideCard(g));
    });

    currentPage = Math.ceil((currentlyShown + nextBatch.length) / perPage);
    isLoadingMore = false;
}

function getFilteredGuides() {
    let filtered = [...currentDesignData.guides];
    if (searchQuery) filtered = filtered.filter(g => g.spacer_sequence.includes(searchQuery.toUpperCase()));
    if (strandFilter !== 'all') filtered = filtered.filter(g => g.strand === strandFilter);
    if (efficiencyMin > 0 || efficiencyMax < 100) {
        filtered = filtered.filter(g => {
            const eff = g.efficiency_score * 100;
            return eff >= efficiencyMin && eff <= efficiencyMax;
        });
    }
    return filtered;
}

// ===== FILTERS =====
function applySearch() {
    searchQuery = document.getElementById('searchInput').value;
    currentPage = 1;
    renderGuides();
}

function applyStrandFilter(strand) {
    strandFilter = strand;
    currentPage = 1;
    renderGuides();
}

function applyEfficiencyFilter() {
    efficiencyMin = parseInt(document.getElementById('effMinSlider').value) || 0;
    efficiencyMax = parseInt(document.getElementById('effMaxSlider').value) || 100;
    document.getElementById('effMinVal').textContent = efficiencyMin + '%';
    document.getElementById('effMaxVal').textContent = efficiencyMax + '%';
    currentPage = 1;
    renderGuides();
}

// ===== EXPAND FEATURES =====
function expandEfficiency() {
    document.getElementById('efficiencyModal').classList.add('active');
    // Initialize fullscreen chart
    setTimeout(() => {
        const container = document.getElementById('efficiency-chart-fullscreen');
        if (container && efficiencySketch) {
            // Create new sketch for fullscreen
            new p5((p) => {
                let bins = [];
                p.setup = () => {
                    p.createCanvas(container.offsetWidth - 20, 450);
                    bins = new Array(10).fill(0);
                    currentDesignData.guides.forEach(g => {
                        const idx = Math.min(9, Math.floor((g.efficiency_score || 0) * 10));
                        bins[idx]++;
                    });
                };
                p.draw = () => {
                    p.background(255);
                    const maxCount = Math.max(...bins, 1);
                    const barW = (p.width - 80) / bins.length;
                    const maxH = p.height - 80;

                    p.stroke(226, 232, 240);
                    p.line(50, 20, 50, p.height - 40);
                    p.line(50, p.height - 40, p.width - 20, p.height - 40);

                    p.noStroke();
                    for (let i = 0; i < bins.length; i++) {
                        const h = (bins[i] / maxCount) * maxH;
                        const x = 60 + i * barW;
                        const pct = (i + 0.5) * 10;
                        if (pct >= 80) p.fill(5, 150, 105);
                        else if (pct >= 60) p.fill(8, 145, 178);
                        else p.fill(217, 119, 6);
                        p.rect(x, p.height - 40 - h, barW - 4, h, 3);

                        // Count label on bar
                        p.fill(71, 85, 105);
                        p.textAlign(p.CENTER);
                        p.textSize(11);
                        if (bins[i] > 0) p.text(bins[i], x + barW / 2 - 2, p.height - 45 - h);
                    }

                    p.fill(100, 116, 139);
                    p.textSize(12);
                    for (let i = 0; i <= 10; i += 2) {
                        p.text(i * 10 + '%', 60 + i * barW, p.height - 20);
                    }
                };
            }, 'efficiency-chart-fullscreen');
        }
    }, 100);
}

function closeEfficiencyModal() {
    document.getElementById('efficiencyModal').classList.remove('active');
}

function expandGenomeMap() {
    document.getElementById('genomeModal').classList.add('active');
    setTimeout(() => {
        const container = document.getElementById('genome-map-fullscreen');
        if (container) {
            container.innerHTML = '';
            new p5((p) => {
                let guides = currentDesignData.guides.slice(0, 50);
                let range = { start: currentDesignData.target?.start || 0, end: currentDesignData.target?.end || 1000000 };

                p.setup = () => {
                    p.createCanvas(container.offsetWidth - 20, 420);
                };

                p.draw = () => {
                    p.background(255);
                    const yMid = 210;
                    const barH = 28;

                    p.fill(241, 245, 249);
                    p.stroke(200, 210, 220);
                    p.rect(50, yMid - barH / 2, p.width - 100, barH, barH / 2);

                    p.fill(100, 116, 139);
                    p.noStroke();
                    p.textSize(12);
                    p.textAlign(p.LEFT);
                    p.text(range.start.toLocaleString(), 50, yMid + barH / 2 + 30);
                    p.textAlign(p.RIGHT);
                    p.text(range.end.toLocaleString(), p.width - 50, yMid + barH / 2 + 30);

                    const rng = range.end - range.start;
                    const positions = guides.map(g => ({
                        x: 50 + ((g.cut_site - range.start) / rng) * (p.width - 100),
                        g
                    })).sort((a, b) => a.x - b.x);

                    let lastX = -50;
                    let alternateUp = false;

                    positions.forEach(pos => {
                        let { x, g } = pos;
                        let dotY = yMid - barH - 15;
                        if (x - lastX < 15) {
                            alternateUp = !alternateUp;
                            dotY = alternateUp ? yMid - barH - 30 : yMid - barH - 8;
                        } else {
                            alternateUp = false;
                        }
                        lastX = x;

                        p.stroke(49, 130, 206);
                        p.strokeWeight(2);
                        p.line(x, dotY + 5, x, yMid + barH / 2 - 2);

                        p.noStroke();
                        p.fill(g.strand === '+' ? [34, 197, 94] : [239, 68, 68]);
                        p.ellipse(x, dotY, 12, 12);

                        p.fill(26, 54, 93);
                        p.textAlign(p.CENTER);
                        p.textSize(10);
                        p.text(g.rank, x, dotY - 15);
                    });

                    p.noLoop();
                };
            }, 'genome-map-fullscreen');
        }
    }, 100);
}

function closeGenomeModal() {
    document.getElementById('genomeModal').classList.remove('active');
}

function expandSequence() {
    document.getElementById('sequenceModal').classList.add('active');
    setTimeout(() => {
        const container = document.getElementById('sequence-canvas-fullscreen');
        if (container && currentDesignData.guides[0]) {
            container.innerHTML = '';
            const seq = currentDesignData.guides[0].spacer_sequence + currentDesignData.guides[0].pam_sequence;
            new p5((p) => {
                p.setup = () => {
                    p.createCanvas(container.offsetWidth - 20, 120);
                    p.textFont('monospace');
                };
                p.draw = () => {
                    p.background(250, 250, 248); // Warm neutral
                    const charW = 22;
                    const y = 60;
                    const pamStart = seq.length - 3;

                    for (let i = 0; i < seq.length; i++) {
                        const ch = seq[i];
                        const x = 20 + i * charW;

                        if (i >= pamStart) {
                            p.fill(168, 85, 247, 40);
                            p.noStroke();
                            p.rect(x - 2, y - 20, charW, 32, 3);
                        }

                        if (ch === 'A') p.fill('#22c55e');
                        else if (ch === 'T') p.fill('#ef4444');
                        else if (ch === 'G') p.fill('#eab308');
                        else if (ch === 'C') p.fill('#3b82f6');

                        p.textSize(24);
                        p.textAlign(p.CENTER);
                        p.noStroke();
                        p.text(ch, x + charW / 2, y);

                        p.fill(100);
                        p.textSize(10);
                        p.text(i + 1, x + charW / 2, y + 25);
                    }

                    p.textSize(14);
                    p.fill(168, 85, 247);
                    p.textAlign(p.RIGHT);
                    p.text('PAM', p.width - 20, 25);
                    p.noLoop();
                };
            }, 'sequence-canvas-fullscreen');
        }
    }, 100);
}

function closeSequenceModal() {
    document.getElementById('sequenceModal').classList.remove('active');
}

// ===== COMPARE MODE =====
function openCompare() {
    if (selectedGuides.size < 2) return;

    const selected = currentDesignData.guides.filter(g => selectedGuides.has(g.rank));

    // Build alignment
    let alignmentHtml = selected.map(g => `
        <div class="compare-alignment-row">
            <span class="rank">#${g.rank}</span>
            <span>${colorize(g.spacer_sequence)}</span>
        </div>
    `).join('');

    // Build table
    let tableHtml = `
        <table class="compare-table">
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Efficiency</th>
                    <th>GC%</th>
                    <th>Strand</th>
                    <th>Position</th>
                </tr>
            </thead>
            <tbody>
                ${selected.map(g => `
                    <tr>
                        <td>#${g.rank}</td>
                        <td>${(g.efficiency_score * 100).toFixed(0)}%</td>
                        <td>${calculateGC(g.spacer_sequence)}%</td>
                        <td>${g.strand}</td>
                        <td>${g.genomic_start.toLocaleString()}</td>
                    </tr>
                `).join('')}
            </tbody>
        </table>
    `;

    document.getElementById('compareAlignment').innerHTML = alignmentHtml;
    document.getElementById('compareTable').innerHTML = tableHtml;
    document.getElementById('compareModal').classList.add('active');
}

function closeCompareModal() {
    document.getElementById('compareModal').classList.remove('active');
}

// ===== EXPORT =====
function openExportDialog() {
    document.getElementById('exportModal').classList.add('active');
    document.getElementById('exportCount').max = currentDesignData.guides.length;
    document.getElementById('exportCount').value = Math.min(100, currentDesignData.guides.length);
}

function closeExportModal() {
    document.getElementById('exportModal').classList.remove('active');
}

function setExportCount(n) {
    document.getElementById('exportCount').value = n === 'all' ? currentDesignData.guides.length : n;
    document.querySelectorAll('.preset-btn').forEach(b => b.classList.remove('active'));
    event.target.classList.add('active');
}

async function exportPDF() {
    const count = parseInt(document.getElementById('exportCount').value) || 10;
    const minEff = parseInt(document.getElementById('exportMinEff').value) || 0;
    // Use allGuides to ensure we have the complete unfiltered list
    let guides = [...(currentDesignData.allGuides || currentDesignData.guides)];
    console.log(`Exporting PDF: Requested ${count}, Available ${guides.length}`);
    if (minEff > 0) guides = guides.filter(g => g.efficiency_score >= minEff / 100);
    guides = guides.slice(0, count);

    const content = generatePDFContent(guides);
    const blob = new Blob([content], { type: 'text/html' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${currentDesignData.target.gene || 'guides'}_report.html`;
    a.click();
    closeExportModal();
}

function exportCSV() {
    const count = parseInt(document.getElementById('exportCount').value) || 10;
    const minEff = parseInt(document.getElementById('exportMinEff').value) || 0;
    // Use allGuides to ensure we have the complete unfiltered list
    let guides = [...(currentDesignData.allGuides || currentDesignData.guides)];
    console.log(`Exporting CSV: Requested ${count}, Available ${guides.length}`);
    if (minEff > 0) guides = guides.filter(g => g.efficiency_score >= minEff / 100);
    guides = guides.slice(0, count);

    const headers = ['Rank', 'Sequence', 'PAM', 'Strand', 'Cut Site', 'Position', 'Efficiency', 'GC%'];
    const rows = guides.map(g =>
        [g.rank, g.spacer_sequence, g.pam_sequence, g.strand, g.cut_site, g.genomic_start,
        (g.efficiency_score * 100).toFixed(1) + '%', calculateGC(g.spacer_sequence) + '%'].join(',')
    );
    const csv = [headers.join(','), ...rows].join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${currentDesignData.target.gene || 'guides'}.csv`;
    a.click();
    closeExportModal();
}

function generatePDFContent(guides) {
    const t = currentDesignData.target;
    return `<!DOCTYPE html><html><head><title>SCALPEL Report - ${t.gene}</title>
<style>body{font-family:Georgia,serif;max-width:900px;margin:auto;padding:2rem}
h1{color:#1a365d}table{width:100%;border-collapse:collapse;margin:1rem 0}
th,td{border:1px solid #ccc;padding:10px;text-align:left}th{background:#f1f5f9}
.seq{font-family:monospace;font-size:14px}</style></head><body>
<h1>SCALPEL gRNA Report</h1>
<p><strong>Gene:</strong> ${t.gene} | <strong>Modality:</strong> ${t.modality} | <strong>Genome:</strong> ${t.genome}</p>
<p><strong>Coordinates:</strong> ${t.chromosome}:${t.start?.toLocaleString()}-${t.end?.toLocaleString()}</p>
<p><strong>Total Candidates:</strong> ${(currentDesignData.allGuides || currentDesignData.guides).length.toLocaleString()} | <strong>Exported:</strong> ${guides.length}</p>
<h2>Guides</h2>
<table><thead><tr><th>#</th><th>Sequence</th><th>PAM</th><th>Strand</th><th>Cut Site</th><th>Efficiency</th><th>GC%</th></tr></thead>
<tbody>${guides.map(g => `<tr><td>${g.rank}</td><td class="seq">${g.spacer_sequence}</td><td>${g.pam_sequence}</td><td>${g.strand}</td><td>${g.cut_site.toLocaleString()}</td><td>${(g.efficiency_score * 100).toFixed(0)}%</td><td>${calculateGC(g.spacer_sequence)}%</td></tr>`).join('')}</tbody></table>
<p style="color:#666;font-size:12px;margin-top:2rem">Generated by SCALPEL - ${new Date().toISOString()}</p></body></html>`;
}

// ===== OFF-TARGET & PLAN =====
async function analyzeOffTarget(spacer) {
    const c = document.getElementById('results-offtarget');
    c.innerHTML = '<div class="loading active"><div class="spinner"></div><p>Analyzing off-targets...</p></div>';
    try {
        const r = await fetch('/api/offtarget', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ spacer, genome: document.getElementById('genome').value }) });
        const d = await r.json();
        c.innerHTML = `<div style="padding:1rem"><h3 style="margin-bottom:1rem">Off-Target Analysis: ${spacer}</h3><p>Total Sites: ${d.total_sites} | Specificity: ${d.specificity_score}%</p></div>`;
    } catch (e) { c.innerHTML = `<p style="padding:1rem;color:#dc2626">Error: ${e.message}</p>`; }
}

async function generatePlan() {
    const c = document.getElementById('results-plan');
    const t = currentDesignData.target;
    const guides = currentDesignData.guides;
    const topGuides = guides.slice(0, 5);
    const avgEff = guides.length > 0 ? (guides.reduce((a, g) => a + g.efficiency_score, 0) / guides.length * 100).toFixed(1) : 0;

    // Generate comprehensive experiment plan
    c.innerHTML = `
        <div style="padding:1.5rem;max-width:900px;">
            <h2 style="color:#1a365d;font-family:Georgia,serif;margin-bottom:0.5rem;">Experiment Plan</h2>
            <p style="color:#64748b;margin-bottom:1.5rem;">CRISPR ${t.modality || 'knockout'} experiment for <strong>${t.gene || 'target gene'}</strong></p>
            
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;margin-bottom:1.5rem;">
                <div style="background:#f8fafc;padding:1rem;border-radius:8px;border:1px solid #e2e8f0;">
                    <h4 style="color:#1a365d;margin-bottom:0.75rem;font-size:0.9rem;text-transform:uppercase;">Target Information</h4>
                    <table style="font-size:0.85rem;color:#475569;">
                        <tr><td style="padding:0.3rem 0;">Gene:</td><td style="font-weight:600;">${t.gene || 'N/A'}</td></tr>
                        <tr><td style="padding:0.3rem 0;">Genome:</td><td>${t.genome || 'GRCh38'}</td></tr>
                        <tr><td style="padding:0.3rem 0;">Chromosome:</td><td>${t.chromosome || 'N/A'}</td></tr>
                        <tr><td style="padding:0.3rem 0;">Position:</td><td>${t.start?.toLocaleString() || 'N/A'} - ${t.end?.toLocaleString() || 'N/A'}</td></tr>
                    </table>
                </div>
                <div style="background:#f8fafc;padding:1rem;border-radius:8px;border:1px solid #e2e8f0;">
                    <h4 style="color:#1a365d;margin-bottom:0.75rem;font-size:0.9rem;text-transform:uppercase;">Design Summary</h4>
                    <table style="font-size:0.85rem;color:#475569;">
                        <tr><td style="padding:0.3rem 0;">Modality:</td><td style="font-weight:600;">${t.modality || 'knockout'}</td></tr>
                        <tr><td style="padding:0.3rem 0;">Total Candidates:</td><td>${guides.length.toLocaleString()}</td></tr>
                        <tr><td style="padding:0.3rem 0;">Avg Efficiency:</td><td>${avgEff}%</td></tr>
                        <tr><td style="padding:0.3rem 0;">Cas Variant:</td><td>SpCas9 (NGG PAM)</td></tr>
                    </table>
                </div>
            </div>
            
            <h3 style="color:#1a365d;border-bottom:2px solid #e2e8f0;padding-bottom:0.5rem;margin:1.5rem 0 1rem;">Recommended Guides</h3>
            <table style="width:100%;border-collapse:collapse;font-size:0.85rem;">
                <thead>
                    <tr style="background:#f1f5f9;">
                        <th style="padding:0.6rem;text-align:left;border:1px solid #e2e8f0;">Rank</th>
                        <th style="padding:0.6rem;text-align:left;border:1px solid #e2e8f0;">Sequence (5' → 3')</th>
                        <th style="padding:0.6rem;text-align:left;border:1px solid #e2e8f0;">PAM</th>
                        <th style="padding:0.6rem;text-align:center;border:1px solid #e2e8f0;">Strand</th>
                        <th style="padding:0.6rem;text-align:center;border:1px solid #e2e8f0;">Efficiency</th>
                        <th style="padding:0.6rem;text-align:center;border:1px solid #e2e8f0;">GC%</th>
                    </tr>
                </thead>
                <tbody>
                    ${topGuides.map(g => `
                        <tr>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;font-weight:600;">#${g.rank}</td>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;font-family:monospace;font-size:0.9rem;">${g.spacer_sequence}</td>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;font-family:monospace;">${g.pam_sequence || 'NGG'}</td>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;text-align:center;">${g.strand}</td>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;text-align:center;color:${g.efficiency_score > 0.7 ? '#059669' : '#d97706'};">${(g.efficiency_score * 100).toFixed(0)}%</td>
                            <td style="padding:0.6rem;border:1px solid #e2e8f0;text-align:center;">${calculateGC(g.spacer_sequence)}%</td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
            
            <h3 style="color:#1a365d;border-bottom:2px solid #e2e8f0;padding-bottom:0.5rem;margin:1.5rem 0 1rem;">Materials Required</h3>
            <ul style="color:#475569;font-size:0.85rem;line-height:1.8;padding-left:1.5rem;">
                <li>SpCas9 protein or expression plasmid</li>
                <li>sgRNA (synthesized or expressed from plasmid) - see sequences above</li>
                <li>Cell line or organism for ${t.gene || 'target gene'}</li>
                <li>Transfection reagent (Lipofectamine, RNP delivery, or electroporation)</li>
                <li>Selection/screening reagents (if applicable)</li>
                <li>PCR primers for genotyping the target locus</li>
            </ul>
            
            <h3 style="color:#1a365d;border-bottom:2px solid #e2e8f0;padding-bottom:0.5rem;margin:1.5rem 0 1rem;">Protocol Overview</h3>
            <ol style="color:#475569;font-size:0.85rem;line-height:1.8;padding-left:1.5rem;">
                <li><strong>Day 0:</strong> Prepare cells at 70-80% confluence</li>
                <li><strong>Day 1:</strong> Transfect with Cas9/sgRNA (RNP or plasmid)</li>
                <li><strong>Day 2-3:</strong> Change media, monitor cell health</li>
                <li><strong>Day 4-7:</strong> Harvest cells for genomic DNA extraction</li>
                <li><strong>Day 7-10:</strong> PCR and sequencing to confirm editing</li>
            </ol>
        </div>
    `;
}

// ===== UTILS =====
function showError(msg) {
    const s = document.getElementById('status');
    s.textContent = msg;
    s.className = 'status-indicator';
    s.style.cssText = 'display:inline-flex;background:rgba(220,38,38,0.1);border:1px solid rgba(220,38,38,0.2);color:#dc2626';
}
