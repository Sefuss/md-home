# MD Ensemble Approach - Setup Plan
**Date:** 2025-12-16
**Goal:** Generate conformational ensemble of AviTag peptide for binder design

---

## What We're Doing

**The Strategy:**
1. Run MD simulation of free AviTag peptide (100-200 ns)
2. Sample all conformations it naturally adopts
3. Cluster to find dominant states
4. Design binders for those real conformations using RFD3

**Why:** Avoids designing for artificial structures

---

## Software Stack Needed

### 1. GROMACS (Molecular Dynamics)
**Purpose:** Run MD simulation
**Installation:** Docker (easiest for Windows)
**Image:** `gromacs/gromacs:latest` (official)

### 2. MDAnalysis (Trajectory Analysis)
**Purpose:** Analyze MD trajectory, cluster conformations
**Installation:** pip via conda Python
**Command:** `pip install MDAnalysis`

### 3. BioPython (Structure Building)
**Purpose:** Build initial peptide structure
**Status:** ✅ Already installed

### 4. Sklearn (Clustering)
**Purpose:** RMSD-based conformational clustering
**Status:** Should be in anaconda (check)

---

## Workflow Overview

```
Step 1: Build Peptide Structure (BioPython)
   └─> avitag_initial.pdb (15 aa extended)

Step 2: GROMACS Preparation
   ├─> Add hydrogens (pdb2gmx)
   ├─> Solvate in water box
   ├─> Add ions (neutralize)
   └─> Energy minimize

Step 3: Equilibration
   ├─> NVT (100 ps) - temperature
   └─> NPT (100 ps) - pressure

Step 4: Production MD (100-200 ns)
   └─> Save frames every 10 ps
   └─> Output: trajectory.xtc (~10,000 frames)

Step 5: Analysis (MDAnalysis)
   ├─> Calculate RMSD
   ├─> Cluster by conformation (K-means or DBSCAN)
   └─> Extract cluster representatives

Step 6: Select for RFD3
   ├─> Top 1-3 dominant conformations
   └─> Create PDB files for each

Step 7: RFD3 Design
   ├─> Design binder for each conformation
   └─> Or ensemble-based design
```

---

## Timeline Estimate

| Step | Time (CPU) | Time (GPU) |
|------|------------|------------|
| Setup & building | 1 hour | 1 hour |
| Energy minimization | 10 min | 5 min |
| Equilibration | 30 min | 10 min |
| Production MD (100 ns) | 12-24 hours | 2-4 hours |
| Analysis | 30 min | 30 min |
| **Total** | **~1-2 days** | **~4-6 hours** |

**Note:** For peptide (15 aa), even CPU should be reasonable

---

## Hardware Requirements

**For 15 aa peptide in water:**
- System size: ~5,000-10,000 atoms (peptide + water)
- RAM: 4-8 GB
- Storage: 1-2 GB for trajectory
- CPU: Works (slow but doable)
- GPU: Much faster if available

**Your system (i3-10100):**
- Should handle this fine
- Expect ~24 hours for 100 ns on CPU
- Consider shorter run (50 ns) for initial test

---

## Parameters We'll Use

**Force Field:** AMBER99SB-ILDN or CHARMM36m
- Good for peptides/proteins
- Well-validated

**Water Model:** TIP3P
- Standard, fast

**Box Size:** 1.0 nm padding
- Small system, fast simulation

**Temperature:** 300 K (physiological)

**Pressure:** 1 bar

**Time Step:** 2 fs (with constraints)

**Production Length:** 100-200 ns
- 100 ns for initial test
- 200 ns if need more sampling

---

## Expected Outcomes

**What MD will tell us:**

1. **Conformational States**
   - How many distinct shapes does AviTag adopt?
   - Are they well-defined or continuous?

2. **Dominant Conformation**
   - Is there a single preferred shape? (best case)
   - Or equally distributed? (harder case)

3. **Structural Features**
   - Any transient secondary structure?
   - Turns, bends, or extended regions?
   - Hydrophobic collapse?

**Possible results:**

**Scenario A: Single dominant state (40%+ population)**
```
Perfect! Design binder for that state
Confidence: High
```

**Scenario B: 2-3 major states (20-30% each)**
```
Good. Design for top 2-3
Or try multi-specific binder
Confidence: Medium
```

**Scenario C: Highly dynamic (no dominant state)**
```
Challenging. Options:
- Design for transition states
- Use ensemble-average
- Consider experimental route
Confidence: Low
```

---

## Files We'll Generate

**Structure files:**
```
avitag_initial.pdb          # Starting structure
avitag_solvated.gro         # After solvation
avitag_minimized.gro        # After energy min
avitag_equilibrated.gro     # After equilibration
```

**Trajectory files:**
```
md_trajectory.xtc           # Full trajectory (compressed)
md_trajectory.trr           # Full precision (optional)
md.log                      # Energy, temp, pressure log
```

**Analysis outputs:**
```
rmsd_plot.png               # RMSD over time
cluster_populations.json    # Cluster statistics
cluster_0.pdb               # Dominant conformation
cluster_1.pdb               # 2nd most common
cluster_2.pdb               # 3rd most common
ensemble_summary.txt        # Analysis report
```

**For RFD3:**
```
02_scaffolds/approach_C_peptide_ensemble/
├── cluster_representatives/
│   ├── cluster_0_dominant.pdb
│   ├── cluster_1.pdb
│   └── cluster_2.pdb
├── rfd3_configs/
│   ├── design_cluster0.json
│   ├── design_cluster1.json
│   └── design_ensemble.json
└── analysis/
    └── conformational_analysis.md
```

---

## Next Steps

1. Pull GROMACS Docker image
2. Install MDAnalysis
3. Build initial AviTag structure
4. Set up MD simulation files
5. Run test (short 10 ns) to validate setup
6. If good → run full 100 ns
7. Analyze and extract conformations
8. Feed to RFD3

---

**Created:** 2025-12-16
**Status:** Ready to begin setup
