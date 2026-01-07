# Next Steps: From MD Analysis to Binder Design

**Status:** MD complete (22.7 ns), ready for clustering ✅
**Date:** January 2, 2026

---

## Current Status

✅ **Completed:**
- MD simulation (22.7 ns, high-quality data)
- Quality assessment (all metrics excellent)
- Energy/temperature analysis

📋 **Ready to Start:**
- Trajectory clustering
- Conformation extraction
- RFDiffusion3 binder design

---

## Workflow Overview

```
MD Trajectory (22.7 ns, 2,270 frames)
          ↓
    Clustering (Group similar conformations)
          ↓
    Extract Representatives (5-10 structures)
          ↓
    Visual Quality Check (PyMOL)
          ↓
    Prepare RFD3 Inputs (JSON configs)
          ↓
    Generate Binders (RFDiffusion3)
          ↓
    Validate Designs (ColabFold)
```

---

## Step-by-Step Guide

### Step 1: Stop the MD Simulation

**Check if still running:**
```bash
docker ps
```

**If running, stop it:**
```bash
docker stop <container_id>
```

**Verify final trajectory:**
```bash
cd C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation

ls -lh md.xtc
# Should be ~85 MB
```

---

### Step 2: Cluster the Trajectory

**Goal:** Group similar conformations and extract representatives

**Option A: GROMACS Clustering (Recommended)**

```bash
cd 01_md_preparation

# Method 1: GROMOS algorithm (distance-based)
export MSYS_NO_PATHCONV=1
echo -e "Backbone\nBackbone" | docker run --rm -i \
  -v "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation:/work" \
  -w /work \
  gromacs/gromacs:latest \
  gmx cluster -f md.xtc -s npt.gro \
  -method gromos -cutoff 0.25 \
  -cl clusters.pdb -o cluster_rmsd.xpm \
  -g cluster.log -clid cluster_over_time.xvg

# This will create:
# - clusters.pdb: PDB with all cluster centers
# - cluster.log: Clustering statistics
# - cluster_over_time.xvg: Which cluster each frame belongs to
```

**Parameters:**
- `-cutoff 0.25`: RMSD cutoff (Ångstroms) for grouping
  - Lower = more clusters (more diversity)
  - Higher = fewer clusters (more similar)
  - 0.2-0.3 Å is typical for peptides

**Option B: Python Clustering (More control)**

```python
# Use MDAnalysis + scikit-learn
import MDAnalysis as mda
from sklearn.cluster import KMeans
import numpy as np

# Load trajectory
u = mda.Universe('npt.gro', 'md.xtc')
backbone = u.select_atoms('backbone')

# Align all frames
from MDAnalysis.analysis import align
align.AlignTraj(u, u, select='backbone', in_memory=True).run()

# Extract coordinates
coords = []
for ts in u.trajectory:
    coords.append(backbone.positions.flatten())
coords = np.array(coords)

# K-means clustering (5-10 clusters)
kmeans = KMeans(n_clusters=8, random_state=42)
clusters = kmeans.fit_predict(coords)

# Find frame closest to each cluster center
representatives = []
for i in range(8):
    cluster_frames = np.where(clusters == i)[0]
    cluster_coords = coords[cluster_frames]
    center = kmeans.cluster_centers_[i]

    # Find closest frame to center
    distances = np.linalg.norm(cluster_coords - center, axis=1)
    rep_idx = cluster_frames[np.argmin(distances)]
    representatives.append(rep_idx)

    print(f"Cluster {i}: {len(cluster_frames)} frames, representative: frame {rep_idx}")

# Save representatives
for i, frame_idx in enumerate(representatives):
    u.trajectory[frame_idx]
    with mda.Writer(f'cluster_{i}_rep.pdb') as W:
        W.write(backbone.residues.atoms)
```

---

### Step 3: Extract Representative Structures

**From GROMACS clustering:**

```bash
# clusters.pdb contains all cluster centers
# Need to split into individual PDBs

# Method 1: Manual splitting (if clusters.pdb has multiple MODELs)
# Use PyMOL or text editor to separate

# Method 2: Use gmx trjconv to extract specific frames
# First, find which frames are representatives from cluster.log

# Example: Extract frames 100, 500, 1200, etc.
export MSYS_NO_PATHCONV=1
echo "Backbone" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx trjconv -f md.xtc -s npt.gro -o cluster_1.pdb -dump 100
```

**Expected Output:**
- 5-10 PDB files, each representing a distinct conformation
- Files named: `cluster_0.pdb`, `cluster_1.pdb`, etc.

---

### Step 4: Visual Quality Check

**Load in PyMOL:**

```bash
pymol cluster_0.pdb cluster_1.pdb cluster_2.pdb cluster_3.pdb cluster_4.pdb
```

**What to check:**
- ✅ Peptide looks reasonable (not broken)
- ✅ No weird atom clashes
- ✅ Different conformations look distinct
- ✅ Structures are compact (not exploded)

**PyMOL commands for comparison:**
```python
# In PyMOL
align cluster_1, cluster_0
align cluster_2, cluster_0
# etc.

# Show cartoon representation
hide everything
show cartoon
```

**Document:**
- Which conformations look most different?
- Are there extended vs compact states?
- Any interesting secondary structure?

---

### Step 5: Analyze Cluster Statistics

**Questions to answer:**
1. How many unique conformations did we find? (5-10 expected)
2. Are they evenly distributed or are some major/minor?
3. What's the RMSD between different clusters?

**Check cluster populations:**
```bash
# From cluster.log or cluster_over_time.xvg
grep "Cluster" cluster.log
```

**Ideal scenario:**
- 5-10 distinct clusters
- Reasonably balanced populations (no single cluster dominates)
- RMSD between clusters > 2-3 Å (distinct enough)

---

### Step 6: Prepare for RFDiffusion3

**Create JSON configs for ensemble design:**

For each cluster representative, create a JSON file:

**Example: `rfd3_cluster_0.json`**
```json
{
  "inference": {
    "input_pdb": "../../../02_scaffolds/approach_C_peptide_ensemble/structures/cluster_0.pdb",
    "contigmap": {
      "contigs": ["A1-15/0 50-100"],
      "inpaint_chains": ["B"]
    },
    "num_designs": 100
  },
  "diffusion": {
    "partial_T": 10,
    "num_steps": 50
  }
}
```

**Explanation:**
- `A1-15`: AviTag (chain A, residues 1-15) - FIXED
- `50-100`: Design 50-100 residue binder (chain B) - DESIGN THIS
- `num_designs`: Generate 100 designs per conformation

**Create one config per cluster:**
```bash
# Directory structure
02_scaffolds/approach_C_peptide_ensemble/
  ├── structures/
  │   ├── cluster_0.pdb
  │   ├── cluster_1.pdb
  │   ├── ...
  │   └── cluster_7.pdb
  └── rfd3_configs/
      ├── cluster_0.json
      ├── cluster_1.json
      ├── ...
      └── cluster_7.json
```

---

### Step 7: Run RFDiffusion3 (Ensemble Design)

**For each conformation:**

```bash
# Run RFD3 design
cd 02_scaffolds/approach_C_peptide_ensemble

for i in {0..7}; do
  export MSYS_NO_PATHCONV=1
  docker run --rm \
    -v "$(pwd):/workspace" \
    -w /workspace \
    rfdiffusion3:latest \
    rfd3 design \
    --inputs rfd3_configs/cluster_${i}.json \
    --out_dir designs/cluster_${i} \
    --n_batches 2 \
    --diffusion_batch_size 50
done
```

**Expected output:**
- 8 clusters × 100 designs = 800 total designs
- Each cluster generates binders for that specific AviTag conformation

---

### Step 8: Filter and Validate Designs

**Step 8A: Initial filtering (pLDDT, pAE)**

```python
# Sort by pLDDT (prediction confidence)
# Keep top 10-20 designs per cluster
# Total: ~80-160 designs for validation
```

**Step 8B: ColabFold validation**

```bash
# Run ColabFold on top designs
for design in top_designs/*.pdb; do
  colabfold_batch $design output/
done
```

**Step 8C: Final selection**

Criteria:
- pLDDT > 80 (high confidence)
- pAE < 5 Å (good interface prediction)
- Visual inspection: Binder makes good contacts

**Expected:** 5-10 final candidate binders

---

## Timeline Estimate

| Step | Time | Cumulative |
|------|------|------------|
| 1. Stop MD | 5 min | 5 min |
| 2. Clustering | 30 min - 2 hours | 2 hours |
| 3. Extract structures | 30 min | 2.5 hours |
| 4. Visual QC | 30 min | 3 hours |
| 5. Cluster analysis | 30 min | 3.5 hours |
| 6. Prepare RFD3 configs | 1 hour | 4.5 hours |
| 7. Run RFD3 (8 × 100 designs) | 6-12 hours | 16.5 hours |
| 8. Filter & validate | 4-8 hours | 24 hours |

**Total:** ~1-2 days from MD stop to validated binders ✅

---

## Alternative: Quick Test Run

**If you want to test the workflow first:**

1. Extract just 1-2 cluster representatives
2. Generate 10 designs each (not 100)
3. Validate with ColabFold
4. Verify pipeline works end-to-end
5. Then scale up to full ensemble

**Time:** ~4-6 hours for complete test

---

## Success Criteria

**By the end, you should have:**

✅ 5-10 distinct AviTag conformations (from clustering)
✅ 5-10 validated binder designs (from RFD3 + ColabFold)
✅ Structural analysis of binder-peptide interfaces
✅ Ready for experimental testing or further refinement

---

## Troubleshooting

### Issue: Clustering gives too many/few clusters

**Too many (>15):** Increase `-cutoff` to 0.3 or 0.4
**Too few (<3):** Decrease `-cutoff` to 0.15 or 0.2

### Issue: RFD3 designs look weird

**Check:**
- Input PDB is clean (no broken chains)
- JSON config specifies correct residue ranges
- Contigmap syntax is correct

### Issue: ColabFold predictions fail

**Check:**
- PDB files are valid
- Sequences are reasonable length (<500 residues)
- No special characters in filenames

---

## Questions to Ask Before Proceeding

1. **Do you want to use all ~2,270 frames or subsample?**
   - Subsampling (e.g., every 10th frame) speeds up clustering
   - May miss rare conformations

2. **How many cluster representatives do you want?**
   - More = better coverage, but more RFD3 runs
   - 5-8 is a good balance

3. **How many designs per conformation?**
   - 100 is standard
   - Could do 50 for faster testing

4. **Do you want to run RFD3 with or without GPU?**
   - GPU: Much faster (6-8x)
   - CPU: Works, but slower

---

**Ready to start clustering?** Let me know and I'll help you run it! 🚀
