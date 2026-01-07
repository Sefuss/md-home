# MD Ensemble Setup Progress
**Date:** 2025-12-17
**Status:** System preparation complete, ready for minimization

---

## Completed Steps

### 1. Structure Generation ✅
**Method:** ColabFold (LocalColabFold via Docker)
- **Input:** AviTag sequence (GLNDIFEAQKIEWHE, 15 residues)
- **Model:** AlphaFold2 PTM
- **Runtime:** 3m 44s total
  - Weight download: ~3 minutes (one-time, now cached)
  - Prediction: 21.8 seconds (3 recycles)
  - AMBER relaxation: 2.4 seconds
- **Quality:**
  - pLDDT: 72.7 (moderate, typical for flexible peptides)
  - pTM: 0.037 (expected for monomeric peptide)
- **Output:** `AviTag_peptide_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb`
  - All-atom structure with proper histidine ring
  - 130 heavy atoms + hydrogens

### 2. GROMACS Topology Generation ✅
**Tool:** `gmx pdb2gmx`
- **Force Field:** AMBER99SB-ILDN
- **Water Model:** TIP3P
- **Input:** 130 atoms (ColabFold structure)
- **Output:** 248 atoms (hydrogens added)
- **Histidine Protonation:** HISE (epsilon protonated, auto-assigned)
- **Generated Files:**
  - `avitag_processed.gro` - Coordinates with hydrogens
  - `topol.top` - Molecular topology
- **System Charge:** -3 e (before neutralization)

### 3. System Solvation ✅
**Tool:** `gmx editconf` + `gmx solvate`

**Box Creation:**
- Type: Cubic
- Size: 4.824 × 4.824 × 4.824 nm
- Volume: 112.26 nm³
- Padding: 1.0 nm from peptide to box edge

**Solvation:**
- Water molecules added: 3,636 (SOL)
- Water atoms: 10,908
- Density: 997.6 g/L (proper liquid water)
- Output: `avitag_solvated.gro`

### 4. System Neutralization ✅
**Tool:** `gmx genion`
- **Ions Added:** 3 Na+ (sodium)
- **Replaced:** 3 water molecules
- **Final Charge:** 0 (neutral system)
- **Output:** `avitag_ions.gro`

---

## Final System Composition

| Component | Count | Atoms | Notes |
|-----------|-------|-------|-------|
| AviTag peptide | 1 | 248 | 15 residues, all atoms + H |
| Na+ ions | 3 | 3 | Neutralizing ions |
| Water (SOL) | 3,633 | 10,899 | TIP3P model |
| **Total** | **3,637** | **11,150** | Neutral, physiological |

**System Properties:**
- Box: 4.824 nm cubic
- Volume: 112.26 nm³
- Net charge: 0
- Temperature: Will equilibrate to 300 K
- Pressure: Will equilibrate to 1 bar

---

## Issues Resolved

### Histidine Ring Error (SOLVED)
**Problem:** Initial BioPython builder created backbone-only structure. GROMACS failed on histidine (residue 14) with "Incomplete ring" error.

**Root Cause:** Histidine requires sidechain atoms (CB, CG, ND1, CD2, CE1, NE2) for protonation state determination.

**Solution:** Used ColabFold to generate proper all-atom structure with complete histidine ring.

**Why ColabFold?**
- Generates all-atom structures (backbone + sidechains)
- AMBER relaxation ensures proper geometry
- Already installed locally (no TPU needed for 15-mer peptide)
- Fast for small peptides (~22s prediction time)

---

## Next Steps

### 5. Energy Minimization (COMPLETED)
**Purpose:** Relax system, remove steric clashes
- **Algorithm:** Steepest descent
- **Convergence:** Max force < 1000 kJ/mol/nm
- **Actual Time:** ~3 minutes (CPU)
- **Result:** Converged at step 373, final energy -172,562 kJ/mol

### 6. NVT Equilibration (COMPLETED)
**Purpose:** Temperature equilibration (constant volume)
- **Duration:** 100 ps
- **Target:** 300 K (physiological)
- **Thermostat:** V-rescale
- **Actual Time:** 5.8 minutes
- **Performance:** 24.743 ns/day

### 7. NPT Equilibration (COMPLETED)
**Purpose:** Pressure equilibration (constant pressure)
- **Duration:** 100 ps
- **Target:** 1 bar
- **Barostat:** Parrinello-Rahman
- **Actual Time:** 5.4 minutes
- **Performance:** 26.636 ns/day

### 8. Production MD (READY TO START)
**Purpose:** Sample conformational ensemble
- **Duration:** 100 ns
- **Save Frequency:** Every 10 ps (~10,000 frames)
- **Expected Time:** ~4 days continuous (CPU, i3-10100)
- **Output Size:** ~573 MB compressed trajectory
- **Checkpoint Frequency:** Every 100 ps (0.1 ns)

### 9. Trajectory Analysis (PENDING)
**Tools:** MDAnalysis + scikit-learn
- Calculate RMSD over time
- Cluster conformations (K-means or DBSCAN)
- Identify 1-3 dominant states
- Extract representative structures

### 10. RFD3 Binder Design (PENDING)
**Input:** Cluster representative structures
- Design binders for top 1-3 conformations
- Compare designs across states
- Select best candidates for validation

---

## File Structure

```
02_scaffolds/approach_C_peptide_ensemble/
├── 00_initial_structure/
│   ├── avitag.fasta                    # Input sequence
│   ├── colabfold_output/               # ColabFold results
│   │   ├── AviTag_peptide_relaxed_rank_001_*.pdb  # Main structure
│   │   ├── AviTag_peptide_plddt.png    # Confidence plot
│   │   └── ...
│   └── build_avitag_peptide.py         # Original builder (not used)
├── 01_md_preparation/
│   ├── avitag_processed.gro            # After pdb2gmx
│   ├── avitag_boxed.gro                # After box creation
│   ├── avitag_solvated.gro             # After solvation
│   ├── avitag_ions.gro                 # After neutralization (CURRENT)
│   ├── topol.top                       # System topology
│   └── ions.mdp                        # Ion addition parameters
├── 02_md_production/                   # (To be created)
├── 03_trajectory_analysis/             # (To be created)
├── 04_cluster_representatives/         # (To be created)
└── 05_rfd3_configs/                    # (To be created)
```

---

## Software Used

| Tool | Version | Purpose |
|------|---------|---------|
| ColabFold | 1.5.5 | Structure prediction |
| GROMACS | 2022.2 | MD simulation |
| AMBER99SB-ILDN | - | Force field |
| TIP3P | - | Water model |
| Docker | - | Containerization |

---

## Performance Notes

**Structure Generation (ColabFold):**
- First run: ~3m 44s (includes weight download)
- Future runs: ~30-40s (weights cached)
- CPU-only mode (no GPU in WSL)

**System Setup (GROMACS):**
- pdb2gmx: <1 second
- Solvation: <5 seconds
- Ion addition: <5 seconds

**Expected MD Timings:**
- Energy minimization: 5-10 minutes
- NVT equilibration: 10-15 minutes
- NPT equilibration: 10-15 minutes
- Production MD (100 ns): 12-24 hours

**Hardware:** Intel i3-10100 (4 cores, 8 threads), 32 GB RAM

---

## Key Decisions

1. **ColabFold over manual building:** Proper all-atom structures
2. **AMBER99SB-ILDN force field:** Well-validated for peptides
3. **TIP3P water:** Standard, fast, compatible with AMBER
4. **100 ns initial run:** Balance between sampling and time
5. **1.0 nm padding:** Small but sufficient for peptide

---

**Last Updated:** 2025-12-17 10:30 AM
**Status:** Ready for energy minimization
