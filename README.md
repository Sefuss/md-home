# AviTag Binder Design: Peptide Ensemble Approach

**Project:** AviTag-Specific Binder Design via Conformational Ensemble Sampling
**Start Date:** December 17, 2025
**Updated:** January 5, 2026
**Status:** Free peptide MD complete (23 ns), Fusion MD running (100 ns)

## Quick Navigation

🏠 **[REPRODUCING_AT_HOME.md](REPRODUCING_AT_HOME.md)** - ⭐ Reproduce MD work at home using checkpoints
📖 **[SETUP_GUIDE.md](SETUP_GUIDE.md)** - How to set up GROMACS, ChimeraX, and run MD simulations
🎬 **[ANALYSIS_GUIDE.md](ANALYSIS_GUIDE.md)** - How to analyze trajectories and create movies
🔧 **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** - Solutions to common issues
📚 **[docs/](docs/)** - Technical references and archived reports

---

## Project Overview

This approach aims to design high-affinity binders to the AviTag peptide (sequence: GLNDIFEAQKIEWHE) by first sampling its conformational ensemble through molecular dynamics simulation, then using representative structures for RFDiffusion3-based binder design.

### Experimental Design

**Phase 1: Free Peptide Ensemble (COMPLETE)**
- 23 ns MD simulation of isolated AviTag in solution
- Result: Dynamic ensemble with 22 compact⟷extended transitions
- RMSD: 8.6 Å average (high flexibility expected for intrinsically disordered peptide)
- Status: COMPLETE ✅

**Phase 2: Fusion Protein Context (IN PROGRESS)**
- 100 ns MD simulation of LaG16-AviTag fusion protein
- Purpose: Test if protein context changes peptide behavior
- ETA: ~16 days (started Jan 2, 2026)
- Status: RUNNING ⏳

**Phase 3: Comparison and Decision (PENDING)**
- Compare free vs fusion AviTag conformational ensembles
- Decide which ensemble to use for binder design
- Decision matrix in: `FUSION_PROTEIN_MD_PLAN.md`

**Phase 4: Binder Design (PENDING)**
- Cluster selected trajectory
- Extract 5-10 representative conformations
- Run RFDiffusion3 with partial diffusion on peptide conformations
- Validate designs with AlphaFold2/ColabFold

---

## Directory Structure

```
approach_C_peptide_ensemble/
├── README.md                           # This file - Master index
│
├── 00_initial_structure/               # Starting structures
│   ├── avitag.fasta                    # Sequence
│   ├── avitag_initial.pdb              # ColabFold prediction
│   └── build_avitag_peptide.py         # Structure generation script
│
├── 01_md_preparation/                  # Free peptide MD
│   ├── npt.gro                         # NPT equilibrated structure
│   ├── md.xtc                          # Raw trajectory (23 ns, 2396 frames)
│   ├── md_centered.xtc                 # CLEANED trajectory (PBC-corrected) ✅
│   └── md.tpr                          # Run input file
│
├── 03_fusion_md/                       # Fusion protein MD (LaG16-AviTag)
│   ├── fusion_processed.gro            # Starting structure (149 residues)
│   ├── md.xtc                          # Trajectory (IN PROGRESS)
│   └── ... (equilibration files)
│
├── visualization/                      # ChimeraX scripts (5 scripts + README)
│   ├── README.md                       # Visualization guide
│   ├── create_avitag_movie.py          # ⭐⭐⭐ ONE-CLICK movie maker
│   ├── chimerax_load_free_peptide_CLEAN.py  # ⭐ Load PBC-corrected
│   ├── chimerax_load_free_peptide.py   # Quick viewing
│   ├── chimerax_export_movie.py        # Manual export
│   └── chimerax_side_by_side_comparison.py  # Free vs fusion
│
├── analysis/                           # Python analysis scripts (2 scripts + README)
│   ├── README.md                       # Analysis guide
│   ├── analyze_trajectory.py           # Statistical analysis
│   └── create_powerpoint_plots.py      # Generate plots
│
├── utilities/                          # Helper scripts (1 script + README)
│   ├── README.md                       # Utilities guide
│   └── monitor_md.py                   # Real-time MD monitoring
│
├── scripts/                            # Legacy (will be deprecated)
│   └── visualize_benchmark.py          # Performance visualization
│
└── Documentation (13 markdown files organized below)
```

---

## Documentation Index

### 📚 COMPREHENSIVE GUIDES (Start Here!)

**1. HOLLYWOOD_QUALITY_MOVIE_GUIDE.md** ⭐ NEW!
   - Complete guide to creating cinematic molecular movies
   - Based on TheVisualHub's VisualFactory philosophy
   - Covers: PBC artifacts, when to show solvent, quality settings
   - Includes: Gleb's smoothing techniques, camera work, color schemes
   - Your complete reference for presentation movies!

**2. CHIMERAX_SETUP_GUIDE.md**
   - Installation guide for ChimeraX and VisualFactory
   - Step-by-step setup (no Docker needed)
   - Usage examples and command reference
   - Links to official resources

**3. MD_SIMULATION_EXPLAINED.md**
   - Ground-up explanation of entire MD workflow
   - Docker containerization explained
   - GROMACS 12-step process
   - Physics being calculated (force fields, equations)
   - Perfect for learning/teaching MD simulations

**4. VISUALIZING_MD_SIMULATION.md**
   - Beginner-friendly visualization guide
   - How to open and view trajectories
   - PyMOL and ChimeraX basics
   - No assumptions about prior knowledge

---

### 📊 ANALYSIS REPORTS

**5. MD_ANALYSIS_REPORT.md**
   - Complete quality assessment of free peptide MD
   - Temperature stability: 300.022 K ✅
   - Energy conservation: 0.87 kJ/mol/ns drift ✅
   - RMSD analysis: 8.6 Å average (good for IDP)
   - Artifact check: None detected
   - Conclusion: Scientifically valid, ready for clustering

**6. MD_MONITORING_REPORT_FINAL.md**
   - Executive summary and recommendation
   - Progress: 23 ns / 100 ns target
   - Performance: 13.6 ns/day
   - Recommendation: Stop at 23 ns (diminishing returns)
   - Decision: Proceed with clustering

---

### 🎯 PROJECT PLANNING

**7. FUSION_PROTEIN_MD_PLAN.md**
   - Critical design question: Free vs bound conformations
   - Experimental design for fusion MD
   - Comparison metrics (RMSD, Rg distributions)
   - Decision matrix for choosing ensemble
   - Addresses: Conformational selection vs induced fit

**8. NEXT_STEPS.md**
   - Action items after free peptide MD completion
   - Clustering strategy
   - RFDiffusion3 setup
   - Timeline and milestones

**9. PRODUCTION_MD_READY.md**
   - Pre-launch checklist for production MD
   - System specifications
   - Parameter validation
   - Launch instructions

---

### 🔧 TECHNICAL SETUP DOCS

**10. MD_SETUP_PLAN.md**
    - Initial planning document for free peptide MD
    - Force field selection (AMBER99SB-ILDN)
    - System size calculations
    - Expected runtime estimates

**11. MD_SETUP_PROGRESS.md**
    - Step-by-step setup tracking
    - Topology generation → Energy minimization → Equilibration
    - Validation at each step

**12. GROMACS_PERFORMANCE_METRICS.md**
    - Benchmarking results
    - Performance optimization
    - Hardware utilization

**13. MD_MONITORING_REPORT.md**
    - Early-stage monitoring report
    - Predecessor to MD_MONITORING_REPORT_FINAL.md

---

## Scripts Index

All scripts are now organized into three subdirectories, each with its own README:

### 🎬 visualization/ - ChimeraX Scripts (5 scripts)

**See: visualization/README.md for complete guide**

**create_avitag_movie.py** ⭐⭐⭐ ONE-CLICK MOVIE MAKER
   - Automated professional movie creation
   - Based on VisualFactory principles
   - Output: 1920x1080 MP4, PowerPoint-ready
   - Just drag into ChimeraX and wait!

**chimerax_load_free_peptide_CLEAN.py** ⭐ RECOMMENDED
   - Load PBC-corrected trajectory (no wrapping artifacts!)
   - Professional styling pre-configured
   - Better for presentation movies

**chimerax_load_free_peptide.py**
   - Load original free peptide trajectory
   - Basic styling and controls
   - Use for quick viewing

**chimerax_export_movie.py**
   - Manual movie export with custom settings
   - For advanced users who want full control
   - Adjustable resolution, framerate, quality

**chimerax_side_by_side_comparison.py** (For later use)
   - Load both free and fusion trajectories
   - Split-screen synchronized playback
   - Use after fusion MD completes

---

### 📈 analysis/ - Python Analysis Scripts (2 scripts)

**See: analysis/README.md for complete guide**

**analyze_trajectory.py**
   - Complete trajectory analysis
   - RMSD, Radius of Gyration, state transitions
   - Requires: MDAnalysis, matplotlib
   - Run: `python ../analysis/analyze_trajectory.py`

**create_powerpoint_plots.py**
   - Generate 4 publication-quality plots
   - 300 DPI PNG format
   - PowerPoint-ready
   - Outputs:
     - MD_Analysis_4Panel.png
     - MD_TimeSeries_States.png
     - MD_StateSpace_2D.png
     - MD_Summary_Simple.png

---

### 🛠️ utilities/ - Helper Scripts (1 script)

**See: utilities/README.md for complete guide**

**monitor_md.py**
   - Real-time monitoring of running MD simulation
   - Tracks progress, temperature, energy
   - Use while simulation is running

---

## Quick Start Guides

### How to: View MD Trajectory in ChimeraX

**Step 1:** Open ChimeraX application

**Step 2:** Drag this file into ChimeraX window:
```
chimerax_load_free_peptide_CLEAN.py
```

**Step 3:** Press Spacebar to play trajectory

**Result:** You'll see AviTag peptide dynamically exploring conformations (no wrapping artifacts!)

---

### How to: Create Presentation Movie

**Option A: One-Click Automated (Recommended)**

1. Open ChimeraX
2. Drag: `create_avitag_movie.py`
3. Wait ~3-5 minutes
4. Find movie on Desktop: `AviTag_Dynamic_Ensemble.mp4`

**Option B: Manual with Custom Settings**

1. Open ChimeraX
2. Load cleaned trajectory: Drag `chimerax_load_free_peptide_CLEAN.py`
3. (Optional) Apply smoothing: Drag `VisualFactory/UltimateSmoothMD/UltimateSmoothMD5.py`
4. Record movie: Drag `chimerax_export_movie.py`

**See:** `HOLLYWOOD_QUALITY_MOVIE_GUIDE.md` for complete tutorial

---

### How to: Analyze MD Trajectory

**Step 1:** Navigate to analysis folder
```bash
cd "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation"
```

**Step 2:** Run analysis script
```bash
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" analyze_trajectory.py
```

**Step 3:** Generate PowerPoint plots
```bash
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" create_powerpoint_plots.py
```

**Outputs:**
- Terminal: Statistical summary
- PNG files: 4 high-resolution plots

---

### How to: Monitor Fusion MD Progress

**Check if running:**
```bash
docker ps -a --filter "name=fusion_md_production"
```

**View latest output:**
```bash
docker logs --tail 50 fusion_md_production
```

**Check progress continuously:**
```bash
docker logs -f fusion_md_production
```

(Press Ctrl+C to stop viewing)

---

## Key Results Summary

### Free Peptide MD (23 ns)

**Simulation Quality:**
- Temperature: 300.022 K ± 2.8 K ✅
- Energy drift: 0.87 kJ/mol/ns ✅
- No artifacts detected ✅

**Conformational Dynamics:**
- RMSD: 0.2 - 25 Å range (average 8.6 Å)
- Radius of Gyration: 7.7 - 18.5 Å
- State transitions: 22 compact⟷extended switches
- Conclusion: Dynamic equilibrium (NOT progressive unfolding)

**Scientific Interpretation:**
- AviTag is intrinsically disordered peptide (IDP)
- High flexibility is expected and desirable
- Explores diverse conformational space
- Suitable for ensemble-based binder design

**Visualization Outputs:**
- 4 PowerPoint plots (300 DPI)
- Cleaned trajectory for ChimeraX (PBC-corrected)
- Ready for movie creation

---

## Troubleshooting

### Issue: Trajectory shows wrapping artifacts in ChimeraX

**Solution:**
Use the CLEANED trajectory instead:
- File: `01_md_preparation/md_centered.xtc`
- Script: `chimerax_load_free_peptide_CLEAN.py`

This trajectory was pre-processed with GROMACS to remove periodic boundary artifacts.

### Issue: Python script not running

**Solution:**
Use full Anaconda Python path:
```bash
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" script_name.py
```

### Issue: Docker container not found

**Check if running:**
```bash
docker ps -a
```

**Restart if stopped:**
```bash
docker start fusion_md_production
```

### Issue: ChimeraX movie too jittery

**Solution 1:** Apply VisualFactory smoothing
```
Drag: VisualFactory/UltimateSmoothMD/UltimateSmoothMD5.py
```

**Solution 2:** Increase supersampling
```python
movie record supersample 4  # Instead of 3
```

**Solution 3:** Play subset of frames
```python
coordset #1 1,end,2  # Every 2nd frame
```

---

## External Resources

### Software Documentation
- **ChimeraX:** https://www.rbvi.ucsf.edu/chimerax/docs/user/index.html
- **GROMACS:** https://manual.gromacs.org/
- **MDAnalysis:** https://docs.mdanalysis.org/
- **RFDiffusion3:** https://github.com/baker-laboratory/RoseTTAFold-All-Atom

### VisualFactory (TheVisualHub)
- **Repository:** https://github.com/TheVisualHub/VisualFactory
- **YouTube:** https://www.youtube.com/@TheVisualHub
- **Philosophy:** "Automated 🔮 Precise 🎯 Beautiful 🌺"

### Project-Specific
- **ChimeraX Setup:** See `CHIMERAX_SETUP_GUIDE.md`
- **Movie Making:** See `HOLLYWOOD_QUALITY_MOVIE_GUIDE.md`
- **MD Explanation:** See `MD_SIMULATION_EXPLAINED.md`

---

## Current Status (January 2, 2026)

### ✅ Completed

1. Free AviTag peptide MD (23 ns)
2. Quality analysis and validation
3. PowerPoint visualization plots
4. PBC-corrected trajectory for movies
5. ChimeraX + VisualFactory installation
6. Comprehensive movie-making guide
7. Fusion protein MD setup

### ⏳ In Progress

1. Fusion protein MD production run
   - Container: `fusion_md_production`
   - Target: 100 ns
   - ETA: ~16 days
   - Status: Running in background

### 📋 Next Steps

1. Test cleaned trajectory in ChimeraX
2. Create professional presentation movie
3. Wait for fusion MD completion (~16 days)
4. Compare free vs fusion ensembles
5. Decide on ensemble approach
6. Cluster trajectory
7. Extract representatives
8. RFDiffusion3 binder design

---

## File Format Reference

### MD Trajectory Files

| Extension | Description | Software |
|-----------|-------------|----------|
| `.xtc` | Compressed trajectory | GROMACS |
| `.gro` | Structure file | GROMACS |
| `.tpr` | Run input file | GROMACS |
| `.edr` | Energy file | GROMACS |
| `.log` | Log file | GROMACS |

### Visualization Files

| Extension | Description | Software |
|-----------|-------------|----------|
| `.pdb` | Protein structure | Universal |
| `.mp4` | Movie file | ChimeraX output |
| `.png` | Image file | Analysis plots |
| `.py` | Python script | ChimeraX/Analysis |

---

## Command Quick Reference

### ChimeraX Commands

```python
# Loading
open structure.gro
open trajectory.xtc structureModel #1

# Styling
cartoon                  # Cartoon representation
hide solvent            # Hide water
color bfactor #1 palette rainbow  # Color by flexibility
lighting soft           # Soft lighting
graphics silhouettes true  # Edge outlines

# Camera
view                    # Center molecule
turn y 45               # Rotate
zoom 0.9                # Zoom out

# Recording
windowsize 1920 1080    # Set HD resolution
movie record supersample 3  # Start recording
coordset #1 1,end       # Play all frames
wait 100                # Wait for playback
movie encode file.mp4   # Save movie
```

### GROMACS Commands

```bash
# Fix PBC wrapping
echo "1 0" | gmx trjconv -s md.tpr -f md.xtc -o md_centered.xtc -center -pbc mol -ur compact

# Extract energy
gmx energy -f md.edr -o temperature.xvg

# Calculate RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg

# Calculate Rg
gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
```

### Docker Commands

```bash
# Check running containers
docker ps

# View logs
docker logs fusion_md_production

# Follow logs
docker logs -f fusion_md_production

# Stop container
docker stop fusion_md_production

# Start container
docker start fusion_md_production
```

---

## Project Timeline

| Date | Milestone |
|------|-----------|
| Dec 17, 2025 | Project started, free peptide MD launched |
| Dec 18-31, 2025 | Free peptide MD running (holiday gap) |
| Jan 2, 2026 | Returned, quality assessment, stopped at 23 ns |
| Jan 2, 2026 | ChimeraX + VisualFactory installed |
| Jan 2, 2026 | Fusion protein MD launched |
| Jan 18, 2026 (Est.) | Fusion MD completion |
| TBD | Ensemble comparison and decision |
| TBD | RFDiffusion3 binder design |
| TBD | Experimental validation |

---

## Contact & Support

**Project Lead:** Boe Martinez
**Date Range:** December 2025 - January 2026
**Platform:** Windows with Docker
**Tools:** GROMACS, ChimeraX, VisualFactory, MDAnalysis, Python

**For Issues:**
1. Check troubleshooting section above
2. Review relevant guide (13 markdown files)
3. Check specific script documentation

---

**Last Updated:** January 2, 2026
**Version:** 2.0 (Added VisualFactory integration and movie-making guides)
