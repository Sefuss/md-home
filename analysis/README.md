# Analysis Scripts (Python)

**Purpose:** Analyze MD trajectories and generate publication-quality plots

**Software Required:** Python with MDAnalysis, matplotlib, numpy

---

## Scripts

### analyze_trajectory.py
**Use Case:** Complete statistical analysis of MD trajectory

**What it analyzes:**
- RMSD (Root Mean Square Deviation)
- Radius of Gyration (Rg) - compactness measure
- State transitions (compact ⟷ extended)
- Temperature stability
- Energy conservation

**How to use:**

```bash
# Navigate to MD preparation folder
cd "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation"

# Run analysis
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" ../analysis/analyze_trajectory.py
```

**Output:**
- Terminal: Statistical summary
- Console output with key metrics

**Example output:**
```
Loaded: 2396 frames, 23.95 ns
RMSD range: 0.2 - 25.0 Å (avg: 8.6 Å)
Rg range: 7.7 - 18.5 Å (avg: 12.8 Å)
Transitions: 22 compact⟷extended switches
Conclusion: Dynamic equilibrium (NOT progressive unfolding)
```

**Files it reads:**
- `npt.gro` - Structure file
- `md.xtc` - Trajectory file

---

### create_powerpoint_plots.py
**Use Case:** Generate 4 publication-quality plots for presentations

**What it creates:**
1. **MD_Analysis_4Panel.png** - Complete overview (RMSD + Rg + distributions)
2. **MD_TimeSeries_States.png** - Time evolution with state regions marked
3. **MD_StateSpace_2D.png** - 2D conformational landscape
4. **MD_Summary_Simple.png** - Dual Y-axis for non-experts

**How to use:**

```bash
# Navigate to MD preparation folder
cd "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation"

# Generate plots
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" ../analysis/create_powerpoint_plots.py
```

**Output files:**
- All in same directory as script
- 300 DPI resolution (publication quality)
- PNG format (PowerPoint-compatible)

**Plot specifications:**
- **MD_Analysis_4Panel.png** (708 KB)
  - 2x2 grid layout
  - RMSD time series (top left)
  - Rg time series (top right)
  - RMSD distribution (bottom left)
  - Rg distribution (bottom right)

- **MD_TimeSeries_States.png** (1.8 MB)
  - RMSD colored by Rg (top)
  - Rg with compact/extended regions marked (bottom)
  - Shows state transitions clearly

- **MD_StateSpace_2D.png** (1.5 MB)
  - RMSD vs Rg scatter plot
  - Colored by time (progression)
  - Start/end markers
  - Shows conformational exploration

- **MD_Summary_Simple.png** (560 KB)
  - RMSD and Rg on same plot (dual Y-axis)
  - Simplified for non-expert audiences
  - Annotation box with summary stats

---

## Requirements

**Python packages:**
```
MDAnalysis
matplotlib
numpy
scipy
```

**Install if needed:**
```bash
conda install -c conda-forge mdanalysis matplotlib numpy scipy
```

---

## Typical Workflow

**Step 1: Run statistical analysis**
```bash
cd 01_md_preparation
python ../analysis/analyze_trajectory.py
```

**Step 2: Generate plots**
```bash
python ../analysis/create_powerpoint_plots.py
```

**Step 3: Use in presentation**
- Insert PNG files into PowerPoint
- High resolution ensures clarity
- Use MD_Summary_Simple.png for general audiences
- Use 4-panel for technical talks

---

## Understanding the Results

### Free AviTag Peptide (23 ns)

**RMSD: 8.6 Å average**
- HIGH flexibility
- This is GOOD for intrinsically disordered peptide (IDP)
- Shows dynamic exploration of conformational space

**Radius of Gyration: 7.7 - 18.5 Å**
- Wide range = transitions between compact and extended states
- NOT progressive unfolding (would be one-way increase)

**State Transitions: 22 switches**
- Dynamic equilibrium
- Peptide samples both compact and extended conformations
- Validates ensemble approach for binder design

---

## Adapting for Fusion Protein

**When fusion MD completes (~16 days):**

Edit file paths in scripts to point to fusion trajectory:

```python
# In analyze_trajectory.py and create_powerpoint_plots.py
# Change:
structure_file = '../01_md_preparation/npt.gro'
trajectory_file = '../01_md_preparation/md.xtc'

# To:
structure_file = '../03_fusion_md/npt.gro'
trajectory_file = '../03_fusion_md/md.xtc'
```

**Then extract AviTag residues only:**
```python
# Select only AviTag residues (135-149 in fusion)
avitag = u.select_atoms('resid 135-149')
```

**Compare results:**
- Free peptide RMSD vs Fusion-bound RMSD
- Free peptide Rg vs Fusion-bound Rg
- Make decision on which ensemble to use

---

## Troubleshooting

**Error: ModuleNotFoundError: No module named 'MDAnalysis'**

**Solution:**
```bash
conda install -c conda-forge mdanalysis
```

**Error: Cannot find npt.gro or md.xtc**

**Solution:**
- Make sure you're in `01_md_preparation/` directory
- Or adjust file paths in script

**Error: UnicodeEncodeError**

**Solution:**
- Set environment variable before running:
```bash
export PYTHONIOENCODING=utf-8
python script.py
```

**Warning: Trajectory has X frames but structure has Y atoms**

**This is normal** - Just means topology and trajectory are compatible

---

## Future Analysis Ideas

**Clustering:**
- Use `MDAnalysis.analysis.encore` for clustering
- Extract representative structures
- Use for RFDiffusion3 input

**Secondary structure:**
- Use `MDAnalysis.analysis.dssp`
- Track helix/sheet content over time

**Contact maps:**
- Calculate residue-residue distances
- Identify transient interactions

**Principal Component Analysis:**
- Identify dominant motions
- Project trajectory onto PC space

---

## Further Reading

- **MDAnalysis docs:** https://docs.mdanalysis.org/
- **Matplotlib gallery:** https://matplotlib.org/stable/gallery/
- **Project guide:** `../MD_ANALYSIS_REPORT.md`

---

**Last Updated:** January 2, 2026
