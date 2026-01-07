# ChimeraX + VisualFactory Installation Guide

**Goal:** Install free ChimeraX and VisualFactory to create Hollywood-quality molecular visualizations

**Date:** January 2, 2026

---

## Part 1: Install ChimeraX (Free!)

### Step 1: Download ChimeraX

**Official Website:** https://www.rbvi.ucsf.edu/chimerax/download.html

**For Windows:**
1. Go to the download page
2. Click **"Download ChimeraX 1.7.1 for Windows"** (latest stable version)
3. Download the installer: `ChimeraX-1.7.1-windows.exe` (~500 MB)

**License:** Free for academic, government, non-profit, and personal use!

### Step 2: Install ChimeraX

1. **Run the installer** (`ChimeraX-1.7.1-windows.exe`)
2. **Accept the license agreement**
3. **Choose installation location:**
   - Default: `C:\Program Files\ChimeraX`
   - Or choose custom location
4. **Create desktop shortcut:** Yes (recommended)
5. **Install**
6. **Launch ChimeraX** to verify installation

**Installation size:** ~1.5 GB

### Step 3: Verify Installation

1. Open ChimeraX
2. You should see the main interface with:
   - Graphics window (black background)
   - Command line at bottom
   - Menu bar at top

**Test command:**
```
open 1ubq
```
This should load ubiquitin (a small test protein) and display it in the window.

**If it works:** Installation successful! ✅

---

## Part 2: Install VisualFactory

### Step 1: Clone VisualFactory Repository

**Option A: Using Git (if installed)**
```bash
cd C:\Users\bmartinez\ProteinDesign
git clone https://github.com/TheVisualHub/VisualFactory.git
```

**Option B: Download ZIP (no Git needed)**
1. Go to: https://github.com/TheVisualHub/VisualFactory
2. Click green **"Code"** button
3. Click **"Download ZIP"**
4. Extract to: `C:\Users\bmartinez\ProteinDesign\VisualFactory`

### Step 2: Explore VisualFactory Modules

**Directory structure after download:**
```
VisualFactory/
  ├── FindPerspective/        # Optimal camera positioning
  ├── ImageProcessing/        # Post-render enhancement
  ├── MorphingMaster/         # Smooth transitions
  ├── UltimateSmoothMD/       # Trajectory smoothing
  ├── VideoProcessing/        # Final video composition
  └── README.md
```

### Step 3: No Additional Installation Needed!

**How VisualFactory works:**
- Scripts are written in Python
- They run inside ChimeraX's built-in Python environment
- Just **drag-and-drop** scripts into ChimeraX window
- Or run from command line

**No pip install, no conda, no dependencies!**

---

## Part 3: Test VisualFactory with ChimeraX

### Quick Test: Load and Visualize a Protein

1. **Open ChimeraX**

2. **Load a test structure:**
   ```
   open 1ubq
   ```

3. **Drag a VisualFactory script:**
   - Navigate to: `C:\Users\bmartinez\ProteinDesign\VisualFactory\FindPerspective\`
   - Find a `.py` script
   - **Drag it into ChimeraX window**
   - Script will run automatically

4. **Or run from command line:**
   ```
   runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\FindPerspective\script_name.py
   ```

---

## Part 4: Load Your MD Trajectories

### For Our AviTag MD Data

**File format conversion needed:**
- GROMACS outputs: `.xtc` (trajectory), `.gro` (topology)
- ChimeraX reads: `.xtc` directly! ✅ (with `.gro` or `.pdb` as topology)

**Loading command:**
```bash
# In ChimeraX command line:
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\npt.gro

# Then load trajectory:
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\md.xtc structureModel #1
```

**This will:**
1. Load the structure (from `npt.gro`)
2. Load all 2,270 frames of trajectory
3. Allow you to play through the simulation like a movie!

### Playback Controls

**In ChimeraX:**
- **Play/Pause:** Spacebar
- **Next frame:** Right arrow
- **Previous frame:** Left arrow
- **Jump to frame:** Type frame number in slider

---

## Part 5: Create Your First Molecular Movie

### Basic Workflow

**Step 1: Load trajectory**
```
open npt.gro
open md.xtc structureModel #1
```

**Step 2: Style the molecule**
```
# Color by flexibility (B-factor/RMSF)
color bfactor #1 palette rainbow

# Or color by chain
color bychain

# Cartoon representation
cartoon

# Hide water (if present)
hide solvent
```

**Step 3: Set camera position**
```
# Center view
view

# Rotate to nice angle
turn y 45
turn x 30

# Zoom
zoom 0.8
```

**Step 4: Record movie**
```
# Start recording
movie record

# Play through frames (or subset)
coordset #1 1,100  # Frames 1 to 100

# Wait for playback to finish, then:
movie encode output C:\Users\bmartinez\Desktop\avitag_movie.mp4 format h264 quality high
```

**Done!** You now have an MP4 movie file ready for PowerPoint.

---

## Part 6: Using VisualFactory Modules

### Module 1: UltimateSmoothMD (Smooth Trajectories)

**Purpose:** Remove jitter from MD trajectories for cinematic viewing

**Usage:**
```python
# Drag UltimateSmoothMD script into ChimeraX
# Or run:
runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\UltimateSmoothMD\smooth_trajectory.py
```

**What it does:**
- Applies temporal smoothing filter
- Removes high-frequency noise
- Makes motion look fluid and professional

### Module 2: FindPerspective (Optimal Camera Angles)

**Purpose:** Automatically find best viewing angle

**Usage:**
```python
runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\FindPerspective\find_best_angle.py
```

**What it does:**
- Tests multiple camera positions
- Scores each based on visibility of features
- Returns best angle automatically

### Module 3: MorphingMaster (Smooth Transitions)

**Purpose:** Create smooth interpolations between conformations

**Usage:**
```python
runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\MorphingMaster\morph_states.py
```

**What it does:**
- Takes two conformations (e.g., compact vs extended)
- Creates smooth transition between them
- Better than jumping between frames

### Module 4: VideoProcessing (Final Touches)

**Purpose:** Post-processing effects (lighting, shadows, depth-of-field)

**Usage:**
```python
runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\VideoProcessing\enhance_video.py
```

---

## Part 7: No Docker Needed!

### Why ChimeraX Doesn't Need Docker

**ChimeraX is a native Windows application:**
- Installs directly on Windows (like Microsoft Word)
- No virtualization needed
- Direct GPU access (better rendering performance)
- Simpler to use

**Docker is for:**
- GROMACS (MD simulation) ✅ We use Docker
- RFDiffusion3 (binder design) ✅ We use Docker

**Native Windows is for:**
- ChimeraX (visualization) ✅ Better performance
- PyMOL (alternative visualization)
- Analysis scripts (Python)

---

## Part 8: Complete Workflow for AviTag Project

### After Fusion MD Finishes

**1. Load Free Peptide Trajectory**
```
# In ChimeraX:
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\npt.gro
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\md.xtc structureModel #1
```

**2. Load Fusion Trajectory (when ready)**
```
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\03_fusion_md\npt.gro
open C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\03_fusion_md\md.xtc structureModel #2
```

**3. Set up side-by-side view**
```
# Split screen
split h

# Assign models to screens
view #1 screen left
view #2 screen right
```

**4. Style both identically**
```
# Color by flexibility
color bfactor #1,2 palette rainbow

# Cartoon representation
cartoon #1,2

# Hide water
hide #1,2 solvent
```

**5. Apply VisualFactory smoothing**
```
runscript C:\Users\bmartinez\ProteinDesign\VisualFactory\UltimateSmoothMD\smooth.py
```

**6. Record comparison movie**
```
movie record
coordset #1,2 1,end  # Play all frames
movie encode output C:\Users\bmartinez\Desktop\free_vs_fusion_comparison.mp4
```

---

## Part 9: Quick Reference Commands

### Essential ChimeraX Commands

| Task | Command |
|------|---------|
| Open structure | `open filename.pdb` |
| Open trajectory | `open traj.xtc structureModel #1` |
| Color by chain | `color bychain` |
| Color by B-factor | `color bfactor palette rainbow` |
| Cartoon view | `cartoon` |
| Hide water | `hide solvent` |
| Center view | `view` |
| Rotate | `turn y 45` |
| Record movie | `movie record` |
| Encode movie | `movie encode output file.mp4` |
| Run script | `runscript path/to/script.py` |
| Get help | `help commandname` |

### File Format Support

| Format | Extension | Purpose |
|--------|-----------|---------|
| PDB | `.pdb` | Protein structures |
| GRO | `.gro` | GROMACS structures |
| XTC | `.xtc` | GROMACS trajectories |
| DCD | `.dcd` | NAMD trajectories |
| CIF | `.cif` | Crystallographic data |
| MP4 | `.mp4` | Output movies |

---

## Part 10: Troubleshooting

### Issue: ChimeraX won't open .xtc files

**Solution:**
- Make sure you load topology first (`.gro` or `.pdb`)
- Then load trajectory with: `open file.xtc structureModel #1`

### Issue: Movie is too jittery

**Solution:**
- Use VisualFactory's UltimateSmoothMD
- Or reduce frame rate: `movie encode framerate 15`

### Issue: Can't see protein well

**Solution:**
```
# Better representation
cartoon
color bychain
lighting soft
```

### Issue: File paths with spaces

**Solution:**
- Use quotes: `open "C:\Path With Spaces\file.pdb"`
- Or use forward slashes: `open C:/Path/file.pdb`

---

## Part 11: Next Steps

### Immediate Actions:

1. **Download ChimeraX** ✅
   - Go to https://www.rbvi.ucsf.edu/chimerax/download.html
   - Download Windows installer
   - Install

2. **Clone VisualFactory** ✅
   - Download from https://github.com/TheVisualHub/VisualFactory
   - Extract to your ProteinDesign folder

3. **Test with existing trajectory** ✅
   - Load free peptide MD (23 ns, already done)
   - Play through frames
   - Export test movie

4. **Wait for fusion MD to complete** ⏳
   - NPT running now (~5-10 min)
   - Production MD next (~2 days)
   - Then: side-by-side comparison

5. **Create presentation materials** 📽️
   - Side-by-side comparison movie
   - Key frame stills for PowerPoint
   - Annotated figures

---

## Summary

**What you're getting:**
- **ChimeraX:** Free, professional molecular visualization software
- **VisualFactory:** Hollywood-quality rendering scripts
- **Your data:** MD trajectories of AviTag (free + fusion)

**What you'll create:**
- Side-by-side comparison movies
- High-quality stills for presentations
- Professional visualizations for papers

**Cost:** $0 (everything is free!)

**Time investment:**
- Installation: 30 minutes
- Learning basics: 1-2 hours
- Creating first movie: 30 minutes
- Professional-quality output: Priceless!

---

**Ready to install? Let me know if you hit any issues!** 🚀

**Resources:**
- ChimeraX documentation: https://www.rbvi.ucsf.edu/chimerax/docs/user/index.html
- VisualFactory repo: https://github.com/TheVisualHub/VisualFactory
- ChimeraX tutorials: https://www.rbvi.ucsf.edu/chimerax/docs/user/tutorials.html
