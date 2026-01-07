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
# How to Visualize Your MD Simulation
**Watch water molecules, peptide movement, temperature/pressure changes in real-time!**

---

## Option 1: PyMOL (You Have This!)

### Load and Play Trajectory

**Basic Loading:**
```python
# In PyMOL command line or script:
load 01_md_preparation/npt.gro, structure
load_traj 01_md_preparation/md.xtc, structure

# Play movie
mplay

# Control speed
mset 1 x500  # Show first 500 frames
mplay
```

**Enhanced Visualization:**
```python
# Load structure and trajectory
load 01_md_preparation/npt.gro, peptide
load_traj 01_md_preparation/md.xtc, peptide

# Visual styling
hide everything
show cartoon, peptide and polymer
show spheres, peptide and resn SOL  # Show water
set sphere_scale, 0.3, resn SOL    # Small water spheres

# Color by element
color cyan, elem C
color red, elem O
color blue, elem N
color white, elem H

# Center on peptide
zoom peptide and polymer
orient peptide and polymer

# Play movie
mplay
```

**Real-Time Loading (Watch as MD Runs):**
PyMOL can load XTC files that are being written. Just reload periodically:
```python
load 01_md_preparation/npt.gro, sim
load_traj 01_md_preparation/md.xtc, sim

# Wait a few minutes, then reload trajectory
load_traj 01_md_preparation/md.xtc, sim, state=1
```

### Create Movie/Video

```python
# Load trajectory
load 01_md_preparation/npt.gro, movie
load_traj 01_md_preparation/md.xtc, movie

# Set up view
hide everything
show cartoon, polymer
show sticks, polymer and sidechain
set ray_shadows, 0
set antialias, 2

# Render movie (THIS TAKES TIME!)
set ray_trace_frames, 1
mpng frame_, width=1920, height=1080

# Convert PNGs to video with ffmpeg:
# ffmpeg -framerate 30 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p md_movie.mp4
```

---

## Option 2: VMD (Free, Excellent for MD)

### Install VMD
Download from: https://www.ks.uiuc.edu/Research/vmd/

### Load Trajectory

**Interactive GUI:**
1. File → New Molecule → Load GRO file (npt.gro)
2. File → Load Data Into Molecule → Load XTC file (md.xtc)
3. Graphics → Representations
4. Click "Create Rep" to add visualization styles

**Command Line (Tcl):**
```tcl
mol new npt.gro type gro waitfor all
mol addfile md.xtc type xtc waitfor all

# Representations
mol representation NewCartoon
mol selection "protein"
mol addrep 0

mol representation VDW 0.3
mol selection "water"
mol addrep 0

# Play animation
animate goto start
animate forward
```

### Real-Time Monitoring

VMD has plugins for watching MD as it runs:

```tcl
# Load initial structure
mol new npt.gro
mol addfile md.xtc waitfor all

# Auto-reload trajectory every 60 seconds
proc reload_traj {} {
    mol delete all
    mol new npt.gro
    mol addfile md.xtc waitfor all
    after 60000 reload_traj  # Repeat after 60 seconds
}

reload_traj
```

### Beautiful Rendering (Tachyon)

```tcl
# Set up nice view
display projection orthographic
axes location off
color Display Background white

mol representation NewCartoon
mol selection "protein"
mol color Structure
mol addrep 0

# Render with Tachyon
render Tachyon md_frame.tga
```

---

## Option 3: ChimeraX (Recommended for Movies)

### Install ChimeraX
Download from: https://www.rbvi.ucsf.edu/chimerax/download.html
(Free for non-commercial use - you qualify!)

### Load and Visualize

**Command Line:**
```python
# Open structure and trajectory
open npt.gro
open md.xtc structureModel #1

# Visual styling
preset ghost  # Nice transparent style
surface protein
color protein bychain
set bgColor white

# Play movie
coordset #1 1,last,1
wait

# Record movie
movie record size 3840,2160
coordset #1 1,last,1; wait
movie encode format h264 quality high output md_movie.mp4
```

### Using TheVisualHub Scripts

Remember the VisualFactory scripts we researched? Use them!

1. Clone: `git clone https://github.com/TheVisualHub/VisualFactory.git`
2. Open ChimeraX
3. Run → Script → Select `UltimateSmoothMD5.py`
4. Follow prompts to smooth trajectory
5. Use `FindPerspective.py` for automatic camera angles
6. Export 4K movie!

---

## Option 4: MDAnalysis (Python - For Quantitative Analysis)

### Watch Specific Atoms/Properties

```python
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load trajectory
u = mda.Universe('npt.gro', 'md.xtc')

# Select groups
protein = u.select_atoms('protein')
water = u.select_atoms('resname SOL')

# Calculate RMSD over time
from MDAnalysis.analysis import rms
R = rms.RMSD(protein, protein, select='backbone')
R.run()

# Plot
plt.plot(R.results.time, R.results.rmsd)
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.title('Backbone RMSD Over Time')
plt.show()

# Animate in 3D
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(frame):
    u.trajectory[frame]
    ax.clear()
    pos = protein.positions
    ax.scatter(pos[:,0], pos[:,1], pos[:,2], c='blue', s=10)
    ax.set_title(f'Frame {frame}: {u.trajectory.time:.1f} ps')

ani = FuncAnimation(fig, update, frames=len(u.trajectory), interval=50)
plt.show()
```

---

## What to Watch For

### 1. **Water Molecules**
- Should look like they're "boiling" gently
- Constant motion around peptide
- No systematic drift (indicates stable box)

**How to see:**
- VMD: `mol representation VDW 0.2; mol selection "water"; mol addrep 0`
- PyMOL: `show spheres, resn SOL; set sphere_scale, 0.2`

### 2. **Peptide Movement**
- Backbone should be relatively stable
- Sidechains should wiggle freely
- Whole peptide may rotate/translate slowly

**How to see:**
- VMD: `mol representation NewCartoon; mol selection "protein"`
- PyMOL: `show cartoon, polymer`

### 3. **Temperature/Pressure (Indirect)**
- Faster motion = higher temperature
- Box size changes = pressure equilibration
- Use energy plots for actual values

### 4. **Hydrogen Bonds**
- Look for dashed lines between N-H and C=O
- Should form and break dynamically

**VMD:**
```tcl
mol representation HBonds 3.0 30
mol selection "protein"
mol addrep 0
```

**PyMOL:**
```python
distance hbonds, polymer, polymer, mode=2
```

---

## Real-Time Energy/Temperature Monitoring

### Use GROMACS gmx energy

**Extract temperature:**
```bash
# Create input file
echo "Temperature" > temp_input.txt
echo "0" >> temp_input.txt

# Extract data
docker run --rm -v "$(pwd):/workspace" -w /workspace/01_md_preparation \
  gromacs/gromacs:latest \
  sh -c "gmx energy -f md.edr -o temperature.xvg < /workspace/temp_input.txt"

# Plot with any tool (Python, Excel, etc.)
```

**Extract pressure:**
```bash
echo "Pressure" > press_input.txt
echo "0" >> press_input.txt

docker run --rm -v "$(pwd):/workspace" -w /workspace/01_md_preparation \
  gromacs/gromacs:latest \
  sh -c "gmx energy -f md.edr -o pressure.xvg < /workspace/press_input.txt"
```

**Extract energy:**
```bash
echo "Total-Energy" > energy_input.txt
echo "0" >> energy_input.txt

docker run --rm -v "$(pwd):/workspace" -w /workspace/01_md_preparation \
  gromacs/gromacs:latest \
  sh -c "gmx energy -f md.edr -o energy.xvg < /workspace/energy_input.txt"
```

---

## Quick PyMOL Script for Instant Visualization

Save as `view_md.py`:

```python
from pymol import cmd

# Load files
cmd.load('01_md_preparation/npt.gro', 'system')
cmd.load_traj('01_md_preparation/md.xtc', 'system')

# Protein representation
cmd.hide('everything', 'system')
cmd.show('cartoon', 'system and polymer')
cmd.color('cyan', 'system and polymer')

# Water representation (sampling for performance)
cmd.show('nb_spheres', 'system and resn SOL and id 1-1000')
cmd.set('sphere_scale', 0.2, 'resn SOL')
cmd.color('red', 'system and resn SOL and elem O')

# Ions
cmd.show('spheres', 'system and resn NA')
cmd.color('purple', 'system and resn NA')

# View setup
cmd.zoom('system and polymer', buffer=10)
cmd.bg_color('white')
cmd.set('ray_shadows', 0)
cmd.set('depth_cue', 0)

# Play movie
cmd.mplay()

print("MD trajectory loaded!")
print("Controls:")
print("  Space: Play/Pause")
print("  <- ->: Step frames")
print("  Mouse: Rotate view")
```

Run with: `pymol view_md.py`

---

## Troubleshooting

### "No trajectory frames loaded"
- MD hasn't produced XTC file yet (needs ~10 ps)
- Check file exists: `ls -lh 01_md_preparation/md.xtc`

### "Out of memory"
- Too many water molecules visible
- Show only subset: `select waters, resn SOL and id 1-1000`

### "Trajectory not updating"
- XTC file is buffered, updates every 10 ps
- GROMACS writes data in chunks

### "Everything looks frozen"
- Timesteps are 2 fs - very fast!
- Each frame is 10 ps (5000 steps)
- Peptide motion is subtle at this timescale

---

## Recommended Workflow

**For Quick Check (Now):**
1. Wait until XTC file appears (~10 ps = ~2 minutes)
2. Open in PyMOL with view_md.py script
3. Watch a few frames to confirm it's working

**For Nice Movie (Later):**
1. Let simulation run to ~1-10 ns
2. Use ChimeraX + TheVisualHub scripts
3. Apply UltimateSmoothMD for smooth motion
4. Render 4K movie for presentations

**For Analysis (After completion):**
1. Use MDAnalysis for RMSD, clustering
2. Use gmx energy for thermodynamics
3. Create matplotlib plots for publication

---

## Current Status Check

Check if trajectory file exists:
```bash
cd C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation

# Check for XTC file
ls -lh md.xtc

# If exists, check how many frames:
docker run --rm -v "$(pwd):/workspace" -w /workspace \
  gromacs/gromacs:latest \
  gmx check -f md.xtc
```

---

**TL;DR:**
- **PyMOL (easiest):** `pymol npt.gro md.xtc` then hit Space to play
- **VMD (best features):** Load GRO, load XTC, Graphics → Representations
- **ChimeraX (best movies):** Use TheVisualHub scripts for 4K output
- **MDAnalysis (quantitative):** Python scripts for RMSD, energy analysis

**Right now:** Wait ~10 minutes for XTC file to appear, then load in PyMOL!
