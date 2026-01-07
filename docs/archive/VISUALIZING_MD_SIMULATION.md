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
