# Reproducing and Continuing MD Simulation at Home

**Goal:** Take the MD work from your work PC and seamlessly continue the simulation from checkpoint on your home PC

---

## TL;DR - Quick Start at Home

```bash
# 1. Clone repository
git clone https://github.com/Sefuss/md-home.git
cd md-home

# 2. Pull exact Docker image
docker pull gromacs/gromacs:2022.2

# 3. Continue from checkpoint (Windows Git Bash)
cd 03_fusion_md
export MSYS_NO_PATHCONV=1
docker run --rm -v "$(pwd)/..:/workspace" -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 gmx mdrun -v -deffnm md -cpi md.cpt

# That's it! Simulation continues from step 1,330,000
```

---

## Why This Works: Understanding Checkpoints

### What is a Checkpoint File (.cpt)?

GROMACS checkpoint files contain **the complete state** of your simulation at a specific moment:
- Atomic positions and velocities
- Current simulation step number
- Thermostat and barostat states
- Random number generator state
- Energy conservation tracking

**The magic:** When you run `gmx mdrun -cpi md.cpt`, GROMACS:
1. Reads the checkpoint file
2. Restores the exact state from step 1,330,000
3. Continues as if nothing happened
4. Appends new frames to `md.xtc`
5. Creates new checkpoints every 15 minutes

---

## Exact Software Versions (CRITICAL for Reproducibility)

### Docker Image

**What work PC is using:**
```bash
Image: gromacs/gromacs:latest
Tag pulled on: January 2, 2026
Actual version: 2022.2
```

**For home PC, use exact version:**
```bash
docker pull gromacs/gromacs:2022.2
```

**Verify version matches:**
```bash
docker run --rm gromacs/gromacs:2022.2 gmx --version

# Should show:
# GROMACS version:    2022.2
# Precision:          mixed
# MPI library:        MPI (CUDA-aware)
# GPU support:        CUDA
```

### Why Version Matters

- Different GROMACS versions use different checkpoint formats
- Mixing versions = "Checkpoint file incompatible" error
- Using 2022.2 guarantees checkpoint compatibility

---

## What's in the Repository

### Files You NEED (Included in Git)

These are the essential files to continue the simulation:

| File | Size | Purpose |
|------|------|---------|
| `03_fusion_md/md.cpt` | 2.5 MB | ⭐ Current checkpoint (step 1,330,000) |
| `03_fusion_md/md_prev.cpt` | 2.5 MB | Backup checkpoint |
| `03_fusion_md/md.tpr` | 3.1 MB | Simulation input file (portable) |
| `03_fusion_md/topol.top` | 617 KB | System topology |
| `03_fusion_md/npt.gro` | ~200 KB | Last equilibrated structure |
| `03_fusion_md/*.mdp` | <10 KB | Parameter files |

### Files Excluded from Git (Too Large)

These are NOT in git because they're huge:

| File | Size | Excluded? |
|------|------|-----------|
| `03_fusion_md/md.xtc` | ~500 MB | ❌ Not in git (too large) |
| `03_fusion_md/md.log` | ~50 MB | ❌ Not in git |
| `03_fusion_md/md.edr` | ~100 MB | ❌ Not in git |

**Don't worry!** When you continue the simulation:
- `md.xtc` will be created fresh (or appended if it exists)
- `md.log` will be created fresh
- Energy will be conserved correctly because checkpoint has all state info

---

## Step-by-Step: Setting Up Home PC

### 1. Install Prerequisites

**Windows:**
1. Install Docker Desktop: https://www.docker.com/products/docker-desktop
   - Enable WSL2 backend
   - Allocate at least 8 GB RAM to Docker
2. Install Git: https://git-scm.com/downloads
3. Install Git Bash (comes with Git for Windows)

**Mac:**
1. Install Docker Desktop for Mac
2. Git is pre-installed on macOS

**Linux:**
1. Install Docker: `sudo apt-get install docker.io`
2. Install Git: `sudo apt-get install git`

### 2. Clone Repository

```bash
git clone https://github.com/Sefuss/md-home.git
cd md-home
```

### 3. Pull Exact Docker Image

```bash
docker pull gromacs/gromacs:2022.2
```

**Verify:**
```bash
docker run --rm gromacs/gromacs:2022.2 gmx --version | head -20
```

Should show:
```
                       :-) GROMACS - gmx_mpi, 2022.2 (-:

GROMACS version:    2022.2
Precision:          mixed
Memory model:       64 bit
MPI library:        MPI (CUDA-aware)
GPU support:        CUDA
```

### 4. Verify Checkpoint Files

```bash
cd 03_fusion_md
ls -lh *.cpt *.tpr topol.top
```

You should see:
```
-rw-r--r-- md.cpt         (2.5 MB)
-rw-r--r-- md_prev.cpt    (2.5 MB)
-rw-r--r-- md.tpr         (3.1 MB)
-rw-r--r-- topol.top      (617 KB)
```

---

## Running the Simulation at Home

### Option 1: Foreground (For Testing)

**Windows (Git Bash or PowerShell):**
```bash
cd 03_fusion_md
export MSYS_NO_PATHCONV=1

docker run --rm \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  gmx mdrun -v -deffnm md -cpi md.cpt
```

**Mac/Linux:**
```bash
cd 03_fusion_md

docker run --rm \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  gmx mdrun -v -deffnm md -cpi md.cpt
```

### Option 2: Background (For Long Runs)

**Windows:**
```powershell
cd 03_fusion_md

docker run -d --name fusion_md_home `
  -v "${PWD}\..:/workspace" `
  -w /workspace/03_fusion_md `
  gromacs/gromacs:2022.2 `
  bash -c "gmx mdrun -v -deffnm md -cpi md.cpt"
```

**Mac/Linux:**
```bash
cd 03_fusion_md

docker run -d --name fusion_md_home \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  bash -c "gmx mdrun -v -deffnm md -cpi md.cpt"
```

### Option 3: With GPU (Much Faster!)

If your home PC has an NVIDIA GPU:

**Windows:**
```bash
docker run -d --name fusion_md_home --gpus all \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  bash -c "gmx mdrun -v -deffnm md -cpi md.cpt"
```

**Linux (requires nvidia-docker2):**
```bash
docker run -d --name fusion_md_home --gpus all \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  bash -c "gmx mdrun -v -deffnm md -cpi md.cpt"
```

---

## Monitoring Progress

### Check if Container is Running

```bash
docker ps
```

Should show container `fusion_md_home` with status "Up"

### View Latest Output

```bash
# Last 30 lines
docker logs --tail 30 fusion_md_home

# Last 50 lines
docker logs --tail 50 fusion_md_home

# Follow live (Ctrl+C to stop)
docker logs -f fusion_md_home
```

### Check Current Simulation Step

```bash
docker exec fusion_md_home tail -30 md.log | grep "Step"
```

Example output:
```
           Step           Time
        1500000     3000.00000
```

Meaning: Currently at step 1,500,000 (3 ns)

### Calculate Progress

Total steps: 50,000,000
Current step: (check from log)
Progress %: (current / 50,000,000) × 100

Example:
- Step 1,330,000 → 2.66% (where work PC left off)
- Step 25,000,000 → 50%
- Step 50,000,000 → 100% (done!)

---

## Current Simulation Status

### Where Work PC Left Off

**Progress:** Step 1,330,000 / 50,000,000 (2.66% complete)
**Time simulated:** 2.66 ns / 100 ns target
**Checkpoint saved:** January 7, 2026 at 13:50

**What this means:**
- Home PC will continue from step 1,330,000
- Remaining: 48,670,000 steps (97.34%)
- Remaining time: 97.34 ns

### Performance Expectations

**Work PC (CPU-only, slow):**
- Speed: ~200-500 steps/second
- Time to completion: ~60 days

**Home PC with good GPU:**
- Speed: Could be 2000-5000+ steps/second
- Time to completion: ~6-24 hours (10-50× faster!)

**Check your speed:**
```bash
docker logs --tail 50 fusion_md_home | grep "ns/day"
```

Example output:
```
Performance: 15.2 ns/day
```

---

## Troubleshooting

### Issue: Checkpoint incompatible error

```
Fatal error:
Checkpoint file incompatible
```

**Cause:** GROMACS version mismatch

**Solution 1:** Use exact version 2022.2
```bash
docker pull gromacs/gromacs:2022.2
# Update docker run commands to use :2022.2
```

**Solution 2:** If that fails, restart from npt.gro
```bash
docker run --rm -v "$(pwd)/..:/workspace" -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  gmx grompp -f md.mdp -c npt.gro -p topol.top -o md_new.tpr

docker run --rm -v "$(pwd)/..:/workspace" -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  gmx mdrun -v -deffnm md_new
```

Note: This restarts from equilibration (step 0), losing progress.

### Issue: File path not found (Windows)

```
Error: Cannot find /workspace/03_fusion_md
```

**Cause:** Windows path conversion issue

**Solution:** Set MSYS_NO_PATHCONV=1
```bash
export MSYS_NO_PATHCONV=1
# Then run docker command again
```

### Issue: Docker can't access GPU

**Windows:**
1. Open Docker Desktop
2. Settings → Resources → WSL Integration
3. Enable integration with your WSL distro

**Linux:**
1. Install nvidia-docker2:
```bash
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update
sudo apt-get install -y nvidia-docker2
sudo systemctl restart docker
```

2. Test GPU access:
```bash
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
```

### Issue: Out of disk space

**Check trajectory size:**
```bash
ls -lh 03_fusion_md/md.xtc
```

**100 ns trajectory = ~1-2 GB**

**Solution:** Free up at least 5 GB before starting

### Issue: Simulation too slow

**Check performance:**
```bash
docker logs --tail 100 fusion_md_home | grep "Performance"
```

**If < 5 ns/day:**
- Make sure GPU is being used (check docker logs for "GPU" messages)
- Close other applications
- Allocate more CPU/RAM to Docker (Docker Desktop settings)

**If > 20 ns/day:**
- Great! 100 ns should finish in ~5 days

---

## Controlling the Simulation

### Stop Simulation

```bash
docker stop fusion_md_home
```

**Safe!** Checkpoint is saved automatically every 15 minutes.
- Last checkpoint: `md.cpt`
- Previous checkpoint: `md_prev.cpt`

### Restart from Stopped Container

```bash
docker start fusion_md_home
```

**IMPORTANT:** This continues the same container, not necessarily from the latest checkpoint!

**Better:** Start fresh container with latest checkpoint:
```bash
docker rm fusion_md_home  # Remove old container

# Start new one from latest checkpoint
docker run -d --name fusion_md_home \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 \
  bash -c "gmx mdrun -v -deffnm md -cpi md.cpt"
```

### Pause and Save Checkpoint Manually

```bash
# This sends SIGTERM, GROMACS saves checkpoint before exit
docker stop fusion_md_home

# Verify checkpoint was updated
ls -lh 03_fusion_md/md.cpt
```

---

## After Simulation Completes

### How to Know it's Done

```bash
docker logs --tail 50 fusion_md_home
```

Look for:
```
Finished mdrun
```

Or check final step:
```
           Step           Time
       50000000   100000.00000
```

### Extract Protein-Only Trajectory for Visualization

```bash
cd 01_md_preparation

# Center and fix PBC
printf "1\n0\n" | docker run --rm -i \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/01_md_preparation \
  gromacs/gromacs:2022.2 \
  gmx trjconv -s ../03_fusion_md/md.tpr -f ../03_fusion_md/md.xtc \
  -o md_centered.xtc -center -pbc mol -ur compact

# Extract protein only
printf "1\n" | docker run --rm -i \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/01_md_preparation \
  gromacs/gromacs:2022.2 \
  gmx trjconv -s ../03_fusion_md/md.tpr -f md_centered.xtc \
  -o protein_only.xtc

# Extract protein-only structure
printf "1\n" | docker run --rm -i \
  -v "$(pwd)/..:/workspace" \
  -w /workspace/01_md_preparation \
  gromacs/gromacs:2022.2 \
  gmx trjconv -s ../03_fusion_md/md.tpr -f ../03_fusion_md/md.xtc \
  -o protein_only.gro -dump 0
```

### Create Publication Movie

See `../../TRAJECTORY_MOVIE_GUIDE.md` for complete instructions.

**Quick version:**
```bash
cd visualization

# Windows
"C:\Program Files\ChimeraX 1.11\bin\chimerax.exe" --script lag16_fusion_movie_publication.py

# Mac
/Applications/ChimeraX.app/Contents/bin/chimerax --script lag16_fusion_movie_publication.py
```

Output: `04_movies/lag16_fusion_publication_HQ.mp4`

---

## Transferring Work Back to Work PC

### Option 1: Push to GitHub

```bash
cd approach_C_peptide_ensemble

# Add checkpoint files
git add 03_fusion_md/md.cpt
git add 03_fusion_md/md_prev.cpt

# Commit
git commit -m "Checkpoint at step X,XXX,XXX"

# Push
git push origin main
```

Then on work PC:
```bash
git pull origin main
```

### Option 2: External Drive

Copy these files to external drive:
- `03_fusion_md/md.cpt` (2.5 MB)
- `03_fusion_md/md.xtc` (current trajectory, ~1 GB)
- `03_fusion_md/md.log` (log file, ~50 MB)

Then on work PC, continue from the transferred checkpoint!

---

## Simulation Parameters Reference

### System Information

**Protein:** LaG16-AviTag fusion
- Residues: 149
- Atoms (protein): 2,254
- Atoms (total with water): 106,938

**Box:**
- Type: Cubic
- Size: ~6-7 nm edges
- Water: TIP3P
- Ions: Added to neutralize

### Production MD Settings

From `md.mdp`:

```
integrator  = md                ; Leap-frog MD
nsteps      = 50000000          ; 100 ns total (2 fs timestep)
dt          = 0.002             ; 2 fs per step

; Temperature coupling (V-rescale)
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1    0.1
ref_t       = 300    300        ; 300 K

; Pressure coupling (Parrinello-Rahman)
pcoupl      = Parrinello-Rahman
tau_p       = 2.0
ref_p       = 1.0               ; 1 bar
compressibility = 4.5e-5

; Output control
nstxout-compressed = 5000       ; Save frame every 10 ps
```

### Timeline

**Total simulation:**
- Steps: 50,000,000
- Time: 100 ns
- Frames: 10,000 (one every 10 ps)

**Work PC progress:**
- Steps completed: 1,330,000 (2.66%)
- Time simulated: 2.66 ns
- Frames saved: 266

**Home PC will do:**
- Remaining steps: 48,670,000 (97.34%)
- Remaining time: 97.34 ns
- New frames: ~9,734

---

## Expected Performance

### CPU-Only

**Typical speeds:**
- Desktop CPU: 5-15 ns/day
- Laptop CPU: 2-8 ns/day

**100 ns completion:**
- Fast CPU: ~7 days
- Slow CPU: ~50 days

### With GPU

**Typical speeds:**
- GTX 1660: 15-25 ns/day
- RTX 3060: 25-40 ns/day
- RTX 3080: 40-60 ns/day
- RTX 4090: 80-120 ns/day

**100 ns completion:**
- Good GPU: 2-7 days
- Great GPU: 1-2 days

### Check Your Speed

After 10-30 minutes of running:
```bash
docker logs --tail 100 fusion_md_home | grep "Performance"
```

Example:
```
Performance: 35.2 ns/day
```

→ 100 ns will take: 100 / 35.2 = 2.8 days

---

## Best Practices

### 1. Test First

Run for 10-15 minutes in foreground mode first:
```bash
docker run --rm -v "$(pwd)/..:/workspace" -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 gmx mdrun -v -deffnm md -cpi md.cpt -nsteps 50000
```

This runs 50,000 steps (~100 ps) as a test.

### 2. Check Checkpoint Saves

Checkpoints are saved every 15 minutes by default. Verify:
```bash
watch -n 60 "ls -lh 03_fusion_md/md.cpt"
```

Timestamp should update every 15 min.

### 3. Monitor Temperature

Every few hours:
```bash
docker logs --tail 100 fusion_md_home | grep "Temperature"
```

Should be ~300 K. If drifting significantly, something is wrong.

### 4. Save Checkpoints Externally

Every 24 hours, copy checkpoint to external drive:
```bash
cp 03_fusion_md/md.cpt ~/Desktop/backup_checkpoints/md_$(date +%Y%m%d).cpt
```

### 5. Check Disk Space

Before starting and periodically:
```bash
df -h .
```

Make sure you have at least 5 GB free at all times.

---

## FAQ

**Q: Can I use a different GROMACS version?**
A: Not recommended. Different versions have incompatible checkpoints. Stick to 2022.2.

**Q: What if I stop mid-simulation?**
A: No problem! GROMACS saves checkpoints every 15 minutes. Restart with `-cpi md.cpt`.

**Q: Can I transfer checkpoints between computers?**
A: YES! That's the whole point. As long as both use GROMACS 2022.2, checkpoints are portable.

**Q: Will starting fresh lose my work PC's progress?**
A: No! The checkpoint file `md.cpt` contains all the work PC's progress (step 1,330,000).

**Q: What if home PC is faster? Will it conflict?**
A: No conflict. The simulation doesn't know or care about speed. It just continues from the checkpoint.

**Q: Can I run on multiple PCs simultaneously?**
A: NO! Don't do this. Each run should start from the latest checkpoint from the previous run.

**Q: What's the difference between md.cpt and md_prev.cpt?**
A:
- `md.cpt` = latest checkpoint
- `md_prev.cpt` = previous checkpoint (backup)

Always use `md.cpt`.

---

## Next Steps After Completion

1. **Analyze trajectory**
   - RMSD, RMSF, radius of gyration
   - Secondary structure analysis
   - Clustering

2. **Create movies**
   - Use `visualization/lag16_fusion_movie_publication.py`
   - Try different colors from `PUBLICATION_COLOR_PALETTES.md`

3. **Extract representative structures**
   - Cluster trajectory
   - Get 5-10 representative conformations

4. **Compare to free peptide**
   - Load both trajectories in ChimeraX
   - Use `visualization/chimerax_side_by_side_comparison.py`

5. **RFDiffusion3 binder design**
   - Use extracted conformations as input
   - Design binders specific to each cluster

---

## Resources

- **GROMACS Manual:** https://manual.gromacs.org/
- **GROMACS Tutorials:** http://www.mdtutorials.com/gmx/
- **Docker Documentation:** https://docs.docker.com/
- **This Project's Main README:** See `README.md`
- **Visualization Guide:** See `TRAJECTORY_MOVIE_GUIDE.md`
- **Color Palettes:** See `PUBLICATION_COLOR_PALETTES.md`

---

## Summary

**What you're doing:**
Taking a partially-complete MD simulation (2.66% done on work PC) and continuing it seamlessly on your home PC using GROMACS checkpoints.

**What you need:**
- Docker with GROMACS 2022.2
- Checkpoint files from git repository
- ~5 GB free disk space
- Patience (could take 1-50 days depending on hardware)

**Key command:**
```bash
docker run --rm -v "$(pwd)/..:/workspace" -w /workspace/03_fusion_md \
  gromacs/gromacs:2022.2 gmx mdrun -v -deffnm md -cpi md.cpt
```

**Result:**
A complete 100 ns MD trajectory of LaG16-AviTag fusion, ready for analysis and binder design!

---

**Created:** 2026-01-07
**Work PC Progress:** Step 1,330,000 / 50,000,000 (2.66%)
**Last Updated:** 2026-01-07
