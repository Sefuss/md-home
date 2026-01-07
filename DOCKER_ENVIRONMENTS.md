# Docker Environment Specification

**CRITICAL**: Use these exact versions for reproducibility across machines.

---

## GROMACS (Molecular Dynamics)

### Docker Image
```bash
# Image details
Repository: gromacs/gromacs
Tag:        latest
Digest:     sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3
Pulled:     2026-01-05
Created:    2022-07-04
```

### GROMACS Version
```
GROMACS version:    2022.2
Precision:          mixed
Memory model:       64 bit
MPI library:        MPI (CUDA-aware)
OpenMP support:     enabled (GMX_OPENMP_MAX_THREADS = 128)
GPU support:        CUDA
SIMD instructions:  AVX2_256
CPU FFT library:    fftw-3.3.8-sse2-avx-avx2-avx2_128
GPU FFT library:    cuFFT
```

### Pull Exact Image
```bash
# Option 1: By digest (recommended for reproducibility)
docker pull gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3

# Option 2: Tag image locally for convenience
docker pull gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3
docker tag gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3 gromacs/gromacs:2022.2-avitag

# Option 3: Verify existing image matches
docker images gromacs/gromacs --digests
# Should show: sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3
```

---

## Active Simulations

### Free Peptide MD (Completed)
- **Directory**: `01_md_preparation/`
- **Duration**: 23 ns
- **Status**: ✅ Complete
- **Container**: N/A (finished)

### Fusion Protein MD (Running)
- **Directory**: `03_fusion_md/`
- **Duration**: 100 ns (target)
- **Status**: 🔄 Running
- **Container**: `fusion_md_production_v2`
- **Started**: 2026-01-05
- **ETA**: ~3 weeks (2026-02-10)
- **Command**:
```bash
docker run --name fusion_md_production_v2 -d \
  -v "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3 \
  gmx mdrun -v -deffnm md
```

---

## Resuming at Different Location

### Transfer Files
```bash
# Copy entire approach_C directory to new machine
# Include:
03_fusion_md/md.xtc      # Trajectory (partial if still running)
03_fusion_md/md.cpt      # Checkpoint (resume point)
03_fusion_md/md.log      # Progress log
03_fusion_md/npt.gro     # Structure
03_fusion_md/topol.top   # Topology
03_fusion_md/md.mdp      # MD parameters
03_fusion_md/md.tpr      # Run input file
```

### Resume Simulation
```bash
# On new machine:

# 1. Pull exact same Docker image
docker pull gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3

# 2. Navigate to directory
cd /path/to/03_fusion_md

# 3. Resume from checkpoint
docker run --rm \
  -v "$(pwd):/workspace" \
  -w /workspace \
  gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3 \
  gmx mdrun -v -deffnm md -cpi md.cpt
```

**NOTE**: GROMACS automatically resumes from checkpoint if `md.cpt` exists.

---

## Verification

### Check GROMACS Version
```bash
docker run --rm gromacs/gromacs@sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3 \
  gmx --version
```

Should output:
```
GROMACS version:    2022.2
Precision:          mixed
```

### Check Image Digest
```bash
docker images gromacs/gromacs --digests
```

Should show: `sha256:ed530212050a891afac6ad785d8d1375a277cc1736fda349ee5bbedf87a369f3`

---

## Why This Matters

**Reproducibility**: Scientific simulations must be reproducible. Using exact Docker digests ensures:
- Same GROMACS version (2022.2)
- Same compilation flags
- Same numerical libraries
- Identical results across machines

**Publication**: Document this in methods section:
```
"Molecular dynamics simulations were performed using GROMACS 2022.2
(Docker: gromacs/gromacs@sha256:ed53021...) with mixed precision,
CUDA GPU acceleration, and MPI parallelization."
```

---

## Other Tools (For Reference)

### ColabFold (AlphaFold2)
```bash
# Used for initial structure prediction
ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2
```

### ChimeraX
```bash
# Not dockerized - native installation
# Version: UCSF ChimeraX 1.11 (2025-12-17)
# Location: C:\Program Files\ChimeraX 1.11\bin\chimerax.exe
```

---

**Last Updated**: 2026-01-05
**Maintained By**: bmartinez (BPS Bioscience)
