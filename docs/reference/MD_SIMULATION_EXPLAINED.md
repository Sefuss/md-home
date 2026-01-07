# MD Simulation Explained: Complete Workflow

**Purpose:** Run molecular dynamics simulation on LaG16-AviTag fusion protein to see how the AviTag behaves when attached to a nanobody

**Date:** January 2, 2026

---

## Part 1: The Big Picture

### What Are We Trying To Do?

**Goal:** Simulate 30-50 nanoseconds of the fusion protein moving in water at body temperature (300 K, ~27°C)

**Why?** To answer: Does the AviTag peptide (residues 135-149) behave differently when attached to a protein vs. free in solution?

**Approach:** Use physics equations to calculate how every atom moves, 2 femtoseconds at a time, for billions of timesteps

---

## Part 2: Docker Side (The Container)

### What Is Docker Doing?

**Think of Docker as a portable computer inside your computer.**

```
Your Windows PC
   ├─ Docker Desktop (container manager)
   │    └─ GROMACS Container (portable Linux environment)
   │         ├─ GROMACS 2022.2 installed
   │         ├─ Force field files (physics parameters)
   │         └─ All dependencies pre-configured
   └─ Your files (Windows file system)
```

### Why Use Docker?

1. **No installation headaches** - GROMACS is complex to install on Windows
2. **Reproducibility** - Everyone uses the exact same GROMACS version
3. **Portability** - Same container works on Windows, Mac, Linux
4. **Isolation** - Doesn't mess with your system

### Anatomy of a Docker Command

Let's break down a typical command:

```bash
docker run --rm \
  -v "C:\Users\bmartinez\...\03_fusion_md:/workspace" \
  -w /workspace/03_fusion_md \
  gromacs/gromacs:latest \
  gmx pdb2gmx -f fusion.pdb -o output.gro
```

**What each part means:**

```
docker run              Start a new container
--rm                    Delete container when done (cleanup)
-v "C:\...\:/workspace" Mount your Windows folder inside container
                        (Your files ↔ Container can see them)
-w /workspace/...       Set working directory inside container
gromacs/gromacs:latest  Which container image to use
gmx pdb2gmx ...         The actual GROMACS command to run
```

**Volume mounting explained:**

```
Windows path:     C:\Users\bmartinez\...\03_fusion_md\fusion.pdb
                          ↓ (mounted as)
Container path:   /workspace/03_fusion_md/fusion.pdb
                          ↑
GROMACS sees this path inside the container
```

When GROMACS writes output, it goes to `/workspace/...` inside the container, which is actually your `C:\Users\...` folder on Windows!

---

## Part 3: GROMACS Side (The Simulation Engine)

### The Complete Workflow

```
START: PDB file (protein structure, just atomic coordinates)
   ↓
[Step 1] pdb2gmx: Add physics (topology)
   ↓
[Step 2] editconf: Create simulation box
   ↓
[Step 3] solvate: Add water molecules
   ↓
[Step 4] genion: Add ions (neutralize charge)
   ↓
[Step 5] grompp: Prepare energy minimization
   ↓
[Step 6] mdrun: Energy minimization
   ↓
[Step 7] grompp: Prepare NVT equilibration
   ↓
[Step 8] mdrun: NVT equilibration (constant volume, temp)
   ↓
[Step 9] grompp: Prepare NPT equilibration
   ↓
[Step 10] mdrun: NPT equilibration (constant pressure, temp)
   ↓
[Step 11] grompp: Prepare production MD
   ↓
[Step 12] mdrun: Production MD (the actual simulation!)
   ↓
END: Trajectory file (movie of all atomic positions over time)
```

Let's go through each step in detail:

---

### STEP 1: pdb2gmx (Add Physics)

**Input:** `lag16_avitag_fusion.pdb` (149 residues, just coordinates)

**What it does:**
1. Reads amino acid sequence
2. Assigns force field parameters (AMBER99SB-ILDN)
   - Each atom type gets: mass, charge, bond lengths, angles
   - This is the "physics" - how atoms attract/repel each other
3. Adds hydrogen atoms (PDB files often don't include them)
4. Identifies special features (disulfide bonds, histidine protonation)

**Output files:**
- `fusion_processed.gro` - Coordinates WITH hydrogens (2,254 atoms)
- `topol.top` - Topology file (the physics parameters)
- `posre.itp` - Position restraints (for equilibration)

**Key decision:** Force field choice
- We chose AMBER99SB-ILDN (well-validated for proteins)
- Alternative: CHARMM, GROMOS, OPLS-AA
- Different force fields = slightly different physics

**Analogy:** Like going from a sketch of a building to a blueprint with engineering specs (materials, weights, structural connections)

---

### STEP 2: editconf (Create Simulation Box)

**Input:** `fusion_processed.gro` (protein in space)

**What it does:**
1. Centers the protein
2. Creates a cubic box around it
3. Distance from protein edge to box edge: 1.2 nm
   - Why? Protein shouldn't "see itself" through periodic boundaries

**Output:** `box.gro` (protein in a defined cubic box)

**Box dimensions:** ~10.3 nm × 10.3 nm × 10.3 nm

**Periodic boundaries explained:**
```
Imagine the box repeats infinitely in all directions

  [Box] [Box] [Box]
  [Box] [YOU] [Box]   ← Protein "sees" copies of itself
  [Box] [Box] [Box]

If protein gets too close to box edge, it will interact
with its own periodic image (BAD!)
```

**Analogy:** Like setting the stage size for a play - need enough room for actors to move without bumping into the wings.

---

### STEP 3: solvate (Add Water)

**Input:** `box.gro` (empty box with protein)

**What it does:**
1. Fills box with pre-equilibrated water molecules (TIP3P model)
2. Removes water molecules that overlap with protein
3. Adds ~34,896 water molecules (104,688 atoms)

**Output:** `solvated.gro` (protein + water, 106,942 atoms total)

**Water model (TIP3P):**
- 3-point model: Oxygen + 2 Hydrogens
- Each water has partial charges (+/- for electrostatics)
- Optimized to reproduce real water properties (density, diffusion)

**Why so much water?**
- Proteins don't float in vacuum - need realistic environment
- Water participates in dynamics (hydrogen bonds, solvation)
- Provides proper dielectric environment

**Analogy:** Like filling a fish tank with water after placing decorations. Fish (protein) needs water to swim!

---

### STEP 4: genion (Add Ions)

**Input:** `solvated.gro` (protein + water, net charge -2e)

**What it does:**
1. Calculates total system charge (-2.0 from amino acids)
2. Replaces 2 water molecules with 2 Na+ ions
3. System is now electrically neutral (total charge = 0)

**Output:** `system.gro` (neutral system ready for simulation)

**Why neutralize?**
- Ewald electrostatics (used in simulation) assumes neutral system
- Non-neutral systems get artifacts (ions drift to low-dielectric regions)
- Mimics physiological conditions (salt in solution)

**Could add more salt (NaCl):**
- Physiological concentration: ~150 mM NaCl
- For simplicity, we just neutralize (0 mM)

**Analogy:** Like adding salt to water to match ocean salinity - proteins evolved in salty environments!

---

### STEP 5-6: Energy Minimization

#### STEP 5: grompp (Prepare Minimization)

**Inputs:**
- `minim.mdp` - Parameters for minimization
- `system.gro` - Starting coordinates
- `topol.top` - Force field parameters

**What it does:**
1. Reads all input files
2. Checks for errors
3. Generates binary run file (`.tpr`)

**Output:** `em.tpr` (portable run input)

**Key parameters in minim.mdp:**
```
integrator  = steep         (Steepest descent algorithm)
emtol       = 1000.0        (Stop when forces < 1000 kJ/mol/nm)
nsteps      = 50000         (Max 50,000 steps)
```

#### STEP 6: mdrun (Run Minimization)

**Input:** `em.tpr`

**What it does:**
1. Calculates forces on every atom
2. Moves atoms downhill (toward lower energy)
3. Repeats until forces are small

**Output files:**
- `em.gro` - Energy-minimized coordinates
- `em.edr` - Energy data
- `em.log` - Log file

**Why minimize?**
- Initial structure has bad contacts (atoms too close)
- Would cause simulation to explode if not fixed
- Minimize = relax structure to nearest energy minimum

**Analogy:** Like ironing out wrinkles in fabric before sewing. Need smooth starting point!

**Result for our system:**
- Took 773 steps to converge
- Starting energy: -959,041 kJ/mol
- Final energy: -1,670,000 kJ/mol (much lower = more stable!)

---

### STEP 7-8: NVT Equilibration (Constant Volume, Temperature)

#### STEP 7: grompp (Prepare NVT)

**Inputs:**
- `nvt.mdp` - NVT parameters
- `em.gro` - Minimized coordinates
- `topol.top` - Topology

**Key parameters in nvt.mdp:**
```
integrator  = md            (Leap-frog molecular dynamics)
dt          = 0.002         (2 femtosecond timesteps)
nsteps      = 50000         (100 ps total: 50000 × 0.002 ps)
tcoupl      = V-rescale     (Temperature coupling algorithm)
ref_t       = 300           (Target temperature: 300 K)
gen_vel     = yes           (Generate initial velocities)
```

**Output:** `nvt.tpr`

#### STEP 8: mdrun (Run NVT)

**Input:** `nvt.tpr`

**What it does:**
1. Generates random velocities for all atoms (300 K distribution)
2. Runs MD for 100 ps
3. Constantly adjusts velocities to maintain 300 K

**Why NVT first?**
- Protein starts at 0 K (no motion) after minimization
- Need to "heat up" system to 300 K
- Keep volume constant so water doesn't expand too fast
- Position restraints on protein backbone (let solvent relax first)

**Position restraints:**
- Protein heavy atoms held in place with weak springs
- Water is free to move and equilibrate around protein
- Prevents protein from unfolding during heating

**Output files:**
- `nvt.gro` - Coordinates after 100 ps
- `nvt.edr` - Energy, temperature data
- `nvt.cpt` - Checkpoint (can restart if crashes)
- `nvt.log` - Log with performance metrics

**Typical performance:** 5-10 minutes on CPU

**Analogy:** Like warming up a car engine - start slow, let everything reach operating temperature before driving!

---

### STEP 9-10: NPT Equilibration (Constant Pressure, Temperature)

#### STEP 9: grompp (Prepare NPT)

**Inputs:**
- `npt.mdp` - NPT parameters
- `nvt.gro` - NVT-equilibrated coordinates
- `topol.top` - Topology

**Key parameters in npt.mdp:**
```
continuation    = yes       (Continue from NVT)
gen_vel         = no        (Keep existing velocities)
tcoupl          = V-rescale (Temperature coupling)
ref_t           = 300       (300 K)
pcoupl          = Berendsen (Pressure coupling algorithm)
ref_p           = 1.0       (1 bar = atmospheric pressure)
```

**Output:** `npt.tpr`

#### STEP 10: mdrun (Run NPT)

**Input:** `npt.tpr`

**What it does:**
1. Runs MD for 100 ps at 300 K
2. Adjusts box size to maintain 1 bar pressure
3. Still has position restraints on protein

**Why NPT after NVT?**
- NVT heated system, but pressure might be wrong
- NPT lets box expand/contract to reach proper density
- Water density should be ~1000 kg/m³ at 300 K, 1 bar

**Pressure coupling:**
- Barostat (pressure-stat, like thermostat for pressure)
- If pressure too high → box expands
- If pressure too low → box contracts
- Box size equilibrates to match 1 bar

**Output files:**
- `npt.gro` - Final equilibrated coordinates
- `npt.edr` - Energy, temperature, pressure, density data
- `npt.cpt` - Checkpoint
- `npt.log` - Log

**What to check:**
- Temperature: Should average 300 K
- Pressure: Should fluctuate around 1 bar (±100 bar is normal)
- Density: Should be ~1000 kg/m³

**Analogy:** Like filling a tire to proper pressure - adjust volume until you hit target pressure!

---

### STEP 11-12: Production MD (The Real Simulation!)

#### STEP 11: grompp (Prepare Production)

**Inputs:**
- `md.mdp` - Production MD parameters
- `npt.gro` - Equilibrated coordinates
- `topol.top` - Topology

**Key parameters in md.mdp:**
```
integrator  = md
dt          = 0.002 ps       (2 fs timesteps)
nsteps      = 25000000       (50,000 ps = 50 ns total)
tcoupl      = V-rescale      (Temperature control)
pcoupl      = Parrinello-Rahman  (More accurate pressure coupling)
ref_t       = 300 K
ref_p       = 1.0 bar

; Output control
nstxout-compressed = 5000    (Save coordinates every 10 ps)
nstenergy          = 5000    (Save energies every 10 ps)
```

**Output:** `md.tpr`

**Why different pressure coupling?**
- Berendsen (used in NPT): Fast equilibration, but not accurate
- Parrinello-Rahman (production): Slower, but gives correct ensemble
- Only use P-R after system is equilibrated

#### STEP 12: mdrun (Production MD)

**Input:** `md.tpr`

**What it does - THE ACTUAL SIMULATION:**

```
For each timestep (0, 2 fs, 4 fs, 6 fs, ...):
    1. Calculate forces on all 106,942 atoms
       - Bonded forces (bonds, angles, dihedrals)
       - Non-bonded (Van der Waals, electrostatics)
       - Constraints (keep bonds rigid)

    2. Update velocities using forces (F = ma)

    3. Update positions using velocities

    4. Apply thermostat (adjust velocities to maintain 300 K)

    5. Apply barostat (adjust box size to maintain 1 bar)

    6. Every 5000 steps (10 ps):
       - Save coordinates to trajectory file
       - Save energies, temperature, pressure

    7. Repeat for 25,000,000 timesteps (50 ns)
```

**Output files:**
- `md.xtc` - **Trajectory file** (compressed coordinates over time)
- `md.edr` - Energy, temperature, pressure time series
- `md.log` - Performance, timing information
- `md.cpt` - Checkpoint (restart if needed)

**Trajectory file explained:**
```
md.xtc contains:
  Frame 0:    Coordinates at 0 ps     (106,942 atoms)
  Frame 1:    Coordinates at 10 ps    (106,942 atoms)
  Frame 2:    Coordinates at 20 ps    (106,942 atoms)
  ...
  Frame 5000: Coordinates at 50 ns    (106,942 atoms)

Total: 5001 frames = "movie" of simulation
```

**Performance estimate:**
- Our CPU (i3-10100): ~27 ns/day
- 50 ns will take: ~2 days
- Larger system (more atoms) = slower
- GPU would be ~6-8x faster (150-200 ns/day)

---

## Part 4: What Happens During Simulation

### The Physics Being Calculated

**Force Field = Physics Equations**

For every atom, at every timestep, GROMACS calculates:

#### 1. Bonded Forces

**Bonds (stretching):**
```
E_bond = k_b × (r - r0)²
  k_b = bond force constant
  r   = current bond length
  r0  = equilibrium bond length
```
Example: C-C bond wants to be ~0.15 nm, spring pulls it back if stretched

**Angles (bending):**
```
E_angle = k_θ × (θ - θ0)²
  k_θ = angle force constant
  θ   = current angle
  θ0  = equilibrium angle
```
Example: H-O-H angle in water wants to be 104.5°

**Dihedrals (rotation):**
```
E_dihedral = k_φ × [1 + cos(n×φ - φ0)]
  Determines rotation around bonds
  Important for protein backbone flexibility
```

#### 2. Non-Bonded Forces

**Van der Waals (Lennard-Jones):**
```
E_LJ = 4ε × [(σ/r)¹² - (σ/r)⁶]
  ε = well depth (how strongly atoms attract)
  σ = distance at zero energy
  r = distance between atoms

  r¹² term: Strong repulsion at close range (atoms don't overlap)
  r⁶ term: Weak attraction at medium range (dispersion forces)
```

**Electrostatics (Coulomb):**
```
E_elec = (q1 × q2) / (4πε0 × r)
  q1, q2 = partial charges on atoms
  r = distance between atoms

Opposite charges attract, like charges repel
```

For 106,942 atoms, this means calculating interactions between:
- All bonded pairs: ~2,277 bonds
- All angles: ~4,074 angles
- All dihedrals: ~6,149 dihedrals
- All atom pairs within cutoff (~500,000 pairs)

**Every 2 femtoseconds!**

That's why it takes days to simulate nanoseconds!

---

### What We Can Learn From The Simulation

After 50 ns of simulation, we'll have data on:

#### 1. AviTag Conformations
- How flexible is the tag when attached?
- Does it prefer specific shapes (compact vs extended)?
- How many distinct conformations does it visit?

#### 2. Comparison to Free Peptide
We already ran 23 ns of free AviTag, so we can compare:

| Metric | Free Peptide | Fusion (To Measure) |
|--------|--------------|---------------------|
| RMSD range | 0-25 Å | ? |
| Average Rg | 12.5 Å | ? |
| Flexibility | High | ? |
| Transitions | 22 (compact ↔ extended) | ? |

**If similar:** Tag is intrinsically flexible regardless of context
**If different:** Fusion partner constrains or influences tag

#### 3. Design Decision
- **Free = Fusion** → Use free peptide ensemble for binder design (easier)
- **Free ≠ Fusion** → Use fusion ensemble (more realistic)

---

## Part 5: File Summary

### Files Created (In Order)

```
Input:
  lag16_avitag_fusion.pdb         Initial structure (149 residues)

Step 1 (pdb2gmx):
  fusion_processed.gro            Coordinates with hydrogens
  topol.top                       Force field parameters
  posre.itp                       Position restraints

Step 2 (editconf):
  box.gro                         Protein in simulation box

Step 3 (solvate):
  solvated.gro                    Protein + water

Step 4 (genion):
  ions.tpr                        Temporary file for ion addition
  system.gro                      Neutralized system

Step 5-6 (Energy Minimization):
  em.tpr                          Minimization input
  em.gro                          Minimized coordinates
  em.edr                          Energy data
  em.log                          Log file

Step 7-8 (NVT Equilibration):
  nvt.tpr                         NVT input
  nvt.gro                         After 100 ps NVT
  nvt.edr                         Temperature, energy data
  nvt.cpt                         Checkpoint
  nvt.log                         Log

Step 9-10 (NPT Equilibration):
  npt.tpr                         NPT input
  npt.gro                         After 100 ps NPT
  npt.edr                         Pressure, density data
  npt.cpt                         Checkpoint
  npt.log                         Log

Step 11-12 (Production MD):
  md.tpr                          Production input
  md.xtc                          TRAJECTORY (the gold!)
  md.edr                          Time series data
  md.cpt                          Checkpoint
  md.log                          Performance log

Analysis (after simulation):
  We'll extract residues 135-149 from md.xtc
  Calculate RMSD, Rg, visualize conformations
```

**Most important file:** `md.xtc` (the trajectory)
**Second most important:** `npt.gro` (starting structure for analysis)

---

## Part 6: Common Questions

### Q1: Why so many equilibration steps?

**A:** Think of it like launching a rocket:
1. **Energy minimization:** Remove bad contacts (pre-flight check)
2. **NVT:** Heat system to 300 K (warm up engines)
3. **NPT:** Stabilize pressure/density (reach cruising altitude)
4. **Production:** Actual science happens (orbit achieved!)

Skip equilibration → Simulation explodes or gives garbage results

### Q2: What if simulation crashes?

**A:** Use checkpoint files (`.cpt`)
```bash
# Restart from last checkpoint
gmx mdrun -v -deffnm md -cpi md.cpt
```

Checkpoints save full state every 15 minutes, can resume exactly where you left off.

### Q3: How do I know if it's working?

**Check these:**
1. **Temperature:** Should fluctuate around 300 K (±5 K)
2. **Pressure:** Should fluctuate around 1 bar (±100 bar is fine)
3. **Energy:** Should be stable (no drifts or explosions)
4. **RMSD:** Protein shouldn't unfold (RMSD < 5 Å for stable proteins)
   - For AviTag: High RMSD is GOOD (flexible peptide)

### Q4: Why 2 femtosecond timesteps?

**A:** Fastest motions in proteins: bond vibrations (~10 fs period)

Nyquist sampling: Need 2× resolution to capture motion
- H-H bond vibrates at ~100 fs period
- 2 fs timesteps safely capture this
- Could use 1 fs (more accurate, 2× slower)
- Can't use 5 fs (misses fast vibrations, unstable)

### Q5: Why does this take so long?

**Calculating forces for 106,942 atoms:**
- Each atom interacts with ~500 neighbors (within cutoff)
- ~50 million force calculations per timestep
- 25 million timesteps for 50 ns
- Total: ~1.25 × 10¹⁵ force calculations!

CPUs do ~10⁹ calculations per second
→ Takes days even on fast computers

---

## Part 7: The End Result

### When Simulation Finishes

You'll have:

**1. Trajectory file (`md.xtc`)**
- 5,001 frames of all 106,942 atoms
- Can visualize in PyMOL (watch protein move!)
- Can calculate RMSD, Rg, hydrogen bonds, etc.

**2. Energy file (`md.edr`)**
- Temperature every 10 ps
- Pressure every 10 ps
- Total energy every 10 ps
- Use for quality control

**3. Analysis-ready data**
- Extract residues 135-149 (AviTag only)
- Calculate metrics for just the tag
- Compare to free peptide simulation

### Analysis Workflow

```
md.xtc (full trajectory)
   ↓
Extract residues 135-149 (AviTag)
   ↓
Calculate RMSD over time (flexibility)
Calculate Rg over time (compactness)
Cluster conformations (find representatives)
   ↓
Compare to free peptide:
  - Are RMSD distributions similar?
  - Are Rg distributions similar?
  - Do both explore similar conformational space?
   ↓
DECISION:
  Similar? → Use free ensemble for design
  Different? → Use fusion ensemble for design
```

---

## Summary: The Journey

```
START: LaG16-AviTag fusion.pdb (atomic coordinates)
         ↓
      Docker: Create portable Linux environment
         ↓
   GROMACS: Add physics (force field parameters)
         ↓
   GROMACS: Build simulation box with water & ions
         ↓
   GROMACS: Minimize energy (remove bad contacts)
         ↓
   GROMACS: Heat to 300 K (NVT equilibration)
         ↓
   GROMACS: Stabilize pressure (NPT equilibration)
         ↓
   GROMACS: Run production MD (50 ns of dynamics)
         ↓ (wait 2 days...)
         ↓
END: Trajectory file (movie of atomic motion)
         ↓
   Python: Extract AviTag residues, calculate metrics
         ↓
ANSWER: Does fusion change tag behavior? Design accordingly!
```

---

**Current Status:** NVT equilibration running in background (100 ps, ~5-10 min)

**Next:** NPT equilibration (100 ps), then production MD (50 ns, ~2 days)

**Questions?** Ask about any step - I can go deeper on any part!
