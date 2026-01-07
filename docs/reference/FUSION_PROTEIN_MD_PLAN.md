# Fusion Protein MD Plan: Testing Context-Dependent Behavior

**Date:** January 2, 2026
**Purpose:** Address fundamental design question - Does AviTag behavior change when attached to a protein?
**Status:** Planning phase

---

## The Critical Question

**User's insight:** "If I take some valuable frames and design binders for it, won't that only bind those poses and not the actual peptide shape?"

**The Gap Identified:**
- Free peptide MD samples unbound conformations
- When attached to a protein, the tag may behave differently:
  - Constrained by fusion partner
  - Stabilized in specific conformations
  - Different conformational preferences
- Binders designed to free ensemble may not work if bound ensemble is different

**The Test:** Run MD on nanobody-AviTag fusion and compare tag behavior to free peptide

---

## Experimental Design

### System 1: Free AviTag (COMPLETED ✅)
- **What:** Isolated 15-residue AviTag peptide in water
- **Duration:** 23 ns
- **Result:** Extensive conformational sampling (RMSD 8.6 Å, Rg 12.5 Å)
- **Conformations:** Dynamic equilibrium between compact and extended states
- **Data:** `01_md_preparation/md.xtc` (2,304 frames)

### System 2: Nanobody-AviTag Fusion (PLANNED)
- **What:** LaG16 nanobody with C-terminal AviTag (GLNDIFEAQKIEWHE)
- **Duration:** 30-50 ns (smaller than 100 ns, since we just need to characterize tag)
- **Analysis:** Track ONLY the AviTag residues (same as System 1)
- **Comparison metrics:**
  - RMSD distribution
  - Rg distribution
  - Conformational space overlap
  - Flexibility (RMSF per residue)

---

## Comparison Metrics

### Question 1: Is the tag more constrained when attached?
**Metric:** RMSD range and standard deviation
- **Free peptide:** RMSD 0-25 Å (high flexibility)
- **Fusion:** RMSD X-Y Å (to be measured)
- **If similar:** Tag remains flexible regardless of fusion
- **If lower:** Tag is constrained by nanobody

### Question 2: Are the conformations similar or different?
**Metric:** Conformational space overlap
- Extract cluster representatives from both simulations
- Calculate RMSD between free and fusion clusters
- **If <3 Å overlap:** Conformations are similar → Free ensemble is valid
- **If >5 Å separation:** Conformations are different → Need fusion ensemble

### Question 3: Does the tag prefer specific states when attached?
**Metric:** Radius of gyration distribution
- **Free peptide:** Rg 7.7-18.5 Å (compact ↔ extended equilibrium)
- **Fusion:** Rg distribution (to be measured)
- **If similar distribution:** No preference, same equilibrium
- **If shifted:** Fusion stabilizes specific state(s)

---

## Decision Matrix

Based on comparison results:

### Scenario A: Free and Fusion Are Similar
**Finding:** Tag behavior is similar whether free or attached
**Interpretation:** Tag is intrinsically flexible, fusion partner doesn't constrain it
**Action:** Proceed with FREE peptide ensemble for binder design
**Rationale:** Free ensemble captures all relevant conformations

### Scenario B: Fusion Is More Constrained
**Finding:** Tag has narrower conformational distribution when attached
**Interpretation:** Fusion partner stabilizes specific conformations
**Action:** Use FUSION ensemble for binder design
**Rationale:** These are the conformations present in real applications

### Scenario C: Fusion Has Different Conformations
**Finding:** Tag adopts distinct conformations not seen in free peptide
**Interpretation:** Fusion partner influences tag structure
**Action:** Use FUSION ensemble AND test with multiple fusion proteins
**Rationale:** Need to capture context-dependent behavior

---

## Building the Fusion Model

### Step 1: Choose Fusion Partner
**Primary choice:** LaG16 nanobody (from Approach B)
- **Why:** We already have the structure
- **AviTag position:** C-terminus (typical biotinylation strategy)
- **Linker:** May need short flexible linker (e.g., GGS or GGGGS)

**Alternative choices:**
- VHH72 (another nanobody)
- Simplified scaffold (just a stable protein domain)

### Step 2: Model Construction
**Option A: Manual PDB editing**
```
LaG16 structure (chain A, residues 1-120)
        ↓
Add C-terminal AviTag (residues 121-135)
        ↓
GLNDIFEAQKIEWHE appended to C-term
```

**Option B: ColabFold prediction**
- Input: LaG16 sequence + AviTag sequence (concatenated)
- Predict full fusion structure
- Advantage: Includes any potential interactions between domains
- Disadvantage: May over-predict interactions (AlphaFold bias toward contacts)

**Option C: PyMOL modeling**
- Load LaG16 structure
- Build AviTag from C-terminus
- Energy minimize to remove clashes
- Advantage: Fast, controlled
- Disadvantage: Initial conformation is arbitrary

### Step 3: System Preparation (Same as Free Peptide)
1. **pdb2gmx:** Generate topology (AMBER99SB-ILDN)
2. **editconf:** Create simulation box (1.2 nm from protein)
3. **solvate:** Add water (TIP3P)
4. **genion:** Neutralize with Na+/Cl-
5. **Energy minimization:** Remove bad contacts
6. **NVT equilibration:** 100 ps at 300 K
7. **NPT equilibration:** 100 ps at 1 bar, 300 K
8. **Production MD:** 30-50 ns

---

## Analysis Protocol

### Step 1: Load Both Trajectories
```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np

# Free peptide
u_free = mda.Universe('free_peptide/npt.gro', 'free_peptide/md.xtc')
peptide_free = u_free.select_atoms('protein')

# Fusion (select ONLY AviTag residues)
u_fusion = mda.Universe('fusion/npt.gro', 'fusion/md.xtc')
# If AviTag is residues 121-135 in fusion
avitag_fusion = u_fusion.select_atoms('resid 121-135')
```

### Step 2: Calculate Metrics for Both
```python
# Align each trajectory independently
align.AlignTraj(u_free, u_free, select='backbone', in_memory=True).run()
align.AlignTraj(u_fusion, u_fusion, select='resid 121-135 and backbone',
                in_memory=True).run()

# RMSD for both
rmsd_free = rms.RMSD(u_free, select='backbone').run().results.rmsd[:, 2]
rmsd_fusion = rms.RMSD(u_fusion, select='resid 121-135 and backbone').run().results.rmsd[:, 2]

# Rg for both
rg_free = [peptide_free.radius_of_gyration() for ts in u_free.trajectory]
rg_fusion = [avitag_fusion.radius_of_gyration() for ts in u_fusion.trajectory]
```

### Step 3: Statistical Comparison
```python
from scipy import stats

# Compare distributions
ks_stat_rmsd, p_value_rmsd = stats.ks_2samp(rmsd_free, rmsd_fusion)
ks_stat_rg, p_value_rg = stats.ks_2samp(rg_free, rg_fusion)

print(f"RMSD distributions: KS stat = {ks_stat_rmsd:.3f}, p = {p_value_rmsd:.3e}")
print(f"Rg distributions: KS stat = {ks_stat_rg:.3f}, p = {p_value_rg:.3e}")

# If p > 0.05: Distributions are NOT significantly different (free = fusion)
# If p < 0.05: Distributions ARE significantly different (free ≠ fusion)
```

### Step 4: Visual Comparison
Create overlay plots:
- RMSD time series (free vs fusion)
- Rg time series (free vs fusion)
- Conformational space (RMSD vs Rg scatter, colored by source)
- Distribution histograms (overlaid)

---

## Timeline Estimate

| Task | Duration | Cumulative |
|------|----------|------------|
| Build fusion model | 1-2 hours | 2 hours |
| Prepare MD system | 1 hour | 3 hours |
| Run equilibration | 15 min | 3.25 hours |
| Run production (30 ns) | 24-36 hours* | 27-39 hours |
| Analysis & comparison | 2-3 hours | 30-42 hours |

*At 27 ns/day (CPU performance from previous run)
*Could be faster if using GPU: ~2-4 hours for 30 ns

**Total:** ~2 days wall time (mostly waiting for MD)

---

## Expected Outcomes

### Most Likely Outcome (Based on Literature)
**Prediction:** Tag remains flexible in fusion context

**Reasoning:**
- AviTag is intrinsically disordered
- C-terminal tags have minimal constraints from protein
- Literature shows C-term tags generally maintain flexibility

**Result:** Free ensemble is representative → Proceed with free ensemble design

### Alternative Outcome
**Prediction:** Tag is partially constrained by fusion

**Reasoning:**
- Some transient interactions with nanobody surface
- Reduced conformational space

**Result:** Use fusion ensemble → More realistic for applications

---

## Addressing the Broader Design Question

Even if free and fusion ensembles are similar, your concern about **induced fit** upon binder binding remains valid.

### The Conformational Selection vs Induced Fit Problem

**Conformational Selection (ensemble approach works):**
```
Free peptide: [A, B, C, D, E] ← Multiple conformations
      ↓
Binder binds: Selects conformation C
      ↓
Result: C already existed in free ensemble ✅
```

**Induced Fit (ensemble approach fails):**
```
Free peptide: [A, B, C, D, E]
      ↓
Binder approaches: Induces new conformation F
      ↓
Result: F was NOT in free ensemble ❌
```

### Mitigation Strategy: Post-Design MD Validation

**After designing binders to ensemble:**
1. Select top 50 designs (by pLDDT, pAE)
2. For each design, build binder-peptide complex
3. Run short MD (10-20 ns) of the complex
4. Check if:
   - Peptide maintains favorable conformation
   - Interface remains stable
   - Binding energy is favorable (MM-PBSA)
5. Filter designs where peptide undergoes large unfavorable changes

**This catches induced fit problems before experiments!**

---

## Recommended Next Steps

### Immediate (Today):
1. ✅ Document the design approach (this file)
2. Build LaG16-AviTag fusion model
3. Prepare fusion MD system

### Short-term (This week):
1. Run fusion protein MD (30-50 ns)
2. Analyze and compare to free peptide
3. Make decision on which ensemble to use

### Medium-term (Next week):
1. Cluster selected ensemble
2. Design binders with RFD3
3. Add MD validation step for top designs

---

## Success Criteria

By the end of this fusion protein MD experiment, we will know:

✅ Whether free peptide behavior is representative of fusion context
✅ If tag conformations change when attached to protein
✅ Which ensemble to use for binder design
✅ Whether our approach is scientifically sound

**This directly addresses your identified gap!**

---

## Files to Create

### Fusion Model Directory Structure
```
02_scaffolds/approach_C_peptide_ensemble/
  ├── 00_initial_structure/
  │   └── avitag.pdb (free peptide) ✅
  ├── 01_md_preparation/
  │   └── md.xtc (free peptide trajectory, 23 ns) ✅
  ├── 02_fusion_model/
  │   ├── lag16_avitag_fusion.pdb (model to build)
  │   ├── lag16_avitag_fusion.fasta (sequence)
  │   └── notes.md (modeling decisions)
  └── 03_fusion_md/
      ├── system preparation files (gro, top)
      ├── mdp files (minim, nvt, npt, md)
      └── md.xtc (fusion trajectory, 30-50 ns)
```

---

**Ready to build the fusion model?** This experiment will definitively answer whether your concerns about free vs bound conformations are justified for this specific system.
