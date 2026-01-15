# MD Simulation Execution Log
## LaG16-AviTag Fusion - 100ns Production Run

This log tracks which machine ran which trajectory segments for later analysis/concatenation.

---

## Trajectory Segments

### Part 0001 - Work PC
- **Machine:** Work PC (unknown GPU)
- **Date:** 2026-01-14
- **Steps:** 0 → 12,615,000
- **Time:** 0 ns → 25.23 ns
- **File:** `md.part0001.xtc` (not in git - stayed on work PC)
- **Status:** Completed

### Part 0002 - Home PC
- **Machine:** Home PC (GTX 1080 Ti, i7-7700K)
- **Date:** 2026-01-14 5:30 PM → 2026-01-15 (ongoing)
- **Steps:** 12,615,000 → 50,000,000 (target)
- **Time:** 25.23 ns → 100 ns (target)
- **File:** `md.part0002.xtc` (not in git - stays on home PC)
- **Status:** Running (~52% complete as of 5:30 AM Jan 15)
- **Performance:** ~102 ns/day

---

## Checkpoints in Git

The checkpoint files (`md.cpt`, `md_prev.cpt`) are committed to git and contain:
- Current simulation state
- Exact step number for restart
- All velocities, positions, forces

**Latest checkpoint:** Step 26,120,000 (52.24 ns) as of 2026-01-15 5:30 AM

---

## For Complete Trajectory Analysis

When analysis requires the full 0-100 ns trajectory, concatenate parts:

```bash
# Copy trajectory files to same location (if on different machines)
# Then combine:
gmx trjcat -f md.part0001.xtc md.part0002.xtc -o md_complete.xtc

# Or analyze individual segments separately
```

---

## Notes

- Trajectory files (.xtc) are excluded from git (too large, ~1 GB each)
- Each machine keeps its trajectory segment locally
- Checkpoint files in git allow seamless continuation across machines
- Part numbers increment automatically when using `-noappend` flag
