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
- **File:** `md.part0002.xtc` (not in git - stays on home PC, currently 1.9 GB)
- **Status:** Running (75% complete as of 4:00 PM Jan 15, at 74.82 ns)
- **Performance:** ~54 ns/day sustained
- **Expected completion:** ~3:00 AM Jan 16

---

## Checkpoints in Git

The checkpoint files (`md.cpt`, `md_prev.cpt`) are committed to git and contain:
- Current simulation state
- Exact step number for restart
- All velocities, positions, forces

**Latest checkpoint:** Step 37,410,000 (74.82 ns) as of 2026-01-15 4:00 PM

---

## For Complete Trajectory Analysis (0-100 ns)

### IMPORTANT: Retrieve md.part0001.xtc from Work PC!

To get the **full 0-100 ns trajectory** for movie/analysis, you need both parts:
- **Part 0001 (Work PC):** 0 → 25.23 ns - **RETRIEVE THIS AT WORK!**
- **Part 0002 (Home PC):** 25.23 → 100 ns - Already have it

### Steps to Create Complete Trajectory:

1. **At work PC tomorrow:** Copy `md.part0001.xtc` from the MD folder to USB/cloud
2. **At home:** Place both files in same directory
3. **Concatenate:**
   ```bash
   gmx trjcat -f md.part0001.xtc md.part0002.xtc -o md_complete_0-100ns.xtc
   ```

**Alternative:** If you only need 25-100 ns analysis, use `md.part0002.xtc` alone (75 ns of data)

---

## Notes

- Trajectory files (.xtc) are excluded from git (too large, ~1 GB each)
- Each machine keeps its trajectory segment locally
- Checkpoint files in git allow seamless continuation across machines
- Part numbers increment automatically when using `-noappend` flag
