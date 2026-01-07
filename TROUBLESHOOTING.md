# Troubleshooting Guide

**Common issues and solutions encountered during project**

**Last Updated:** January 2, 2026

---

## Visualization Issues

### Issue 1: Peptide Wraps Around Screen (PBC Artifacts)

**Symptom:**
- When viewing trajectory in ChimeraX, peptide stretches unnaturally
- Molecule appears to "teleport" to opposite side of screen
- Bonds look broken across the visualization

**Cause:**
- MD simulations use periodic boundary conditions (PBC)
- This creates infinite space by wrapping at box edges
- Perfect for physics, terrible for visualization

**Solution:**

**Step 1:** Pre-process trajectory with GROMACS `trjconv`

```bash
cd "01_md_preparation"

# Create PBC-corrected trajectory
echo "1 0" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx trjconv -s md.tpr -f md.xtc -o md_centered.xtc \
  -center -pbc mol -ur compact
```

**What this does:**
- `-center` → Centers protein in simulation box
- `-pbc mol` → Makes molecules "whole" (no broken bonds across boundaries)
- `-ur compact` → Uses compact representation

**Step 2:** Use cleaned trajectory in ChimeraX

```
Drag: visualization/chimerax_load_free_peptide_CLEAN.py
```

This script automatically loads `md_centered.xtc` (the corrected version)

**Files created:**
- `01_md_preparation/md_centered.xtc` - Cleaned trajectory
- `visualization/chimerax_load_free_peptide_CLEAN.py` - Loads cleaned version

**Status:** FIXED ✅ (Jan 2, 2026)

---

### Issue 2: ChimeraX Syntax Error - "Expected a number or a keyword"

**Symptom:**
```
chimerax.core.errors.UserError: Expected a number or a keyword
Error in line: graphics quality higher
```

**Cause:**
- Invalid ChimeraX command syntax
- `graphics quality higher` does not exist in ChimeraX
- Attempted to set rendering quality, but used wrong syntax

**Solution:**

Remove the problematic line from scripts:

**Before (WRONG):**
```python
run(session, "graphics quality higher")  # ❌ Invalid syntax
```

**After (CORRECT):**
```python
# Line removed - default quality is sufficient
# (Alternatively could use proper ChimeraX quality settings if needed)
```

**Files fixed:**
- `visualization/chimerax_load_free_peptide_CLEAN.py` ✅
- `visualization/create_avitag_movie.py` ✅
- `visualization/chimerax_export_movie.py` ✅

**Status:** FIXED ✅ (Jan 2, 2026)

**Note:** If high quality rendering is needed, use:
```python
run(session, "set maxFrameRate 60")  # Smoother playback
run(session, "graphics silhouettes true")  # Better depth perception
```

---

### Issue 3: ChimeraX coordset Error - "Expected a keyword"

**Symptom:**
```
chimerax.core.errors.UserError: Expected a keyword
Error in line: coordset #1 1,end
```

**Cause:**
- Invalid ChimeraX coordset syntax
- `coordset #1 1,end` is not valid (ChimeraX doesn't recognize "end")
- Correct syntax uses asterisk (*) for all frames

**Solution:**

Replace `1,end` with `*`:

**Before (WRONG):**
```python
run(session, "coordset #1 1,end")  # ❌ Invalid - "end" not recognized
```

**After (CORRECT):**
```python
run(session, "coordset #1 *")  # ✅ Asterisk means "all frames"
```

**Files fixed:**
- `visualization/create_avitag_movie.py` ✅
- `visualization/chimerax_export_movie.py` ✅
- `visualization/chimerax_load_free_peptide.py` ✅
- `visualization/chimerax_load_free_peptide_CLEAN.py` ✅
- `visualization/chimerax_side_by_side_comparison.py` ✅

**Alternative syntax:**
```python
coordset #1 1,2396  # Play frames 1 to 2396 (specific)
coordset #1 1,100,2 # Play frames 1 to 100, every 2nd frame
coordset #1 *       # Play all frames (recommended)
```

**Status:** FIXED ✅ (Jan 2, 2026)

---

### Issue 4: ChimeraX "Already recording a movie"

**Symptom:**
```
chimerax.core.errors.UserError: Already recording a movie
```

**Cause:**
- Previous script run started recording but crashed
- Recording session never stopped
- ChimeraX only allows one recording at a time

**Solution:**

**Option 1: Stop recording manually**
```python
# In ChimeraX command line:
movie stop
movie reset
```

**Option 2: Restart ChimeraX**
- Close ChimeraX completely
- Reopen
- Try script again

**Option 3: Use updated scripts (automatic fix)**
- All scripts now auto-clear previous recordings
- Try-except blocks handle this gracefully

**Prevention:**
- Updated scripts now include:
```python
try:
    run(session, "movie stop")
    run(session, "movie reset")
except:
    pass  # No previous recording
```

**Files updated:**
- `visualization/create_avitag_movie.py` ✅
- `visualization/chimerax_export_movie.py` ✅

**Status:** FIXED ✅ (Jan 2, 2026)

---

## Python Script Issues

### Issue 3: UnicodeEncodeError with Checkmark Characters

**Symptom:**
```
UnicodeEncodeError: 'charmap' codec can't encode character '\u2713' in position 2
```

**Cause:**
- Used unicode characters (✓ checkmarks) in print statements
- Windows console (cmd.exe) uses cp1252 encoding
- Doesn't support many unicode characters

**Solution:**

**Option 1:** Remove unicode characters (USED)
```python
# Before
print("✅ Saved: file.png")  # ❌ Fails on Windows

# After
print("Saved: file.png")     # ✅ Works everywhere
```

**Option 2:** Set UTF-8 encoding
```bash
export PYTHONIOENCODING=utf-8
python script.py
```

**Files fixed:**
- `analysis/create_powerpoint_plots.py` - Rewrote inline without unicode

**Status:** FIXED ✅ (Dec 18, 2025)

---

## Docker Issues

### Issue 4: Docker Path Expansion with $(pwd)

**Symptom:**
```
docker: Error response from daemon: create $(pwd): "$(pwd)" includes invalid characters
```

**Cause:**
- Using `$(pwd)` in Git Bash on Windows
- Docker on Windows doesn't expand shell variables properly
- Needs absolute Windows paths

**Solution:**

Use absolute paths instead:

**Before (WRONG):**
```bash
docker run -v "$(pwd):/work" ...  # ❌ Doesn't work on Windows
```

**After (CORRECT):**
```bash
docker run -v "C:\Users\bmartinez\...\01_md_preparation:/work" ...  # ✅ Works
```

**Or use environment variable:**
```bash
export MSYS_NO_PATHCONV=1  # Disable Git Bash path conversion
docker run -v "$(pwd):/work" ...
```

**Status:** FIXED ✅ (Dec 18, 2025)

---

### Issue 5: Python Command Not Found

**Symptom:**
```
Python was not found
```

**Cause:**
- Python not in system PATH
- Multiple Python installations
- Using `python` instead of full path

**Solution:**

Use full Anaconda Python path:

```bash
"C:\Users\bmartinez\AppData\Local\anaconda3\python.exe" script.py
```

**Or add to PATH (permanent fix):**
1. Open: System Properties → Environment Variables
2. Edit PATH
3. Add: `C:\Users\bmartinez\AppData\Local\anaconda3`
4. Restart terminal

**Status:** WORKAROUND ✅ (Using full path in documentation)

---

## MD Simulation Issues

### Issue 6: Container Not Found

**Symptom:**
```
docker: Error response from daemon: No such container: fusion_md_production
```

**Cause:**
- Container was stopped or removed
- Wrong container name

**Solution:**

```bash
# List all containers
docker ps -a

# Check if container exists but stopped
docker ps -a | grep fusion

# Restart if stopped
docker start fusion_md_production

# If removed, need to re-launch MD
```

**Prevention:**
- Use `--name` flag when creating containers
- Don't use `--rm` for long-running simulations
- Document container names

**Status:** DOCUMENTED ✅

---

## Quick Reference

| Issue | Quick Fix | Documentation |
|-------|-----------|---------------|
| PBC wrapping | Use `md_centered.xtc` | This file, Issue #1 |
| ChimeraX graphics quality | Remove invalid command | This file, Issue #2 |
| ChimeraX coordset | Use `coordset #1 *` | This file, Issue #3 |
| Already recording movie | `movie stop; movie reset` | This file, Issue #4 |
| Unicode errors | Remove unicode chars | This file, Issue #5 |
| Docker paths | Use absolute paths | This file, Issue #6 |
| Python not found | Use full path | This file, Issue #7 |
| Container missing | Check `docker ps -a` | This file, Issue #8 |

---

## Reporting New Issues

**When encountering new issues:**

1. **Document the error:**
   - Full error message
   - Command that caused it
   - Environment (Windows/Linux, software versions)

2. **Try basic debugging:**
   - Check file paths
   - Verify containers running (`docker ps`)
   - Check file permissions

3. **Add to this file:**
   - Clear symptom description
   - Root cause
   - Solution that worked
   - Files modified

---

## Prevention Tips

**Avoid common pitfalls:**

1. **Always read files first** before editing (use Read tool)
2. **Use absolute paths** on Windows with Docker
3. **Check ChimeraX commands** against official docs
4. **Test scripts** on small datasets first
5. **Keep backups** of working configurations
6. **Document changes** immediately

---

## External Resources

**ChimeraX:**
- Command reference: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/
- Troubleshooting: https://www.rbvi.ucsf.edu/chimerax/docs/user/troubleshooting.html

**GROMACS:**
- Manual: https://manual.gromacs.org/
- Common issues: https://manual.gromacs.org/documentation/current/user-guide/faq.html

**Docker:**
- Windows guide: https://docs.docker.com/desktop/windows/
- Troubleshooting: https://docs.docker.com/desktop/troubleshoot/overview/

---

**Last Updated:** January 2, 2026
**Issues Logged:** 8
**Issues Resolved:** 8
