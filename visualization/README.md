# Visualization Scripts (ChimeraX)

**Purpose:** Scripts for visualizing MD trajectories and creating presentation movies

**Software Required:** UCSF ChimeraX (installed)

---

## Quick Start

All scripts are **drag-and-drop** - just open ChimeraX and drag the `.py` file into the window!

---

## Scripts

### ⭐⭐⭐ create_avitag_movie.py (RECOMMENDED - ONE-CLICK)
**Use Case:** Create professional movie automatically

**What it does:**
- Loads PBC-corrected trajectory
- Applies professional styling (B-factor coloring, lighting, silhouettes)
- Sets HD resolution (1920x1080)
- Records and encodes movie automatically
- Output: `C:\Users\bmartinez\Desktop\AviTag_Dynamic_Ensemble.mp4`

**How to use:**
1. Open ChimeraX
2. Drag this file into ChimeraX window
3. Wait ~5 minutes
4. Done! Movie is on your Desktop

**Settings:**
- Resolution: 1920x1080 (Full HD)
- Framerate: 25 fps
- Quality: High (H.264)
- Supersampling: 3x
- Duration: ~96 seconds

---

### ⭐⭐ chimerax_load_free_peptide_CLEAN.py (RECOMMENDED)
**Use Case:** Load trajectory for manual exploration/recording

**What it does:**
- Loads **PBC-corrected** trajectory (no wrapping artifacts!)
- Applies professional styling
- Sets up optimal camera angle
- White background (PowerPoint-ready)

**How to use:**
1. Open ChimeraX
2. Drag this file into window
3. Press Spacebar to play
4. Use mouse to rotate/zoom
5. Manually record if desired

**Advantages:**
- No periodic boundary artifacts
- Professional styling pre-configured
- Ready for presentation

---

### chimerax_load_free_peptide.py
**Use Case:** Quick viewing of original trajectory

**What it does:**
- Loads original trajectory (md.xtc)
- Basic styling
- Good for quick checks

**Note:** May show PBC wrapping artifacts. Use CLEAN version for presentations!

---

### chimerax_export_movie.py
**Use Case:** Manual movie export with custom settings

**What it does:**
- Exports currently loaded trajectory as MP4
- High-quality settings pre-configured
- For advanced users who want manual control

**How to use:**
1. First load trajectory with one of the load scripts above
2. Arrange camera angle as desired
3. Drag this script to export

**Settings:**
- 1920x1080 resolution
- 3x supersampling
- H.264 encoding

---

### chimerax_side_by_side_comparison.py
**Use Case:** Compare free vs fusion trajectories (AFTER fusion MD completes)

**What it does:**
- Loads both free peptide and fusion protein trajectories
- Sets up split-screen view
- Synchronized playback
- Identical styling for fair comparison

**How to use:**
1. Wait for fusion MD to complete (~16 days from Jan 2, 2026)
2. Open ChimeraX
3. Drag this file
4. Watch side-by-side comparison!

**Output idea:**
- Record comparison movie
- Shows how protein context affects peptide dynamics
- Perfect for presentations explaining the design rationale

---

## Tips for Manual Movie Recording

After loading with `chimerax_load_free_peptide_CLEAN.py`:

```python
# Set resolution
windowsize 1920 1080

# Start recording
movie record supersample 3

# Play trajectory
coordset #1 1,end

# Wait for completion
wait 100

# Encode
movie encode C:\Users\bmartinez\Desktop\my_movie.mp4 format h264 quality high framerate 25
```

---

## Apply VisualFactory Smoothing (Optional)

For even smoother, more cinematic motion:

**After loading trajectory:**
1. Drag: `C:\Users\bmartinez\ProteinDesign\VisualFactory\UltimateSmoothMD\UltimateSmoothMD5.py`
2. Script will smooth the trajectory
3. Then record movie

**Result:** Removes thermal jitter, creates fluid motion

---

## Troubleshooting

**Issue: Molecule wraps around screen**
- Use `chimerax_load_free_peptide_CLEAN.py` (PBC-corrected!)
- Not the regular version

**Issue: Movie is jittery**
- Apply VisualFactory smoothing (see above)
- Increase supersampling: `movie record supersample 4`
- Play fewer frames: `coordset #1 1,end,2`

**Issue: Can't see details**
- Zoom in: `zoom 1.5`
- Increase quality: `graphics quality higher`

**Issue: File too large**
- Reduce resolution: `windowsize 1280 720`
- Lower framerate: `framerate 15`

---

## File Paths Reference

**Free peptide trajectory (PBC-corrected):**
```
Structure:  C:\Users\bmartinez\...\01_md_preparation\npt.gro
Trajectory: C:\Users\bmartinez\...\01_md_preparation\md_centered.xtc
```

**Fusion protein trajectory (when ready):**
```
Structure:  C:\Users\bmartinez\...\03_fusion_md\npt.gro
Trajectory: C:\Users\bmartinez\...\03_fusion_md\md.xtc
```

---

## Further Reading

- **Complete guide:** `../HOLLYWOOD_QUALITY_MOVIE_GUIDE.md`
- **ChimeraX setup:** `../CHIMERAX_SETUP_GUIDE.md`
- **ChimeraX docs:** https://www.rbvi.ucsf.edu/chimerax/docs/user/index.html

---

**Last Updated:** January 2, 2026
