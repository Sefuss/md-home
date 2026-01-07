# Creating VisualHub-Quality Movies in ChimeraX

## The Problem

**Current movie**: Frame-by-frame jittery playback, wrong zoom level
**VisualHub standard**: Smooth, cinematic, professional-quality molecular dynamics visualization

## Why Gleb's Movies Look Better

TheVisualHub (Gleb Kuznetsov) uses several techniques:

### 1. Frame Interpolation
- **Current**: Playing raw trajectory frames (2,396 frames → 96 seconds at 25 fps)
- **VisualHub**: Interpolates between frames to create smoother motion
- **Fix**: Use `coordset #1 1,-1,1` with slower playback OR render at 60 fps and interpolate in post

### 2. Camera Work
- **Current**: Static camera, basic zoom
- **VisualHub**: Dynamic camera angles, smooth rotations, optimal framing
- **Fix**: Multiple camera positions, smooth transitions

### 3. Temporal Smoothing
- **Current**: No smoothing - shows every trajectory jitter
- **VisualHub**: Applies temporal averaging to reduce noise
- **Fix**: Pre-process trajectory with GROMACS `gmx filter` or smooth in ChimeraX

### 4. Rendering Quality
- **Current**: Standard rendering with 3x supersampling
- **VisualHub**: Ray-tracing or high-quality ambient occlusion
- **Fix**: Use ChimeraX's ray-tracing renderer

### 5. Visual Style
- **Current**: Basic cartoon + color
- **VisualHub**: Sophisticated lighting, shadows, depth of field
- **Fix**: Advanced lighting settings, silhouettes, shadows

---

## Solution: Multi-Step Workflow

### Step 1: Smooth the Trajectory (GROMACS)

```bash
# Average over 5-frame window for smoother motion
cd "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation"

docker run --rm -v "C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble:/workspace" \
  -w /workspace/01_md_preparation gromacs/gromacs:latest \
  gmx filter -f md_centered.xtc -s npt.gro -nf 5 -o md_smooth.xtc
```

This creates `md_smooth.xtc` with temporal averaging.

### Step 2: Updated ChimeraX Script (VisualHub-Style)

See `create_avitag_movie_VISUALHUB.py` for implementation.

Key changes:
1. **Zoom adjustment**: `zoom 1.2` instead of `0.85`
2. **Slower playback**: Play every 2nd frame for smoother appearance
3. **Ray-tracing**: Use high-quality renderer
4. **Advanced lighting**: Shadows and ambient occlusion
5. **Higher framerate**: 60 fps output for smoother playback

### Step 3: Optional Post-Processing (FFmpeg)

```bash
# Add motion blur in post-processing
ffmpeg -i input.mp4 -vf "minterpolate=fps=60:mi_mode=mci" output_smooth.mp4
```

---

## Gleb's Signature Techniques

From analyzing VisualFactory workflows:

### Lighting Setup
```python
run(session, "lighting soft")
run(session, "lighting shadows true intensity 0.6")
run(session, "lighting depthCue true start 0.5 end 1.0")
run(session, "lighting quality finer")
```

### Camera Motion
```python
# Smooth rotation during playback
run(session, "turn y 1 120")  # 1 degree per frame for 120 frames
run(session, "wait 120")
```

### Color Grading
```python
# Use color schemes that pop on white backgrounds
run(session, "color bychain palette Paired-12")  # High-contrast palette
run(session, "color byhetero")  # Highlight ligands/special residues
```

### Depth of Field (Advanced)
```python
# Requires ChimeraX 1.3+
run(session, "camera depthOfField true")
run(session, "camera focalPlane 0")
run(session, "camera fieldOfView 30")
```

---

## Recommended Settings for AviTag Movie

### Basic Quality (Quick)
- Smooth trajectory: `gmx filter -nf 5`
- ChimeraX: 3x supersampling, 30 fps
- Zoom: `1.2` to show full peptide dynamics

### Professional Quality (Slow)
- Smooth trajectory: `gmx filter -nf 10`
- ChimeraX: Ray-tracing renderer, 60 fps
- Add camera rotation: `turn y 0.5 frames`
- Post-process: FFmpeg interpolation

### VisualHub Quality (TheVisualHub Standard)
- Smooth trajectory: `gmx filter -nf 10`
- ChimeraX: Ray-tracing + shadows + depth cue
- Multiple camera angles (render 3 versions, composite)
- Color grading in DaVinci Resolve or similar
- 4K resolution (3840x2160)

---

## Quick Reference: ChimeraX Commands

```python
# Zoom adjustments
run(session, "zoom 0.8")   # Zoom out (see more)
run(session, "zoom 1.5")   # Zoom in (focus on region)

# Camera positioning
run(session, "view")                  # Reset to fit
run(session, "turn y 30")             # Rotate 30° around Y-axis
run(session, "move x 10")             # Pan 10 Å along X-axis
run(session, "cofr 0,0,0")            # Set center of rotation

# Playback speed
run(session, "coordset #1 1,-1,1")    # Every frame (default)
run(session, "coordset #1 1,-1,2")    # Every 2nd frame (slower, smoother)
run(session, "coordset #1 1,-1,5")    # Every 5th frame (much slower)

# Advanced rendering
run(session, "set silhouetteColor black")
run(session, "set silhouetteWidth 3")
run(session, "lighting shadows true")
run(session, "lighting depthCue true")
```

---

## Troubleshooting

### "Movie still looks jittery"
- Pre-smooth trajectory with `gmx filter -nf 10`
- Use `coordset #1 1,-1,2` (every 2nd frame)
- Post-process with FFmpeg motion interpolation

### "Peptide is cut off at edges"
- Increase zoom out: `zoom 1.5` or `zoom 2.0`
- Check trajectory centering: should use `md_centered.xtc`

### "Doesn't look like Gleb's work"
- Enable ray-tracing (slow but beautiful)
- Add shadows: `lighting shadows true`
- Use professional color palette
- Consider 4K resolution

---

## Next Steps

1. Create smoothed trajectory: `gmx filter`
2. Test new movie script with better zoom/lighting
3. Compare output to VisualHub examples
4. Iterate on camera angles and timing
5. Optional: Learn DaVinci Resolve for color grading
