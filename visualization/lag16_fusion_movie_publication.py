#!/usr/bin/env python3
"""
LaG16-AviTag Fusion MD Movie - Publication Quality
===================================================

Pre-configured for the LaG16-AviTag fusion MD trajectory
with publication-quality transparent surface rendering.

Usage:
    "C:\Program Files\ChimeraX 1.11\bin\chimerax.exe" --script lag16_fusion_movie_publication.py
"""

from chimerax.core.commands import run

# ============================================================================
# CONFIGURATION - Optimized for LaG16-AviTag Fusion
# ============================================================================

# Input files (LaG16-AviTag fusion MD)
STRUCTURE_FILE = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\protein_only.gro"
TRAJECTORY_FILE = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\protein_only.xtc"

# Output settings - Publication quality
OUTPUT_FILE = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\04_movies\lag16_fusion_publication_HQ.mp4"
RESOLUTION = (1920, 1080)  # Full HD
FRAME_RATE = 25
SUPERSAMPLE = 3  # High quality anti-aliasing

# Visual style - Professional teal (popular in Nature/Science papers)
SURFACE_COLOR = "#4ECDC4"  # Teal - clean, professional, publication-quality
SURFACE_OPACITY = 0.65  # 65% transparency
CARTOON_COLOR = "byhetero"  # Rainbow by chain
BACKGROUND_COLOR = "white"

# Alternative publication colors (uncomment to try):
# SURFACE_COLOR = "#5DADE2"  # Light blue - structural biology favorite
# SURFACE_COLOR = "#FA8072"  # Salmon - warm, good contrast
# SURFACE_COLOR = "#9B59B6"  # Purple - elegant, distinctive
# SURFACE_COLOR = "#6C7A89"  # Slate gray - sophisticated, neutral

# Rendering quality - Publication grade
LIGHTING = "soft"
SILHOUETTES = True
SILHOUETTE_WIDTH = 2
AMBIENT_OCCLUSION = True

# Camera settings
CAMERA_ROTATE = False  # Show dynamics without rotation
INITIAL_CAMERA_ANGLE = (30, 20, 0)
ZOOM_FACTOR = 0.85

# Trajectory playback - All 267 frames
PLAY_ALL_FRAMES = True
START_FRAME = 1
END_FRAME = -1
FRAME_STEP = 1

# Title overlay
ADD_TITLE = True
TITLE_TEXT = "LaG16-AviTag Fusion: Conformational Dynamics"
TITLE_POSITION = (0.5, 0.95)
TITLE_SIZE = 20
TITLE_COLOR = "black"

# ============================================================================
# SCRIPT - Same as trajectory_movie_maker.py main function
# ============================================================================

def create_md_movie(session):
    """Create publication-quality MD trajectory movie"""

    print("=" * 70)
    print("CHIMERAX MD TRAJECTORY MOVIE MAKER")
    print("LaG16-AviTag Fusion - Publication Quality")
    print("=" * 70)

    # Step 1: Load structure and trajectory
    print("\nStep 1: Loading structure and trajectory...")

    import os
    if not os.path.exists(STRUCTURE_FILE):
        print(f"ERROR: Structure file not found: {STRUCTURE_FILE}")
        return

    if not os.path.exists(TRAJECTORY_FILE):
        print(f"ERROR: Trajectory file not found: {TRAJECTORY_FILE}")
        return

    try:
        run(session, f"open {STRUCTURE_FILE}")
        print(f"  Loaded: {os.path.basename(STRUCTURE_FILE)}")

        run(session, f"open {TRAJECTORY_FILE} structureModel #1")
        print(f"  Loaded: {os.path.basename(TRAJECTORY_FILE)}")

        from chimerax.core.models import Model
        models = session.models.list()
        coords = [m for m in models if hasattr(m, 'coordset_ids')]
        if coords:
            num_frames = len(coords[0].coordset_ids)
            print(f"  Frames detected: {num_frames}")
        else:
            num_frames = None

    except Exception as e:
        print(f"ERROR loading files: {e}")
        return

    # Step 2: Set up visualization style
    print("\nStep 2: Configuring visualization style...")

    run(session, "hide atoms")
    run(session, "hide bonds")
    run(session, "cartoon")
    print("  Cartoon representation: ON")

    if CARTOON_COLOR == "byhetero":
        run(session, "color byhetero")
        print("  Cartoon color: By chain (rainbow)")
    elif CARTOON_COLOR == "bfactor":
        run(session, "color bfactor palette rainbow")
        print("  Cartoon color: By B-factor")
    else:
        run(session, f"color {CARTOON_COLOR}")
        print(f"  Cartoon color: {CARTOON_COLOR}")

    run(session, f"surface #1")
    run(session, f"color {SURFACE_COLOR} #1 target s")
    run(session, f"transparency {SURFACE_OPACITY * 100} #1 target s")
    print(f"  Surface: {SURFACE_COLOR} @ {SURFACE_OPACITY*100:.0f}% opacity")

    run(session, "surface style #1 mesh false")
    print("  Surface smoothing: ON")

    # Step 3: Configure lighting and rendering quality
    print("\nStep 3: Setting up lighting and rendering...")

    run(session, f"lighting {LIGHTING}")
    print(f"  Lighting: {LIGHTING}")

    if SILHOUETTES:
        run(session, "graphics silhouettes true")
        run(session, f"graphics silhouettes width {SILHOUETTE_WIDTH}")
        print(f"  Silhouettes: ON (width={SILHOUETTE_WIDTH})")
    else:
        run(session, "graphics silhouettes false")

    if AMBIENT_OCCLUSION:
        run(session, "lighting shadows true")
        run(session, "lighting quality finer")
        print("  Ambient occlusion: ON")

    run(session, f"set bgColor {BACKGROUND_COLOR}")
    print(f"  Background: {BACKGROUND_COLOR}")

    # Step 4: Position camera
    print("\nStep 4: Positioning camera...")

    run(session, "view")

    y_rot, x_rot, z_rot = INITIAL_CAMERA_ANGLE
    if y_rot != 0:
        run(session, f"turn y {y_rot}")
    if x_rot != 0:
        run(session, f"turn x {x_rot}")
    if z_rot != 0:
        run(session, f"turn z {z_rot}")

    run(session, f"zoom {ZOOM_FACTOR}")
    print(f"  Camera angle: {INITIAL_CAMERA_ANGLE}")
    print(f"  Zoom: {ZOOM_FACTOR}")

    # Step 5: Set window size
    print("\nStep 5: Setting resolution...")
    width, height = RESOLUTION
    run(session, f"windowsize {width} {height}")
    print(f"  Resolution: {width}x{height}")

    # Step 6: Add title
    if ADD_TITLE:
        print("\nStep 6: Adding title overlay...")
        x_pos, y_pos = TITLE_POSITION
        run(session, f'2dlabels create title text "{TITLE_TEXT}" '
                     f'xpos {x_pos} ypos {y_pos} color {TITLE_COLOR} size {TITLE_SIZE}')
        print(f'  Title: "{TITLE_TEXT}"')

    # Step 7: Prepare movie recording
    print("\n" + "=" * 70)
    print("MOVIE RECORDING")
    print("=" * 70)

    try:
        run(session, "movie stop")
        run(session, "movie reset")
    except:
        pass

    print(f"\nStarting recording (supersample {SUPERSAMPLE}x)...")
    run(session, f"movie record supersample {SUPERSAMPLE}")

    # Step 8: Play trajectory
    print("Playing trajectory...")

    if PLAY_ALL_FRAMES:
        frame_spec = "1,-1"
        if num_frames:
            print(f"  Frames: 1 to {num_frames} (all)")
    else:
        frame_spec = f"{START_FRAME},{END_FRAME},{FRAME_STEP}"
        print(f"  Frames: {START_FRAME} to {END_FRAME} (step={FRAME_STEP})")

    run(session, f"coordset #1 {frame_spec}")

    if num_frames and PLAY_ALL_FRAMES:
        wait_seconds = int(num_frames / FRAME_RATE * 1.2)
    else:
        wait_seconds = 100

    print(f"  Waiting for playback to complete (~{wait_seconds}s)...")
    run(session, f"wait {wait_seconds}")

    if CAMERA_ROTATE:
        print(f"\nAdding {ROTATION_DEGREES}° rotation around {ROTATION_AXIS}-axis...")
        run(session, f"turn {ROTATION_AXIS} {ROTATION_DEGREES} {wait_seconds}")
        run(session, f"wait {wait_seconds}")

    # Step 9: Encode movie
    print("\n" + "=" * 70)
    print("ENCODING MOVIE")
    print("=" * 70)

    output_dir = os.path.dirname(OUTPUT_FILE)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    print(f"\nEncoding to: {OUTPUT_FILE}")
    print(f"  Format: H.264 MP4")
    print(f"  Quality: High")
    print(f"  Frame rate: {FRAME_RATE} fps")
    print("  (This may take several minutes...)")

    run(session, f'movie encode "{OUTPUT_FILE}" format h264 quality high framerate {FRAME_RATE}')

    # Step 10: Summary
    print("\n" + "=" * 70)
    print("MOVIE CREATION COMPLETE!")
    print("=" * 70)

    print(f"\nOutput file: {OUTPUT_FILE}")
    print(f"\nMovie specifications:")
    print(f"  Resolution:     {RESOLUTION[0]}x{RESOLUTION[1]}")
    print(f"  Frame rate:     {FRAME_RATE} fps")
    print(f"  Supersampling:  {SUPERSAMPLE}x")
    print(f"  Surface:        {SURFACE_COLOR} @ {SURFACE_OPACITY*100:.0f}% opacity")
    print(f"  Lighting:       {LIGHTING}")
    print(f"  Silhouettes:    {SILHOUETTES}")
    print(f"  AO shadows:     {AMBIENT_OCCLUSION}")
    print(f"  Background:     {BACKGROUND_COLOR}")

    print("\nReady for publication!")
    print("=" * 70)


# Run the function
create_md_movie(session)
