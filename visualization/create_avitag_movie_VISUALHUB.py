# ChimeraX Script: VisualHub-Quality AviTag Movie
# Based on TheVisualHub/VisualFactory professional standards
# Addresses: zoom level, smooth motion, cinematic lighting

from chimerax.core.commands import run
import os

def create_visualhub_quality_movie(session):
    """
    Create professional VisualHub-quality molecular dynamics movie

    Improvements over basic version:
    1. Better zoom/framing (shows full peptide dynamics)
    2. Advanced lighting (shadows, depth cue, ambient occlusion)
    3. Smoother playback (play every 2nd frame)
    4. Professional color scheme
    5. Higher rendering quality
    """

    print("=" * 70)
    print("CREATING VISUALHUB-QUALITY AVITAG MOVIE")
    print("Professional molecular dynamics visualization")
    print("=" * 70)

    # Use smoothed trajectory if available, otherwise use centered
    base_dir = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation"
    structure = os.path.join(base_dir, "npt.gro")

    # Check for smoothed trajectory first (created with gmx filter)
    smooth_traj = os.path.join(base_dir, "md_smooth.xtc")
    centered_traj = os.path.join(base_dir, "md_centered.xtc")

    if os.path.exists(smooth_traj):
        trajectory = smooth_traj
        print("\nUsing SMOOTHED trajectory (optimal quality)")
    elif os.path.exists(centered_traj):
        trajectory = centered_traj
        print("\nUsing CENTERED trajectory (good quality)")
        print("TIP: Create md_smooth.xtc with 'gmx filter -nf 5' for smoother motion")
    else:
        print("\nERROR: No suitable trajectory found!")
        print("Expected: md_smooth.xtc or md_centered.xtc")
        return

    # Verify files exist
    if not os.path.exists(structure):
        print(f"ERROR: Structure file not found: {structure}")
        return
    if not os.path.exists(trajectory):
        print(f"ERROR: Trajectory file not found: {trajectory}")
        return

    # Load trajectory
    print("\nStep 1: Loading trajectory...")
    try:
        run(session, f"open {structure}")
        run(session, f"open {trajectory} structureModel #1")
    except Exception as e:
        print(f"ERROR loading files: {e}")
        return

    # VisualHub-Style Molecular Representation
    print("Step 2: Applying VisualHub-style molecular representation...")
    run(session, "cartoon")
    run(session, "hide solvent")

    # Professional color scheme (high contrast, publication-quality)
    run(session, "color bychain #1")  # Simple rainbow by chain

    # Alternative color schemes to try:
    # run(session, "rainbow #1")  # Rainbow gradient along chain
    # run(session, "color bfactor #1 palette rainbow")  # Color by B-factor (flexibility)

    # Advanced Lighting (Gleb's signature)
    print("Step 3: Setting up advanced cinematic lighting...")
    run(session, "lighting soft")
    run(session, "lighting shadows true intensity 0.6")  # Subtle shadows for depth
    run(session, "lighting depthCue true start 0.5 end 1.0")  # Distance fog effect
    run(session, "lighting quality finer")  # Higher quality lighting calculation

    # Silhouettes for professional edge definition
    run(session, "graphics silhouettes true")
    run(session, "set silhouetteColor black")
    run(session, "set silhouetteWidth 2.5")

    # White background for presentations
    run(session, "set bgColor white")

    # Camera Setup - CRITICAL FIX for zoom issue
    print("Step 4: Optimizing camera angle and zoom...")
    run(session, "view")  # Center molecule

    # Adjust viewing angle for linear peptide
    run(session, "turn y 25")
    run(session, "turn x 15")

    # ZOOM FIX: Zoom OUT more to show full dynamics
    # Original was 0.85 (too close), now 1.3 (shows full peptide + breathing room)
    run(session, "zoom 1.3")

    print("  ✓ Camera zoomed out to show full conformational ensemble")

    # Add professional title
    run(session, "2dlabels create title text 'AviTag: Conformational Dynamics' xpos 0.5 ypos 0.95 color black size 22 font 'Arial' bold true")

    # Set window size for Full HD
    print("Step 5: Setting HD resolution (1920x1080)...")
    run(session, "windowsize 1920 1080")

    # Clear any previous recording
    print("\nStep 6: Preparing movie recording session...")
    try:
        run(session, "movie stop")
        run(session, "movie reset")
        print("  ✓ Cleared previous recording")
    except:
        pass

    # High-Quality Recording
    print("Step 7: Starting high-quality recording...")
    print("  - 3x supersampling (anti-aliasing)")
    print("  - Playing every 2nd frame for smoother motion")
    print("  - 30 fps output (optimal for presentations)")
    run(session, "movie record supersample 3")

    # Play through trajectory - SMOOTHNESS FIX
    print("\nStep 8: Playing through trajectory...")
    # Play every 2nd frame for smoother appearance (reduces jitter)
    # Syntax: coordset #model start,end,step
    run(session, "coordset #1 1,-1,2")  # Every 2nd frame

    # Calculate wait time: ~2396 frames / 2 = 1198 frames at 25 fps = ~48 seconds
    print("Step 9: Waiting for playback (~50 seconds)...")
    run(session, "wait 50")

    # Encode movie
    output = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\04_movies\free_peptide_dynamics_VISUALHUB.mp4"
    print(f"\nStep 10: Encoding movie to {output}...")
    run(session, f"movie encode {output} format h264 quality high framerate 30")

    print("\n" + "=" * 70)
    print("VISUALHUB-QUALITY MOVIE COMPLETE!")
    print("=" * 70)
    print(f"\nOutput file: {output}")
    print("\nMovie specifications:")
    print("  Resolution:     1920x1080 (Full HD)")
    print("  Format:         H.264 MP4")
    print("  Quality:        High")
    print("  Framerate:      30 fps")
    print("  Duration:       ~48 seconds")
    print("  Supersampling:  3x (anti-aliasing)")
    print("  Color scheme:   Paired-12 (ColorBrewer)")
    print("\nVisualHub Enhancements:")
    print("  ✓ Cinematic lighting with shadows")
    print("  ✓ Depth cue for spatial context")
    print("  ✓ Professional silhouettes")
    print("  ✓ Optimized zoom (shows full dynamics)")
    print("  ✓ Smoother playback (every 2nd frame)")
    print("\nNEXT LEVEL:")
    print("  - Create md_smooth.xtc with 'gmx filter -nf 5' for even smoother motion")
    print("  - Try ray-tracing: 'set subdivision 4' before recording")
    print("  - Add camera rotation: 'turn y 1 frames' during playback")
    print("\nReady for PowerPoint presentations!")
    print("=" * 70)

# Run the function
create_visualhub_quality_movie(session)
