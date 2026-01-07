# ChimeraX Script: One-Click AviTag Movie Creation
# Based on TheVisualHub's VisualFactory principles
# Drag this file into ChimeraX to create a professional movie automatically

from chimerax.core.commands import run

def create_avitag_movie(session):
    """
    Create professional AviTag conformational dynamics movie
    Automated, Precise, Beautiful - TheVisualHub philosophy
    """

    print("=" * 70)
    print("CREATING HOLLYWOOD-QUALITY AVITAG MOVIE")
    print("Based on TheVisualHub's VisualFactory Principles")
    print("=" * 70)

    # Load cleaned trajectory (no PBC artifacts, protein only!)
    structure = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\protein_only.gro"
    trajectory = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\protein_only.xtc"

    # Check if files exist
    import os
    if not os.path.exists(structure):
        print(f"ERROR: Structure file not found: {structure}")
        return
    if not os.path.exists(trajectory):
        print(f"ERROR: Trajectory file not found: {trajectory}")
        print("HINT: Run gmx trjconv to create md_centered.xtc first!")
        return

    print("\nStep 1: Loading cleaned trajectory (PBC-corrected)...")
    try:
        run(session, f"open {structure}")
        run(session, f"open {trajectory} structureModel #1")
    except Exception as e:
        print(f"ERROR loading files: {e}")
        return

    # Professional styling
    print("Step 2: Applying professional styling...")
    run(session, "cartoon")
    run(session, "hide solvent")

    # Color by flexibility (B-factor) - shows dynamic regions!
    run(session, "color bfactor #1 palette rainbow")

    # Cinematic lighting (Gleb's signature)
    print("Step 3: Setting up cinematic lighting...")
    run(session, "lighting soft")
    run(session, "graphics silhouettes true")
    run(session, "graphics silhouettes width 2")

    # White background for PowerPoint presentations
    run(session, "set bgColor white")

    # Optimal camera angle for linear peptide
    print("Step 4: Finding optimal camera perspective...")
    run(session, "view")
    run(session, "turn y 30")
    run(session, "turn x 20")
    run(session, "zoom 0.85")

    # Add title
    run(session, "2dlabels create title text 'AviTag: Dynamic Conformational Ensemble' xpos 0.5 ypos 0.95 color black size 20")

    # Set window size for Full HD
    print("Step 5: Setting HD resolution (1920x1080)...")
    run(session, "windowsize 1920 1080")

    # Stop any existing recording (in case script was run before)
    print("\nStep 6: Preparing movie recording...")
    try:
        run(session, "movie stop")
        run(session, "movie reset")
        print("Cleared previous recording session")
    except:
        pass  # No previous recording, that's fine

    # Record with high quality (3x supersampling for smooth edges)
    print("Step 7: Starting high-quality recording...")
    print("(This will take a few minutes - recording 2,396 frames)")
    run(session, "movie record supersample 3")

    # Play through all frames
    print("Step 8: Playing through trajectory...")
    run(session, "coordset #1 1,-1")  # 1,-1 means "from frame 1 to last frame"

    # Wait for playback to complete
    # 2396 frames at 25 fps = ~96 seconds
    print("Step 9: Waiting for playback completion (~100 seconds)...")
    run(session, "wait 100")

    # Encode movie
    output = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\04_movies\lag16_avitag_fusion_dynamics.mp4"
    print(f"\nStep 10: Encoding movie to {output}...")
    run(session, f"movie encode {output} format h264 quality high framerate 25")

    print("\n" + "=" * 70)
    print("MOVIE CREATION COMPLETE!")
    print("=" * 70)
    print(f"\nOutput file: {output}")
    print("\nMovie specifications:")
    print("  Resolution:     1920x1080 (Full HD)")
    print("  Format:         H.264 MP4")
    print("  Quality:        High")
    print("  Framerate:      25 fps")
    print("  Duration:       ~96 seconds")
    print("  Supersampling:  3x (anti-aliasing)")
    print("  Color scheme:   Flexibility (B-factor rainbow)")
    print("\nFeatures:")
    print("  - No PBC artifacts (pre-processed trajectory)")
    print("  - Soft cinematic lighting")
    print("  - Edge silhouettes for depth")
    print("  - White background (PowerPoint-ready)")
    print("  - Professional camera angle")
    print("\nReady to insert into your presentation!")
    print("=" * 70)

# Run the function when script is executed
create_avitag_movie(session)
