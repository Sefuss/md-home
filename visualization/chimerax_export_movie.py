# ChimeraX Script: Export High-Quality Movie
# USE THIS AFTER LOADING TRAJECTORY
# Drag into ChimeraX or run: chimerax --script chimerax_export_movie.py

from chimerax.core.commands import run

def export_high_quality_movie(session):
    """
    Export currently loaded trajectory as high-quality MP4 movie
    Optimized for PowerPoint presentations
    """

    output_path = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\04_movies\custom_export.mp4"

    print("=== Exporting High-Quality Movie ===\n")

    # Stop any existing recording first
    print("Preparing recording session...")
    try:
        run(session, "movie stop")
        run(session, "movie reset")
        print("Cleared previous recording")
    except:
        pass

    # Optimize view for recording
    print("Optimizing view...")
    run(session, "view")  # Center molecule
    run(session, "lighting soft")
    run(session, "graphics silhouettes true")

    # Set window size for HD video
    run(session, "windowsize 1920 1080")  # Full HD

    print("Starting recording...")
    run(session, "movie record supersample 3")  # 3x supersampling for quality

    # Play through trajectory
    print("Playing through trajectory...")
    print("(This will take a few minutes depending on frame count)")
    run(session, "coordset #1 1,-1")  # 1,-1 means "from frame 1 to last frame"

    # Wait for playback to complete
    # Note: Adjust wait time based on your trajectory length
    # Rule of thumb: (num_frames / 25 fps) + 5 seconds
    # For 2270 frames: ~95 seconds
    run(session, "wait 100")

    # Encode movie
    print(f"\nEncoding movie to: {output_path}")
    run(session, f"movie encode {output_path} format h264 quality high framerate 25")

    print("\n✅ Movie export complete!")
    print(f"Output file: {output_path}")
    print("\nMovie settings:")
    print("  - Resolution: 1920x1080 (Full HD)")
    print("  - Quality: High (H.264)")
    print("  - Framerate: 25 fps")
    print("  - Supersampling: 3x (anti-aliasing)")
    print("\nReady for PowerPoint!")

# Note: This function is designed to be run interactively
# Uncomment the line below if you want it to run automatically:
# export_high_quality_movie(session)

print("Movie export script loaded!")
print("To export movie, run: export_high_quality_movie(session)")
