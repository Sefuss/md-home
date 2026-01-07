# ChimeraX Script: Load CLEANED Free AviTag MD Trajectory (No PBC Wrapping!)
# Drag this file into ChimeraX for presentation-ready visualization

from chimerax.core.commands import run

def load_clean_trajectory(session):
    """
    Load the cleaned, centered trajectory without PBC artifacts
    Perfect for making presentation movies
    """

    # Paths to CLEANED MD files
    structure_file = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\npt.gro"
    trajectory_file = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\md_centered.xtc"

    print("=" * 60)
    print("Loading CLEANED Free AviTag MD Trajectory")
    print("(No periodic boundary artifacts - presentation ready!)")
    print("=" * 60)

    # Load structure
    run(session, f"open {structure_file}")

    # Load CLEANED trajectory
    run(session, f"open {trajectory_file} structureModel #1")

    # Professional styling for presentations
    print("\nApplying professional styling...")

    # Cartoon representation
    run(session, "cartoon")

    # Hide water molecules
    run(session, "hide solvent")

    # Color scheme options - pick what looks best
    # Option 1: Color by secondary structure (default)
    run(session, "color byhet")

    # Option 2: Color by chain (uncomment if preferred)
    # run(session, "color bychain")

    # Option 3: Single professional color (uncomment if preferred)
    # run(session, "color #1 cornflowerblue")

    # Professional lighting and effects
    run(session, "lighting soft")
    run(session, "graphics silhouettes true")  # Adds nice outline
    run(session, "graphics silhouettes width 2")  # Thicker outline

    # Background - white for PowerPoint
    run(session, "set bgColor white")

    # Set up camera
    run(session, "view")  # Center view
    run(session, "turn y 20")  # Slight rotation for depth

    print("\n" + "=" * 60)
    print("Trajectory loaded successfully!")
    print("=" * 60)
    print(f"\nTotal frames: 2,396 frames (23.95 ns)")
    print("\nPlayback Controls:")
    print("  Spacebar:     Play/Pause")
    print("  Right arrow:  Next frame")
    print("  Left arrow:   Previous frame")
    print("  Home:         First frame")
    print("  End:          Last frame")
    print("\nCamera Controls:")
    print("  Mouse drag:   Rotate")
    print("  Scroll:       Zoom")
    print("  Shift+drag:   Translate")
    print("\nQuick Commands:")
    print("  'view'        - Re-center molecule")
    print("  'turn y 90'   - Rotate 90 degrees around Y")
    print("  'movie record' - Start recording")
    print("  'coordset #1 1,-1' - Play all frames")
    print("=" * 60)

# Run the function when script is executed
load_clean_trajectory(session)
