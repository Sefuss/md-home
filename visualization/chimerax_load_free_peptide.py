# ChimeraX Script: Load Free AviTag MD Trajectory
# Drag this file into ChimeraX or run: chimerax --script chimerax_load_free_peptide.py

from chimerax.core.commands import run

def load_free_peptide_trajectory(session):
    """
    Load and visualize the free AviTag peptide MD trajectory (23 ns)
    """

    # Paths to your MD files
    structure_file = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\npt.gro"
    trajectory_file = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\md.xtc"

    print("Loading Free AviTag MD trajectory...")
    print(f"Structure: {structure_file}")
    print(f"Trajectory: {trajectory_file}")

    # Load structure
    run(session, f"open {structure_file}")

    # Load trajectory
    run(session, f"open {trajectory_file} structureModel #1")

    # Style the molecule
    run(session, "cartoon")  # Cartoon representation
    run(session, "hide solvent")  # Hide water molecules
    run(session, "color bychain")  # Color by chain

    # Set up nice view
    run(session, "view")  # Center view
    run(session, "lighting soft")  # Soft lighting
    run(session, "graphics silhouettes true")  # Add outline

    # Set background
    run(session, "set bgColor white")  # White background for presentations

    print("\nTrajectory loaded successfully!")
    print(f"Total frames: Use 'coordset #1 1,-1' to play through all frames")
    print("\nUseful commands:")
    print("  - Spacebar: Play/Pause")
    print("  - Right arrow: Next frame")
    print("  - Left arrow: Previous frame")
    print("  - 'movie record' then 'coordset #1 1,-1' then 'movie encode output movie.mp4'")

# Run the function when script is executed
load_free_peptide_trajectory(session)
