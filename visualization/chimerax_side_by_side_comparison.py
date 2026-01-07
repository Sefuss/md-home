# ChimeraX Script: Side-by-Side Comparison (Free vs Fusion AviTag)
# USE THIS AFTER FUSION MD COMPLETES
# Drag this file into ChimeraX or run: chimerax --script chimerax_side_by_side_comparison.py

from chimerax.core.commands import run

def create_side_by_side_comparison(session):
    """
    Load both free and fusion AviTag trajectories for side-by-side comparison
    """

    # Paths to free peptide MD
    free_structure = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\npt.gro"
    free_trajectory = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\01_md_preparation\md.xtc"

    # Paths to fusion protein MD
    fusion_structure = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\03_fusion_md\npt.gro"
    fusion_trajectory = r"C:\Users\bmartinez\ProteinDesign\AviTag_Binders_2025\02_scaffolds\approach_C_peptide_ensemble\03_fusion_md\md.xtc"

    print("=== Side-by-Side Comparison: Free vs Fusion AviTag ===\n")

    # Load free peptide (will be model #1)
    print("Loading FREE peptide trajectory...")
    run(session, f"open {free_structure}")
    run(session, f"open {free_trajectory} structureModel #1")
    run(session, "name #1 'Free AviTag'")

    # Load fusion protein (will be model #2)
    print("Loading FUSION protein trajectory...")
    run(session, f"open {fusion_structure}")
    run(session, f"open {fusion_trajectory} structureModel #2")
    run(session, "name #2 'LaG16-AviTag Fusion'")

    # Style both identically
    print("\nStyling molecules...")
    run(session, "cartoon #1,2")
    run(session, "hide #1,2 solvent")
    run(session, "color bychain #1,2")

    # Split screen view
    print("Setting up split screen...")
    run(session, "split h")  # Horizontal split

    # Position models in separate screens
    run(session, "view #1 screen left")
    run(session, "view #2 screen right")

    # Add labels
    run(session, "2dlabels text 'Free AviTag' xpos 0.1 ypos 0.95 color black size 24")
    run(session, "2dlabels text 'Fusion Protein' xpos 0.6 ypos 0.95 color black size 24")

    # Lighting and effects
    run(session, "lighting soft")
    run(session, "graphics silhouettes true")
    run(session, "set bgColor white")

    print("\n✅ Side-by-side comparison ready!")
    print("\nTo play synchronized:")
    print("  coordset #1,2 1,-1")
    print("\nTo record movie:")
    print("  movie stop; movie reset  # Clear any previous recording")
    print("  movie record supersample 3")
    print("  coordset #1,2 1,-1")
    print("  wait 100")
    print("  movie encode C:\\Users\\bmartinez\\ProteinDesign\\AviTag_Binders_2025\\02_scaffolds\\approach_C_peptide_ensemble\\04_movies\\free_vs_fusion_comparison.mp4 format h264 quality high framerate 25")

# Run the function
create_side_by_side_comparison(session)
