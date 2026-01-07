#!/usr/bin/env python3
"""
Quick MD Trajectory Analysis
Analyzes stability, artifacts, and quality of ongoing MD simulation
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def analyze_trajectory(gro_file, xtc_file, output_prefix="analysis"):
    """
    Comprehensive trajectory analysis
    """
    print("="*60)
    print("MD Trajectory Analysis")
    print("="*60)

    # Load trajectory
    print(f"\nLoading trajectory...")
    print(f"  Structure: {gro_file}")
    print(f"  Trajectory: {xtc_file}")

    u = mda.Universe(gro_file, xtc_file)
    protein = u.select_atoms('protein')

    n_frames = len(u.trajectory)
    time_ns = u.trajectory[-1].time / 1000.0  # Convert ps to ns

    print(f"  Frames: {n_frames}")
    print(f"  Time: {time_ns:.2f} ns")
    print(f"  Peptide atoms: {len(protein)}")

    # 1. RMSD Analysis
    print("\n" + "="*60)
    print("1. RMSD Analysis (Backbone Stability)")
    print("="*60)

    backbone = u.select_atoms('protein and backbone')

    # Align trajectory to first frame
    aligner = align.AlignTraj(u, u, select='protein and backbone',
                              in_memory=True).run()

    # Calculate RMSD
    R = rms.RMSD(u, u, select='protein and backbone',
                 ref_frame=0).run()

    rmsd_data = R.results.rmsd
    time = rmsd_data[:, 1] / 1000.0  # Convert to ns
    rmsd_values = rmsd_data[:, 2]  # RMSD in Angstroms

    print(f"  Initial RMSD: {rmsd_values[0]:.3f} Å (should be ~0)")
    print(f"  Final RMSD: {rmsd_values[-1]:.3f} Å")
    print(f"  Average RMSD: {np.mean(rmsd_values):.3f} Å")
    print(f"  Max RMSD: {np.max(rmsd_values):.3f} Å")
    print(f"  Std Dev: {np.std(rmsd_values):.3f} Å")

    # Interpretation
    avg_rmsd = np.mean(rmsd_values)
    if avg_rmsd < 2.0:
        print("  ✅ STABLE - Peptide maintaining structure")
    elif avg_rmsd < 4.0:
        print("  ⚠️  MODERATE - Some conformational changes")
    else:
        print("  ❌ UNSTABLE - Significant structural drift")

    # 2. Radius of Gyration
    print("\n" + "="*60)
    print("2. Radius of Gyration (Compactness)")
    print("="*60)

    rg_values = []
    for ts in u.trajectory:
        rg_values.append(protein.radius_of_gyration())

    rg_values = np.array(rg_values)

    print(f"  Initial Rg: {rg_values[0]:.3f} Å")
    print(f"  Final Rg: {rg_values[-1]:.3f} Å")
    print(f"  Average Rg: {np.mean(rg_values):.3f} Å")
    print(f"  Change: {rg_values[-1] - rg_values[0]:.3f} Å")
    print(f"  Std Dev: {np.std(rg_values):.3f} Å")

    # Interpretation
    rg_change = abs(rg_values[-1] - rg_values[0])
    if rg_change < 1.0:
        print("  ✅ COMPACT - Peptide not unfolding")
    elif rg_change < 2.0:
        print("  ⚠️  SLIGHT EXPANSION - Minor structural changes")
    else:
        print("  ❌ EXPANDING - Possible unfolding")

    # 3. Per-residue RMSF (flexibility)
    print("\n" + "="*60)
    print("3. Per-Residue Flexibility (RMSF)")
    print("="*60)

    # Reset trajectory
    u.trajectory[0]
    aligner = align.AlignTraj(u, u, select='protein and backbone',
                              in_memory=True).run()

    # Calculate average positions
    avg_coords = np.zeros((len(protein), 3))
    for ts in u.trajectory:
        avg_coords += protein.positions
    avg_coords /= len(u.trajectory)

    # Calculate RMSF
    rmsf_values = np.zeros(len(protein))
    for ts in u.trajectory:
        rmsf_values += np.sum((protein.positions - avg_coords)**2, axis=1)
    rmsf_values = np.sqrt(rmsf_values / len(u.trajectory))

    # Average RMSF per residue
    residues = protein.residues
    rmsf_per_residue = []
    for residue in residues:
        atoms = residue.atoms
        res_rmsf = np.mean([rmsf_values[i] for i in range(len(protein))
                           if protein[i] in atoms])
        rmsf_per_residue.append(res_rmsf)

    print(f"  Average RMSF: {np.mean(rmsf_per_residue):.3f} Å")
    print(f"  Most flexible residue: {np.max(rmsf_per_residue):.3f} Å")
    print(f"  Least flexible residue: {np.min(rmsf_per_residue):.3f} Å")

    if np.mean(rmsf_per_residue) < 2.0:
        print("  ✅ LOW FLEXIBILITY - Stable structure")
    else:
        print("  ⚠️  HIGH FLEXIBILITY - Dynamic/disordered")

    # 4. Create Plots
    print("\n" + "="*60)
    print("4. Generating Plots")
    print("="*60)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: RMSD over time
    ax1 = axes[0, 0]
    ax1.plot(time, rmsd_values, 'b-', linewidth=1.5)
    ax1.set_xlabel('Time (ns)', fontsize=11)
    ax1.set_ylabel('RMSD (Å)', fontsize=11)
    ax1.set_title('Backbone RMSD vs Time', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=np.mean(rmsd_values), color='r', linestyle='--',
                label=f'Mean: {np.mean(rmsd_values):.2f} Å')
    ax1.legend()

    # Plot 2: Radius of gyration
    ax2 = axes[0, 1]
    ax2.plot(time, rg_values, 'g-', linewidth=1.5)
    ax2.set_xlabel('Time (ns)', fontsize=11)
    ax2.set_ylabel('Rg (Å)', fontsize=11)
    ax2.set_title('Radius of Gyration vs Time', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=np.mean(rg_values), color='r', linestyle='--',
                label=f'Mean: {np.mean(rg_values):.2f} Å')
    ax2.legend()

    # Plot 3: RMSD histogram
    ax3 = axes[1, 0]
    ax3.hist(rmsd_values, bins=50, color='cyan', edgecolor='black', alpha=0.7)
    ax3.set_xlabel('RMSD (Å)', fontsize=11)
    ax3.set_ylabel('Frequency', fontsize=11)
    ax3.set_title('RMSD Distribution', fontsize=12, fontweight='bold')
    ax3.axvline(x=np.mean(rmsd_values), color='r', linestyle='--',
                label=f'Mean: {np.mean(rmsd_values):.2f} Å')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')

    # Plot 4: Per-residue RMSF
    ax4 = axes[1, 1]
    residue_ids = [r.resid for r in residues]
    ax4.bar(residue_ids, rmsf_per_residue, color='orange', alpha=0.7)
    ax4.set_xlabel('Residue Number', fontsize=11)
    ax4.set_ylabel('RMSF (Å)', fontsize=11)
    ax4.set_title('Per-Residue Flexibility', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    plt.suptitle(f'MD Trajectory Analysis ({time_ns:.1f} ns)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()

    plot_file = f"{output_prefix}_stability.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot_file}")

    # 5. Summary Report
    print("\n" + "="*60)
    print("SUMMARY: Quality Assessment")
    print("="*60)

    issues = []

    # Check RMSD plateau
    rmsd_last_quarter = rmsd_values[int(len(rmsd_values)*0.75):]
    rmsd_trend = np.polyfit(range(len(rmsd_last_quarter)), rmsd_last_quarter, 1)[0]

    if abs(rmsd_trend) < 0.01:
        print("✅ RMSD plateaued - System equilibrated")
    else:
        print(f"⚠️  RMSD still changing ({rmsd_trend:+.4f} Å/frame)")
        issues.append("RMSD not fully plateaued")

    # Check for drift
    if np.mean(rmsd_values) < 3.0:
        print("✅ Low RMSD - Structure stable")
    else:
        print("⚠️  High RMSD - Possible structural drift")
        issues.append("High RMSD values")

    # Check compactness
    if np.std(rg_values) < 1.0:
        print("✅ Rg stable - Compact structure maintained")
    else:
        print("⚠️  Rg fluctuating - Possible conformational changes")
        issues.append("Rg fluctuations")

    # Check flexibility
    if np.max(rmsf_per_residue) < 5.0:
        print("✅ Reasonable flexibility - No extreme motions")
    else:
        print("⚠️  High flexibility detected")
        issues.append("High RMSF values")

    print("\n" + "="*60)
    if len(issues) == 0:
        print("VERDICT: ✅ Simulation looks GOOD - Continue running!")
    else:
        print("VERDICT: ⚠️  Some concerns detected:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nRecommendation: Review plots and decide if acceptable")
    print("="*60)

    # Return statistics
    return {
        'time_ns': time_ns,
        'n_frames': n_frames,
        'rmsd_mean': np.mean(rmsd_values),
        'rmsd_max': np.max(rmsd_values),
        'rg_mean': np.mean(rg_values),
        'rg_change': rg_change,
        'rmsf_mean': np.mean(rmsf_per_residue),
        'issues': issues
    }


if __name__ == "__main__":
    # File paths
    gro_file = "npt.gro"
    xtc_file = "md.xtc"

    if not os.path.exists(gro_file):
        print(f"ERROR: Structure file not found: {gro_file}")
        sys.exit(1)

    if not os.path.exists(xtc_file):
        print(f"ERROR: Trajectory file not found: {xtc_file}")
        sys.exit(1)

    # Run analysis
    stats = analyze_trajectory(gro_file, xtc_file)

    print("\nAnalysis complete!")
