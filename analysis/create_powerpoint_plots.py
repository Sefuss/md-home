#!/usr/bin/env python3
"""
Create publication-quality plots for PowerPoint presentation
Shows MD trajectory analysis with clear interpretation
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from datetime import datetime

# Set publication-quality defaults
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

print("="*70)
print("Creating PowerPoint-Ready Visualization Plots")
print("="*70)

# Load trajectory
print("\nLoading trajectory...")
u = mda.Universe('npt.gro', 'md.xtc')
protein = u.select_atoms('protein')
backbone = u.select_atoms('protein and backbone')

n_frames = len(u.trajectory)
time_ns = u.trajectory[-1].time / 1000.0

print(f"  Frames: {n_frames}")
print(f"  Time: {time_ns:.2f} ns")

# Align trajectory
print("\nAligning trajectory...")
aligner = align.AlignTraj(u, u, select='protein and backbone', in_memory=True).run()

# Calculate RMSD
print("Calculating RMSD...")
R = rms.RMSD(u, u, select='protein and backbone', ref_frame=0).run()
rmsd_data = R.results.rmsd
time = rmsd_data[:, 1] / 1000.0  # Convert to ns
rmsd_values = rmsd_data[:, 2]

# Calculate Radius of Gyration
print("Calculating Radius of Gyration...")
rg_values = []
for ts in u.trajectory:
    rg_values.append(protein.radius_of_gyration())
rg_values = np.array(rg_values)

print("\n" + "="*70)
print("Creating Plots...")
print("="*70)

# ============================================================================
# FIGURE 1: Comprehensive 4-panel analysis
# ============================================================================

fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.3)

# Panel A: RMSD over time
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(time, rmsd_values, 'b-', linewidth=2, alpha=0.7)
ax1.axhline(y=np.mean(rmsd_values), color='red', linestyle='--',
            linewidth=2, label=f'Mean: {np.mean(rmsd_values):.1f} Å')
ax1.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
ax1.set_ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
ax1.set_title('A) Backbone RMSD Over Time', fontsize=16, fontweight='bold', loc='left')
ax1.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_xlim(0, time[-1])

# Add interpretation text
if np.std(rmsd_values[int(len(rmsd_values)*0.5):]) < 2.0:
    interpretation = "Plateaued - Sampling stable states"
else:
    interpretation = "Fluctuating - Exploring conformations"
ax1.text(0.95, 0.95, interpretation, transform=ax1.transAxes,
         fontsize=11, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Panel B: Radius of Gyration over time
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(time, rg_values, 'g-', linewidth=2, alpha=0.7)
ax2.axhline(y=np.mean(rg_values), color='red', linestyle='--',
            linewidth=2, label=f'Mean: {np.mean(rg_values):.1f} Å')
ax2.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
ax2.set_ylabel('Radius of Gyration (Å)', fontsize=14, fontweight='bold')
ax2.set_title('B) Radius of Gyration Over Time', fontsize=16, fontweight='bold', loc='left')
ax2.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_xlim(0, time[-1])

# Add state labels
rg_min, rg_max = np.min(rg_values), np.max(rg_values)
ax2.axhspan(rg_min, rg_min + (rg_max-rg_min)*0.3, alpha=0.2, color='blue',
            label='Compact')
ax2.axhspan(rg_max - (rg_max-rg_min)*0.3, rg_max, alpha=0.2, color='orange',
            label='Extended')
ax2.text(time[-1]*0.05, rg_min + (rg_max-rg_min)*0.15, 'Compact',
         fontsize=11, fontweight='bold', color='darkblue')
ax2.text(time[-1]*0.05, rg_max - (rg_max-rg_min)*0.15, 'Extended',
         fontsize=11, fontweight='bold', color='darkorange')

# Panel C: RMSD Distribution
ax3 = fig.add_subplot(gs[1, 0])
counts, bins, patches = ax3.hist(rmsd_values, bins=40, color='skyblue',
                                 edgecolor='black', alpha=0.7)
ax3.axvline(x=np.mean(rmsd_values), color='red', linestyle='--',
            linewidth=2.5, label=f'Mean: {np.mean(rmsd_values):.1f} Å')
ax3.axvline(x=np.median(rmsd_values), color='orange', linestyle='--',
            linewidth=2.5, label=f'Median: {np.median(rmsd_values):.1f} Å')
ax3.set_xlabel('RMSD (Å)', fontsize=14, fontweight='bold')
ax3.set_ylabel('Frequency', fontsize=14, fontweight='bold')
ax3.set_title('C) RMSD Distribution', fontsize=16, fontweight='bold', loc='left')
ax3.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)
ax3.grid(True, alpha=0.3, axis='y', linestyle='--')

# Panel D: Rg Distribution
ax4 = fig.add_subplot(gs[1, 1])
counts, bins, patches = ax4.hist(rg_values, bins=40, color='lightgreen',
                                 edgecolor='black', alpha=0.7)
ax4.axvline(x=np.mean(rg_values), color='red', linestyle='--',
            linewidth=2.5, label=f'Mean: {np.mean(rg_values):.1f} Å')
ax4.axvline(x=np.median(rg_values), color='orange', linestyle='--',
            linewidth=2.5, label=f'Median: {np.median(rg_values):.1f} Å')
ax4.set_xlabel('Radius of Gyration (Å)', fontsize=14, fontweight='bold')
ax4.set_ylabel('Frequency', fontsize=14, fontweight='bold')
ax4.set_title('D) Rg Distribution (Compact vs Extended)', fontsize=16, fontweight='bold', loc='left')
ax4.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)
ax4.grid(True, alpha=0.3, axis='y', linestyle='--')

# Overall title
fig.suptitle('MD Trajectory Analysis: AviTag Conformational Dynamics',
             fontsize=18, fontweight='bold', y=0.995)

# Add timestamp and statistics
stats_text = f"Time: {time_ns:.1f} ns | Frames: {n_frames} | RMSD: {np.mean(rmsd_values):.1f}±{np.std(rmsd_values):.1f} Å | Rg: {np.mean(rg_values):.1f}±{np.std(rg_values):.1f} Å"
fig.text(0.5, 0.01, stats_text, ha='center', fontsize=10, style='italic')

plt.savefig('MD_Analysis_4Panel.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  ✓ Saved: MD_Analysis_4Panel.png (4-panel comprehensive)")

plt.close()

# ============================================================================
# FIGURE 2: Time series with state annotations (PowerPoint focus)
# ============================================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Top: RMSD with state coloring
scatter = ax1.scatter(time, rmsd_values, c=rg_values, cmap='RdYlGn_r',
                     s=30, alpha=0.6, edgecolors='black', linewidth=0.5)
ax1.plot(time, rmsd_values, 'gray', linewidth=0.5, alpha=0.3, zorder=0)
cbar = plt.colorbar(scatter, ax=ax1, label='Rg (Å)')
cbar.set_label('Radius of Gyration (Å)', fontsize=12, fontweight='bold')
ax1.axhline(y=np.mean(rmsd_values), color='red', linestyle='--', linewidth=2,
            label=f'Mean RMSD: {np.mean(rmsd_values):.1f} Å')
ax1.set_ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
ax1.set_title('Conformational Exploration: RMSD colored by Compactness',
              fontsize=16, fontweight='bold')
ax1.legend(fontsize=12, loc='upper right')
ax1.grid(True, alpha=0.3, linestyle='--')

# Bottom: Rg with state regions
ax2.plot(time, rg_values, 'b-', linewidth=2.5, alpha=0.8)
ax2.fill_between(time, rg_values, alpha=0.3, color='skyblue')

# Mark compact and extended regions
compact_threshold = np.percentile(rg_values, 25)
extended_threshold = np.percentile(rg_values, 75)

ax2.axhline(y=compact_threshold, color='blue', linestyle=':', linewidth=2,
            label=f'Compact (<{compact_threshold:.1f} Å)')
ax2.axhline(y=extended_threshold, color='orange', linestyle=':', linewidth=2,
            label=f'Extended (>{extended_threshold:.1f} Å)')
ax2.axhline(y=np.mean(rg_values), color='red', linestyle='--', linewidth=2,
            label=f'Mean: {np.mean(rg_values):.1f} Å')

# Shade regions
ax2.axhspan(np.min(rg_values), compact_threshold, alpha=0.15, color='blue')
ax2.axhspan(extended_threshold, np.max(rg_values), alpha=0.15, color='orange')

ax2.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
ax2.set_ylabel('Radius of Gyration (Å)', fontsize=14, fontweight='bold')
ax2.set_title('Peptide Size: Compact ↔ Extended State Transitions',
              fontsize=16, fontweight='bold')
ax2.legend(fontsize=11, loc='upper right')
ax2.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('MD_TimeSeries_States.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  ✓ Saved: MD_TimeSeries_States.png (state transitions)")

plt.close()

# ============================================================================
# FIGURE 3: State space exploration (2D scatter)
# ============================================================================

fig, ax = plt.subplots(figsize=(10, 8))

# Create 2D density plot
from scipy.stats import gaussian_kde

# Sample data for KDE (use every 5th point for speed)
sample_idx = np.arange(0, len(rmsd_values), 5)
xy = np.vstack([rmsd_values[sample_idx], rg_values[sample_idx]])
z = gaussian_kde(xy)(xy)

scatter = ax.scatter(rmsd_values, rg_values, c=time, cmap='viridis',
                    s=50, alpha=0.6, edgecolors='black', linewidth=0.5)
cbar = plt.colorbar(scatter, ax=ax, label='Time (ns)')
cbar.set_label('Time (ns)', fontsize=12, fontweight='bold')

# Mark start and end
ax.scatter(rmsd_values[0], rg_values[0], s=300, c='green', marker='*',
          edgecolors='black', linewidth=2, label='Start', zorder=10)
ax.scatter(rmsd_values[-1], rg_values[-1], s=300, c='red', marker='*',
          edgecolors='black', linewidth=2, label='End', zorder=10)

ax.set_xlabel('RMSD (Å)', fontsize=14, fontweight='bold')
ax.set_ylabel('Radius of Gyration (Å)', fontsize=14, fontweight='bold')
ax.set_title('Conformational State Space Exploration', fontsize=16, fontweight='bold')
ax.legend(fontsize=12, loc='best', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle='--')

# Add quadrant labels
ax.text(0.95, 0.95, 'Extended\nHigh RMSD', transform=ax.transAxes,
        fontsize=11, ha='right', va='top',
        bbox=dict(boxstyle='round', facecolor='orange', alpha=0.3))
ax.text(0.05, 0.05, 'Compact\nLow RMSD', transform=ax.transAxes,
        fontsize=11, ha='left', va='bottom',
        bbox=dict(boxstyle='round', facecolor='blue', alpha=0.3))

plt.tight_layout()
plt.savefig('MD_StateSpace_2D.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  ✓ Saved: MD_StateSpace_2D.png (conformational space)")

plt.close()

# ============================================================================
# FIGURE 4: Simple summary for PowerPoint slide
# ============================================================================

fig, ax = plt.subplots(figsize=(12, 6))

# Plot both metrics on same axes with different scales
ax1 = ax
ax2 = ax.twinx()

# RMSD on left axis
line1 = ax1.plot(time, rmsd_values, 'b-', linewidth=3, alpha=0.7, label='RMSD')
ax1.set_xlabel('Time (ns)', fontsize=16, fontweight='bold')
ax1.set_ylabel('RMSD (Å)', fontsize=16, fontweight='bold', color='blue')
ax1.tick_params(axis='y', labelcolor='blue', labelsize=12)
ax1.grid(True, alpha=0.3, linestyle='--')

# Rg on right axis
line2 = ax2.plot(time, rg_values, 'g-', linewidth=3, alpha=0.7, label='Radius of Gyration')
ax2.set_ylabel('Radius of Gyration (Å)', fontsize=16, fontweight='bold', color='green')
ax2.tick_params(axis='y', labelcolor='green', labelsize=12)

# Combined legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, fontsize=14, loc='upper left', frameon=True,
          fancybox=True, shadow=True)

ax1.set_title('AviTag Peptide: Dynamic Conformational Sampling',
             fontsize=18, fontweight='bold', pad=20)

# Add interpretation box
interpretation_text = f"""
MD Simulation Results ({time_ns:.1f} ns):
• RMSD fluctuates: {np.min(rmsd_values):.1f}-{np.max(rmsd_values):.1f} Å
• Rg fluctuates: {np.min(rg_values):.1f}-{np.max(rg_values):.1f} Å
• Conclusion: Peptide explores multiple states
"""

ax1.text(0.98, 0.97, interpretation_text, transform=ax1.transAxes,
        fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=1),
        family='monospace')

plt.tight_layout()
plt.savefig('MD_Summary_Simple.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  ✓ Saved: MD_Summary_Simple.png (PowerPoint summary)")

plt.close()

# ============================================================================
# Analysis Summary
# ============================================================================

print("\n" + "="*70)
print("ANALYSIS COMPLETE - Check for State Transitions")
print("="*70)

# Check if Rg oscillates (visits both states multiple times)
compact_visits = np.sum(rg_values < compact_threshold)
extended_visits = np.sum(rg_values > extended_threshold)
middle_visits = np.sum((rg_values >= compact_threshold) & (rg_values <= extended_threshold))

# Count transitions
transitions = 0
current_state = 'middle'
if rg_values[0] < compact_threshold:
    current_state = 'compact'
elif rg_values[0] > extended_threshold:
    current_state = 'extended'

for rg in rg_values[1:]:
    if rg < compact_threshold and current_state != 'compact':
        transitions += 1
        current_state = 'compact'
    elif rg > extended_threshold and current_state != 'extended':
        transitions += 1
        current_state = 'extended'
    elif compact_threshold <= rg <= extended_threshold:
        current_state = 'middle'

print(f"\nState Occupancy:")
print(f"  Compact (<{compact_threshold:.1f} Å): {compact_visits} frames ({compact_visits/len(rg_values)*100:.1f}%)")
print(f"  Middle: {middle_visits} frames ({middle_visits/len(rg_values)*100:.1f}%)")
print(f"  Extended (>{extended_threshold:.1f} Å): {extended_visits} frames ({extended_visits/len(rg_values)*100:.1f}%)")
print(f"\nState Transitions: {transitions}")

if transitions > 5:
    print("\n✓ DYNAMIC EXPLORATION: Peptide visits both compact and extended states multiple times!")
    print("  This is EXCELLENT for ensemble-based binder design.")
elif transitions > 0:
    print("\n⚠ MODERATE EXPLORATION: Some transitions between states")
else:
    print("\n⚠ PROGRESSIVE CHANGE: Peptide may be unfolding in one direction")

# Final statistics
print(f"\nFinal Statistics:")
print(f"  RMSD: {np.mean(rmsd_values):.2f} ± {np.std(rmsd_values):.2f} Å")
print(f"  Rg:   {np.mean(rg_values):.2f} ± {np.std(rg_values):.2f} Å")
print(f"  Range (Rg): {np.min(rg_values):.1f} - {np.max(rg_values):.1f} Å")

print("\n" + "="*70)
print("PowerPoint-Ready Plots Created!")
print("="*70)
print("\nFiles created:")
print("  1. MD_Analysis_4Panel.png          (Comprehensive 4-panel figure)")
print("  2. MD_TimeSeries_States.png        (Time series with state regions)")
print("  3. MD_StateSpace_2D.png             (2D conformational space)")
print("  4. MD_Summary_Simple.png            (Simple overview for slides)")
print("\nAll plots are 300 DPI, suitable for publication/presentation.")
print("="*70)
