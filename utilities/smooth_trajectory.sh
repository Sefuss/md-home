#!/bin/bash
# Smooth MD Trajectory for High-Quality Visualization
# Uses GROMACS filter to apply temporal averaging

# Usage:
#   bash smooth_trajectory.sh <input.xtc> <structure.gro> <output.xtc> [nframes]
#
# Example:
#   bash smooth_trajectory.sh md_centered.xtc npt.gro md_smooth.xtc 5

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "ERROR: Missing arguments"
    echo ""
    echo "Usage: bash smooth_trajectory.sh <input.xtc> <structure.gro> <output.xtc> [nframes]"
    echo ""
    echo "Arguments:"
    echo "  input.xtc    - Input trajectory (e.g., md_centered.xtc)"
    echo "  structure.gro - Structure file (e.g., npt.gro)"
    echo "  output.xtc   - Output smoothed trajectory (e.g., md_smooth.xtc)"
    echo "  nframes      - Number of frames to average (default: 5)"
    echo ""
    echo "Larger nframes = smoother but more blurred motion"
    echo "  nframes=3  - Subtle smoothing (good for fast dynamics)"
    echo "  nframes=5  - Moderate smoothing (recommended for presentation)"
    echo "  nframes=10 - Heavy smoothing (good for very jittery trajectories)"
    exit 1
fi

INPUT_XTC=$1
STRUCTURE=$2
OUTPUT_XTC=$3
NFRAMES=${4:-5}  # Default to 5 frames if not specified

echo "========================================"
echo "Smoothing MD Trajectory for Visualization"
echo "========================================"
echo ""
echo "Input trajectory:  $INPUT_XTC"
echo "Structure file:    $STRUCTURE"
echo "Output trajectory: $OUTPUT_XTC"
echo "Averaging window:  $NFRAMES frames"
echo ""

# Check if files exist
if [ ! -f "$INPUT_XTC" ]; then
    echo "ERROR: Input trajectory not found: $INPUT_XTC"
    exit 1
fi

if [ ! -f "$STRUCTURE" ]; then
    echo "ERROR: Structure file not found: $STRUCTURE"
    exit 1
fi

# Get absolute path to current directory
WORKDIR=$(pwd)

# Convert to Windows path for Docker (if needed)
if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    # We're on Windows Git Bash
    WORKDIR_WIN=$(cygpath -w "$WORKDIR")
    echo "Running on Windows: $WORKDIR_WIN"

    # Run GROMACS filter in Docker
    export MSYS_NO_PATHCONV=1
    docker run --rm -v "$WORKDIR_WIN:/workspace" -w /workspace gromacs/gromacs:latest \
        gmx filter -f "$INPUT_XTC" -s "$STRUCTURE" -nf "$NFRAMES" -o "$OUTPUT_XTC" <<EOF
0
EOF

else
    # Linux/Mac
    docker run --rm -v "$WORKDIR:/workspace" -w /workspace gromacs/gromacs:latest \
        gmx filter -f "$INPUT_XTC" -s "$STRUCTURE" -nf "$NFRAMES" -o "$OUTPUT_XTC" <<EOF
0
EOF
fi

echo ""
echo "========================================"
echo "Trajectory smoothing complete!"
echo "========================================"
echo ""
echo "Output file: $OUTPUT_XTC"
echo ""
echo "Next steps:"
echo "  1. Load smoothed trajectory in ChimeraX"
echo "  2. Run create_avitag_movie_VISUALHUB.py for best quality"
echo ""
