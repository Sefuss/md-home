#!/bin/bash
# Quick energy analysis using GROMACS tools

echo "======================================================================"
echo "Energy Analysis - Temperature, Pressure, Total Energy"
echo "======================================================================"

# Extract temperature
echo ""
echo "Extracting temperature data..."
echo -e "Temperature\n0" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx energy -f md.edr -o temperature.xvg

# Extract pressure
echo ""
echo "Extracting pressure data..."
echo -e "Pressure\n0" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx energy -f md.edr -o pressure.xvg

# Extract total energy
echo ""
echo "Extracting total energy..."
echo -e "Total-Energy\n0" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx energy -f md.edr -o total_energy.xvg

# Extract potential energy
echo ""
echo "Extracting potential energy..."
echo -e "Potential\n0" | docker run --rm -i \
  -v "$(pwd):/work" -w /work \
  gromacs/gromacs:latest \
  gmx energy -f md.edr -o potential.xvg

echo ""
echo "======================================================================"
echo "Energy files created:"
echo "  - temperature.xvg"
echo "  - pressure.xvg"
echo "  - total_energy.xvg"
echo "  - potential.xvg"
echo "======================================================================"
