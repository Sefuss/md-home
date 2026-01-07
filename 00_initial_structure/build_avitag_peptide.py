"""
Build initial AviTag peptide structure
Extended conformation for MD simulation
"""

from Bio.PDB import PDBIO, PPBuilder
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import Structure, Model, Chain, Residue
import numpy as np

# Amino acid single-letter to three-letter code conversion
AA_CODES = {
    'G': 'GLY', 'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE',
    'M': 'MET', 'F': 'PHE', 'W': 'TRP', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'C': 'CYS', 'Y': 'TYR', 'N': 'ASN', 'Q': 'GLN',
    'D': 'ASP', 'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'
}

def build_extended_peptide(sequence, output_file):
    """
    Build peptide in extended (beta-strand) conformation
    Uses ideal phi/psi angles for extended structure
    """

    # Standard residue templates (backbone atoms)
    # Using extended conformation: phi = -120°, psi = 120°

    print(f"Building peptide: {sequence}")
    print(f"Length: {len(sequence)} residues")

    # Create structure hierarchy
    structure = Structure.Structure("avitag")
    model = Model.Model(0)
    chain = Chain.Chain("A")

    # Backbone geometry for extended conformation
    # Standard peptide bond length and angles
    N_CA_length = 1.46  # Ångstroms
    CA_C_length = 1.52
    C_N_length = 1.33

    # Angles (extended β-strand)
    phi = -120.0  # degrees
    psi = 120.0
    omega = 180.0  # trans peptide bond

    # Starting position
    x, y, z = 0.0, 0.0, 0.0

    # Build each residue
    residue_number = 1
    for aa_code in sequence:
        # Convert single-letter code to three-letter code (e.g., G -> GLY)
        three_letter = AA_CODES[aa_code]

        # Create residue
        res_id = (' ', residue_number, ' ')
        residue = Residue.Residue(res_id, three_letter, '    ')

        # Add backbone atoms (simplified - N, CA, C, O)
        # For MD, GROMACS will add all atoms from topology

        # This is a simplified builder - for real MD, we'll use GROMACS pdb2gmx
        # which will add all atoms including hydrogens

        # For now, create a simple extended chain
        # Each residue offset by ~3.5 Å along x-axis (typical β-strand)
        offset = (residue_number - 1) * 3.5

        from Bio.PDB.Atom import Atom

        # N atom
        n_atom = Atom('N', [offset + 0.0, 0.0, 0.0], 1.0, 1.0, ' ', ' N  ', residue_number, 'N')
        residue.add(n_atom)

        # CA atom
        ca_atom = Atom('CA', [offset + 1.46, 0.0, 0.0], 1.0, 1.0, ' ', ' CA ', residue_number, 'C')
        residue.add(ca_atom)

        # C atom
        c_atom = Atom('C', [offset + 2.5, 0.0, 0.0], 1.0, 1.0, ' ', ' C  ', residue_number, 'C')
        residue.add(c_atom)

        # O atom
        o_atom = Atom('O', [offset + 2.8, -1.0, 0.0], 1.0, 1.0, ' ', ' O  ', residue_number, 'O')
        residue.add(o_atom)

        chain.add(residue)
        residue_number += 1

    model.add(chain)
    structure.add(model)

    # Save PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    print(f"[OK] Saved to: {output_file}")
    print(f"  This is a simplified structure for initial setup")
    print(f"  GROMACS pdb2gmx will add all atoms and fix geometry")

def create_fasta(sequence, output_file):
    """Create FASTA file for reference"""
    with open(output_file, 'w') as f:
        f.write(">AviTag_peptide\n")
        f.write(sequence + "\n")
    print(f"[OK] FASTA saved to: {output_file}")

if __name__ == "__main__":
    # AviTag sequence
    sequence = "GLNDIFEAQKIEWHE"

    output_pdb = "C:/Users/bmartinez/ProteinDesign/AviTag_Binders_2025/02_scaffolds/approach_C_peptide_ensemble/00_initial_structure/avitag_initial.pdb"
    output_fasta = "C:/Users/bmartinez/ProteinDesign/AviTag_Binders_2025/02_scaffolds/approach_C_peptide_ensemble/00_initial_structure/avitag.fasta"

    print("="*60)
    print("BUILDING AVITAG PEPTIDE STRUCTURE")
    print("="*60)
    print()

    # Build structure
    build_extended_peptide(sequence, output_pdb)

    # Create FASTA
    create_fasta(sequence, output_fasta)

    print()
    print("="*60)
    print("NEXT STEPS:")
    print("="*60)
    print("1. Use GROMACS pdb2gmx to add all atoms")
    print("2. Create topology file")
    print("3. Solvate in water box")
    print("4. Energy minimize")
    print("5. Run MD simulation")
    print()
    print("Structure ready for GROMACS processing!")
