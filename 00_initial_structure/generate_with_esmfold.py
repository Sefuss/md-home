"""
Generate AviTag peptide structure using ESMFold
This will create a proper all-atom structure that GROMACS can use
"""

import requests
import time

def generate_structure_esmfold(sequence, output_file):
    """
    Generate structure using ESMFold API
    """
    print(f"Generating structure for: {sequence}")
    print(f"Length: {len(sequence)} residues")
    print("Using ESMFold via Hugging Face API...")

    # ESMFold API endpoint (updated router)
    API_URL = "https://router.huggingface.co/models/facebook/esmfold_v1"

    headers = {
        "Content-Type": "application/json",
    }

    payload = {
        "inputs": sequence
    }

    print("Sending request to ESMFold...")
    response = requests.post(API_URL, headers=headers, json=payload)

    if response.status_code == 200:
        # Save PDB file
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print(f"[OK] Structure saved to: {output_file}")
        print("     This structure has all atoms including sidechains")
        return True
    elif response.status_code == 503:
        print("Model is loading, please wait...")
        time.sleep(20)
        # Retry once
        response = requests.post(API_URL, headers=headers, json=payload)
        if response.status_code == 200:
            with open(output_file, 'wb') as f:
                f.write(response.content)
            print(f"[OK] Structure saved to: {output_file}")
            return True
        else:
            print(f"Error: {response.status_code}")
            print(response.text)
            return False
    else:
        print(f"Error: {response.status_code}")
        print(response.text)
        return False

if __name__ == "__main__":
    # AviTag sequence
    sequence = "GLNDIFEAQKIEWHE"

    output_pdb = "C:/Users/bmartinez/ProteinDesign/AviTag_Binders_2025/02_scaffolds/approach_C_peptide_ensemble/00_initial_structure/avitag_esmfold.pdb"

    print("="*60)
    print("GENERATING AVITAG STRUCTURE WITH ESMFOLD")
    print("="*60)
    print()

    success = generate_structure_esmfold(sequence, output_pdb)

    if success:
        print()
        print("="*60)
        print("SUCCESS!")
        print("="*60)
        print()
        print("Structure generated with all atoms (backbone + sidechains)")
        print("Ready for GROMACS processing!")
        print()
        print("NEXT STEPS:")
        print("1. Use this structure with GROMACS pdb2gmx")
        print("2. Should work without histidine issues")
    else:
        print()
        print("="*60)
        print("FAILED")
        print("="*60)
        print("Could not generate structure via ESMFold")
        print("Alternative: Use local structure generation tool")
