
import os
import re
import shutil
import subprocess
import sys

from _isotropy_paths import SCRATCH_DIR, configure_environment, resolve_input_path

configure_environment()

def get_structure_signature(cif_file, index):
    """
    Identifies the standardized symmetry and structure using findsym.
    """
    # Use unique filenames for parallel-like safety
    temp_in = SCRATCH_DIR / f"temp_comp_{index}.in"
    temp_log = SCRATCH_DIR / f"findsym_{index}.log"
    findsym_log = SCRATCH_DIR / "findsym.log"
    
    # 1. Preprocess with findsym_cifinput
    with open(temp_in, "w") as f:
        subprocess.run(["findsym_cifinput", str(cif_file)], stdout=f, cwd=SCRATCH_DIR)
    
    # 2. Run findsym
    # We pass the input file. If findsym works on stdin, we can also pipe.
    res = subprocess.run(["findsym", str(temp_in)], input="\n0.001\n", capture_output=True, text=True, cwd=SCRATCH_DIR)
    out = res.stdout
    
    # Check findsym.log which is usually created in the CWD
    if findsym_log.exists():
        with open(findsym_log, "r") as f:
            out += f.read()
        # Clean up or move
        shutil.move(findsym_log, temp_log)
    
    if "Space Group" not in out:
        # Try one more time with simple stdin if !useKeyWords in temp_in failed
        with open(temp_in, "r") as f:
             res = subprocess.run(["findsym"], input=f.read() + "\n0.001\n", capture_output=True, text=True, cwd=SCRATCH_DIR)
             out = res.stdout
             if findsym_log.exists():
                 with open(findsym_log, "r") as f2:
                     out += f2.read()

    if "Space Group" not in out:
        return None

    # 3. Extract key features for comparison
    sg_num = re.findall(r"Space Group:?\s+(\d+)", out)[0]
    lattice_matches = re.findall(r"a,b,c,alpha,beta,gamma:\s+([\d.]+.*)", out)
    lattice = lattice_matches[-1].strip() if lattice_matches else ""
    
    # Extract positions (standardized)
    # We look for the final Atomic positions section
    pos_matches = re.findall(r"\s+\d+\s+[-]?[\d.]+\s+[-]?[\d.]+\s+[-]?[\d.]+\s+[\d.]+", out)
    # Filter out possible duplicates or headers if any, and keep unique ones
    positions = sorted(list(set(pos_matches)))
    pos_str = "\n".join(positions)
    
    signature = f"SG: {sg_num}\nLattice: {lattice}\nPositions:\n{pos_str}"
    # Debug: save signature for inspection
    with open(SCRATCH_DIR / f"sig_{index}.txt", "w") as f:
        f.write(signature)
        
    return signature

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_cifs.py <cif1> <cif2>")
        return
    
    cif1 = resolve_input_path(sys.argv[1])
    cif2 = resolve_input_path(sys.argv[2])
    
    sig1 = get_structure_signature(cif1, 1)
    sig2 = get_structure_signature(cif2, 2)
    
    if sig1 is None or sig2 is None:
        print("Error: Could not identify symmetry for one or both files.")
        return

    if sig1 == sig2:
        print(f"\n[MATCH] {cif1.name} and {cif2.name} are IDENTICAL.")
    else:
        print(f"\n[DIFFERENT] {cif1.name} and {cif2.name} are DISTINCT.")
        l1, l2 = sig1.split('\n'), sig2.split('\n')
        if l1[0] != l2[0]: print(f"  Reason: {l1[0]} vs {l2[0]}")
        elif l1[1] != l2[1]: print(f"  Reason: Different lattice ({l1[1]} vs {l2[1]})")
        else: print("  Reason: Different atomic arrangements.")

if __name__ == "__main__":
    main()
