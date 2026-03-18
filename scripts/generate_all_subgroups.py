
import os
import re
import shutil
import subprocess
import sys

from _isotropy_paths import SCRATCH_DIR, SUBGROUPS_DIR, configure_environment, resolve_input_path

configure_environment()

def run_iso(commands):
    full_cmd = "PAGE 10000\n" + commands + "\nQUIT\n"
    res = subprocess.run(["iso"], input=full_cmd, capture_output=True, text=True, cwd=SCRATCH_DIR)
    return res.stdout

def get_kpoints(parent_sg):
    cmd = f"VALUE PARENT {parent_sg}\nSHOW KPOINT\nDISPLAY KPOINT"
    out = run_iso(cmd)
    kpoints = []
    lines = out.splitlines()
    found_header = False
    for line in lines:
        if "k vector" in line:
            found_header = True
            continue
        if found_header and line.strip():
            match = re.search(r"([A-Z]+)\s+\(", line.strip())
            if match: kpoints.append(match.group(1))
    return kpoints

def get_irreps(parent_sg, kpoint):
    cmd = f"VALUE PARENT {parent_sg}\nVALUE KPOINT {kpoint}\nSHOW IRREP\nDISPLAY IRREP"
    out = run_iso(cmd)
    irreps = []
    lines = out.splitlines()
    found_header = False
    for line in lines:
        if "Irrep" in line:
            found_header = True
            continue
        if found_header and line.strip():
            parts = line.strip().split()
            if parts: irreps.append(parts[0])
    return irreps

def get_subgroups_detailed(parent_sg, kpoint, irrep):
    kv_cmd = "VALUE KVALUE 1,1/4\n" if kpoint in ["DT", "LD", "SM", "GP"] else ""
    cmd = f"VALUE PARENT {parent_sg}\n{kv_cmd}VALUE KPOINT {kpoint}\nVALUE IRREP {irrep}\nSHOW SUBGROUP\nSHOW DIRECTION\nSHOW BASIS\nSHOW ORIGIN\nDISPLAY ISOTROPY"
    out = run_iso(cmd)
    subgroups = []
    lines = out.splitlines()
    found_header = False
    for line in lines:
        if "Subgroup" in line and "Dir" in line:
            found_header = True
            continue
        if found_header and line.strip() and not ("Enter RETURN" in line or "Quit display" in line):
            # Format: SG_NUM SYMBOL DIR BASIS ORIGIN
            # Example: 123 P4/mmm P1 (0,0,1),(1,0,0),(0,1,0) (0,0,0)
            match = re.search(r"(\d+)\s+([A-Za-z0-9/_-]+)\s+([A-Z0-9\(\)]+)\s+([-\d.,\(\)/ ]+)\s+([-\d.,\(\)/ ]+)", line)
            if match:
                subgroups.append({
                    'sg_num': match.group(1),
                    'symbol': match.group(2),
                    'dir': match.group(3),
                    'basis': match.group(4),
                    'origin': match.group(5)
                })
    return subgroups

def generate_subgroup_cif(parent_cif, subgroup_info, amplitude=0.0):
    # generate ideal subgroup CIF using findsym
    temp_in = SCRATCH_DIR / "temp_fs_ideal.in"
    findsym_cif = SCRATCH_DIR / "findsym.cif"
    try:
        with open(temp_in, "w") as f:
            subprocess.run(["findsym_cifinput", str(parent_cif)], stdout=f, cwd=SCRATCH_DIR)
        
        # We need to tell findsym the subgroup SG, Basis, and Origin
        # findsym supports keyword !targetSpaceGroup, !targetBasis, !targetOrigin
        with open(temp_in, "a") as f:
            f.write(f"\n!targetSpaceGroup {subgroup_info['sg_num']}\n")
            f.write(f"!targetBasis {subgroup_info['basis']}\n")
            f.write(f"!targetOrigin {subgroup_info['origin']}\n")
        
        # Run findsym
        subprocess.run(["findsym", str(temp_in)], input="\n0.001\n", capture_output=True, text=True, cwd=SCRATCH_DIR)
        
        if findsym_cif.exists():
            # If amplitude > 0, we should perturb? 
            # For now, let's just use this ideal CIF to get a signature
            return findsym_cif
    except:
        pass
    return None

def get_structure_signature(cif_file):
    temp_in = SCRATCH_DIR / "temp_sig_in"
    findsym_log = SCRATCH_DIR / "findsym.log"
    try:
        with open(temp_in, "w") as f:
            subprocess.run(["findsym_cifinput", str(cif_file)], stdout=f, cwd=SCRATCH_DIR)
        res = subprocess.run(["findsym", str(temp_in)], input="\n0.001\n", capture_output=True, text=True, cwd=SCRATCH_DIR)
        out = res.stdout
        if findsym_log.exists():
            with open(findsym_log, "r") as f:
                out += f.read()
        if "Space Group" not in out:
            return None
        sg_num = re.findall(r"Space Group:?\s+(\d+)", out)[0]
        lattice = re.findall(r"a,b,c,alpha,beta,gamma:\s+([\d.]+.*)", out)[-1].strip()
        pos_matches = re.findall(r"\s+\d+\s+[-]?[\d.]+\s+[-]?[\d.]+\s+[-]?[\d.]+\s+[\d.]+", out)
        pos_str = "\n".join(sorted(list(set(pos_matches))))
        return f"SG: {sg_num}|Lattice: {lattice}|Positions: {pos_str}"
    except:
        return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 generate_all_subgroups.py <parent_cif>")
        return
    parent_cif = resolve_input_path(sys.argv[1])
    parent_sg = "225" # assume for now or parse
    
    seen_signatures = set()
    total = 0
    kpoints = get_kpoints(parent_sg)
    
    for kp in kpoints:
        if kp not in ['GM', 'X', 'L', 'W']: continue
        for irrep in get_irreps(parent_sg, kp):
            subgroups = get_subgroups_detailed(parent_sg, kp, irrep)
            for sub in subgroups:
                clean_irrep = irrep.replace("*", "star").replace("+", "plus").replace("-", "minus")
                out_name = SUBGROUPS_DIR / f"subgroup_{kp}_{clean_irrep}_{sub['dir']}.cif"
                
                cif_file = generate_subgroup_cif(parent_cif, sub)
                if cif_file:
                    sig = get_structure_signature(cif_file)
                    if sig and sig not in seen_signatures:
                        seen_signatures.add(sig)
                        shutil.move(cif_file, out_name)
                        total += 1
                        print(f"[{total}] UNIQUE: {out_name.name} (SG {sig.split('|')[0]})")
                    else:
                        if cif_file.exists():
                            cif_file.unlink()
                
                findsym_cif = SCRATCH_DIR / "findsym.cif"
                findsym_log = SCRATCH_DIR / "findsym.log"
                if findsym_cif.exists():
                    findsym_cif.unlink()
                if findsym_log.exists():
                    findsym_log.unlink()

    print(f"\nDone! Generated {total} unique subgroup structures.")

if __name__ == "__main__":
    main()
