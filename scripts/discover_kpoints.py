
import subprocess

from _isotropy_paths import SCRATCH_DIR, configure_environment

configure_environment()

def run_iso(commands):
    full_cmd = "PAGE 10000\n" + commands + "\nQUIT\n"
    res = subprocess.run(["iso"], input=full_cmd, capture_output=True, text=True, cwd=SCRATCH_DIR)
    return res.stdout

def get_all_kpoints(sg_num):
    # iso command to show all special k-points
    cmd = f"VALUE PARENT {sg_num}\nSHOW KPOINT\nDISPLAY KPOINT"
    out = run_iso(cmd)
    print(f"ISO K-POINT OUTPUT for Space Group {sg_num}:")
    print(out)

get_all_kpoints("225")
