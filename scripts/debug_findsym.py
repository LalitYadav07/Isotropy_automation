
import subprocess
import re

from _isotropy_paths import SCRATCH_DIR, configure_environment

configure_environment()

def run_iso(commands):
    full_cmd = "PAGE 10000\n" + commands + "\nQUIT\n"
    res = subprocess.run(["iso"], input=full_cmd, capture_output=True, text=True, cwd=SCRATCH_DIR)
    return res.stdout

def get_subgroups_table(parent_sg, kp, irrep):
    cmd = f"VALUE PARENT {parent_sg}\nVALUE KPOINT {kp}\nVALUE IRREP {irrep}\nSHOW SUBGROUP\nSHOW DIRECTION\nSHOW BASIS\nSHOW ORIGIN\nDISPLAY ISOTROPY"
    out = run_iso(cmd)
    print(f"ISO TABLE (with Dir):\n{out}")

get_subgroups_table("225", "X", "X5-")
