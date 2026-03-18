
import subprocess

from _isotropy_paths import SCRATCH_DIR, configure_environment

configure_environment()

def run_iso(commands):
    full_cmd = "PAGE 10000\n" + commands + "\nQUIT\n"
    res = subprocess.run(["iso"], input=full_cmd, capture_output=True, text=True, cwd=SCRATCH_DIR)
    return res.stdout

def debug_l_point():
    cmd = "VALUE PARENT 225\nVALUE KPOINT L\nSHOW IRREP\nSHOW SUBGROUP\nDISPLAY SUBGROUP"
    out = run_iso(cmd)
    print("ISO SUBGROUP OUTPUT for L-point:")
    print(out)

debug_l_point()
