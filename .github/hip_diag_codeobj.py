#!/usr/bin/env python3
# TEMP (throwaway branch): print resource usage of the run_newton kernels embedded in an executable.
# roc-obj-extract did not produce files on the runners, so carve the code objects out by offset.
import re, subprocess, sys
exe = sys.argv[1]
ls = subprocess.run(["roc-obj-ls", exe], capture_output=True, text=True).stdout
keys = ("agpr_count", "vgpr_count", "sgpr_count", "max_flat_workgroup_size", "private_segment_fixed_size",
        "group_segment_fixed_size", "uses_dynamic_stack", "vgpr_spill_count", "sgpr_spill_count",
        "wavefront_size", "kernarg_segment_size", ".name:")
n = 0
with open(exe, "rb") as f:
    for line in ls.splitlines():
        m = re.search(r"(hipv4-amdgcn\S*)\s+file://\S*#offset=(\d+)&size=(\d+)", line)
        if not m:
            continue
        off, size = int(m.group(2)), int(m.group(3))
        f.seek(off); blob = f.read(size)
        if b"run_newton" not in blob:
            continue
        fn = f"/tmp/co_run_newton_{n}.o"; n += 1
        open(fn, "wb").write(blob)
        print(f"=== {fn} ({m.group(1)}, {size} bytes)")
        notes = subprocess.run(["/opt/rocm/lib/llvm/bin/llvm-readobj", "--notes", fn], capture_output=True, text=True).stdout
        sel = [l.strip()[:220] for l in notes.splitlines() if any(k in l for k in keys)]
        for i, l in enumerate(sel):
            if ".name:" in l and "run_newton" in l:
                print("\n".join(sel[max(0, i - 10): i + 10])); print("---")
        for l in notes.splitlines():
            if re.search(r"amdhsa\.(version|target)|code_object", l):
                print(l.strip()[:200])
print(f"carved {n} code objects containing run_newton")
