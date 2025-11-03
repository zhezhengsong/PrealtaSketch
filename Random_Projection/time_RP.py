# import re, statistics, subprocess

# PROG = "./RP_v2_time"
# ARGS = ["1024", "15", "1.0", "96"]
# REPEATS = 10
# pat = re.compile(r"CountSketch_ms\s*=\s*([0-9]*\.?[0-9]+)")

# times = []
# for _ in range(REPEATS):
#     out = subprocess.check_output([PROG, *ARGS], text=True)
#     m = pat.search(out)
#     if m:
#         times.append(float(m.group(1)))

# print(f"runs={len(times)}  mean_ms={statistics.mean(times):.3f}  "
#       f"std_ms={statistics.pstdev(times) if len(times)>1 else 0.0:.3f}  "
#       f"min_ms={min(times):.3f}  max_ms={max(times):.3f}")

import re, statistics, subprocess, collections

PROG = "./RP_v2_time"
ARGS = ["1024", "15", "1.0", "96"]
REPEATS = 10

pat = re.compile(r"([A-Za-z0-9_+\-]+)_ms\s*=\s*([0-9]*\.?[0-9]+)")

buckets = collections.defaultdict(list)

for _ in range(REPEATS):
    out = subprocess.check_output([PROG, *ARGS], text=True)
    for name, val in pat.findall(out):
        buckets[name].append(float(val))

for name, vals in buckets.items():
    print(f"[{name}] runs={len(vals)}  mean_ms={statistics.mean(vals):.3f}  "
          f"std_ms={statistics.pstdev(vals) if len(vals)>1 else 0.0:.3f}  "
          f"min_ms={min(vals):.3f}  max_ms={max(vals):.3f}")
