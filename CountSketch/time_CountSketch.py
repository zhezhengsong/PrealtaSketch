import re, statistics, subprocess

PROG = "./CountSketch_v2"
ARGS = ["1024", "15", "1.0"]
REPEATS = 10
pat = re.compile(r"CountSketch_ms\s*=\s*([0-9]*\.?[0-9]+)")

times = []
for _ in range(REPEATS):
    out = subprocess.check_output([PROG, *ARGS], text=True)
    m = pat.search(out)
    if m:
        times.append(float(m.group(1)))

print(f"runs={len(times)}  mean_ms={statistics.mean(times):.3f}  "
      f"std_ms={statistics.pstdev(times) if len(times)>1 else 0.0:.3f}  "
      f"min_ms={min(times):.3f}  max_ms={max(times):.3f}")
