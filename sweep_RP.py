import itertools, os, re, subprocess, csv, statistics, multiprocessing as mp, math

# Program binary / executable
PROG = os.environ.get("PROG", "./RP_v2")

# Search space
S_LIST = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
K_LIST = [5, 10, 15, 20, 30]
R_LIST = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]   # clustering resolution

# --- Option A: static s2 candidates (recommended baseline) ---
S2_STATIC = [5, 8, 16, 24, 32, 48, 64, 96, 128]

# --- Option B: s2 tied to s (uncomment to use) ---
# def s2_candidates_for(s):
#     alphas = [0.5, 1.0, 2.0, 4.0]  # s2 ≈ alpha * s
#     return [max(1, int(alpha * s)) for alpha in alphas]

# Number of repeats for each (s,k,res,s2) to compute mean/std (override via env REPEATS)
REPEATS = int(os.environ.get("REPEATS", "10"))

OUT_CSV = f"results_{os.getpid()}.csv"
PAT = re.compile(r"Clustering Accuracy\s*=\s*([0-9]*\.?[0-9]+)")

def run_once(s, k, res, s2):
    """Run the RP program once and parse the clustering accuracy from stdout.
    Returns float('nan') on any error or if the pattern isn't found.
    """
    try:
        out = subprocess.check_output([PROG, str(s), str(k), str(res), str(s2)], text=True)
        m = PAT.search(out)
        return float(m.group(1)) if m else float('nan')
    except Exception:
        return float('nan')

def run_avg(args):
    """Repeat run_once REPEATS times and return mean/std/n for a grid point."""
    s, k, res, s2 = args
    vals = [run_once(s, k, res, s2) for _ in range(REPEATS)]
    vals = [v for v in vals if math.isfinite(v)]
    n = len(vals)
    if n == 0:
        mean = float('nan')
        std = float('nan')
    elif n == 1:
        mean = vals[0]
        std = 0.0
    else:
        mean = statistics.fmean(vals) if hasattr(statistics, 'fmean') else statistics.mean(vals)
        # Use population stdev to match many sweeps; if you prefer sample stdev, swap to statistics.stdev
        try:
            std = statistics.pstdev(vals)
        except Exception:
            # Fallback if pstdev not available
            mu = mean
            std = (sum((x - mu) ** 2 for x in vals) / n) ** 0.5
    return (s, k, res, s2, mean, std, n)

if __name__ == "__main__":
    # Build the grid
    # A) Static list:
    grid = list(itertools.product(S_LIST, K_LIST, R_LIST, S2_STATIC))

    # B) Relative-to-s:
    # grid = []
    # for s, k, res in itertools.product(S_LIST, K_LIST, R_LIST):
    #     for s2 in s2_candidates_for(s):
    #         grid.append((s, k, res, s2))

    # Evaluate
    with mp.Pool(mp.cpu_count()) as pool:
        rows = list(pool.map(run_avg, grid))

    # Write CSV: include mean/std and how many successful runs contributed
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s","k","res","s2","acc_mean","acc_std","n_runs"])
        w.writerows(rows)

    # Rank by acc_mean (ignore NaNs)
    rows_valid = [r for r in rows if r[4] == r[4]]
    rows_valid.sort(key=lambda x: x[4], reverse=True)

    print("Best (top 10):")
    for s,k,res,s2,mean,std,n in rows_valid[:10]:
        print(f"s={s:>4}  k={k:>3}  res={res:>4}  s2={s2:>4}  acc_mean={mean:.4f}  (std={std:.4f}, n={n})")
    print(f"\nSaved CSV -> {OUT_CSV}")
