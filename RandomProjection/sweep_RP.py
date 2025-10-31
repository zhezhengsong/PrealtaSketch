import itertools, os, re, subprocess, csv, multiprocessing as mp, math

PROG = os.environ.get("PROG", "./RP_v2")

# Search space
S_LIST = [64, 128, 256, 512, 1024]
K_LIST = [5, 10, 15, 20, 30]
R_LIST = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]   # clustering resolution

# --- Option A: static s2 candidates (recommended baseline) ---
S2_STATIC = [3, 5, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256]

# --- Option B: s2 tied to s (uncomment to use) ---
# def s2_candidates_for(s):
#     alphas = [0.5, 1.0, 2.0, 4.0]  # s2 ≈ alpha * s
#     out = {max(8, int(round(a * s))) for a in alphas}
#     return sorted(out)

OUT_CSV = f"results_{os.getpid()}.csv"
PAT = re.compile(r"Clustering Accuracy\s*=\s*([0-9]*\.?[0-9]+)")

def run_once(s, k, res, s2):
    try:
        out = subprocess.check_output([PROG, str(s), str(k), str(res), str(s2)], text=True)
        m = PAT.search(out)
        return float(m.group(1)) if m else float("nan")
    except Exception:
        return float("nan")

def run_row(args):
    s, k, res, s2 = args
    acc = run_once(s, k, res, s2)
    return (s, k, res, s2, acc)

if __name__ == "__main__":
    # --- choose ONE of these two grids ---

    # A) Static list:
    grid = list(itertools.product(S_LIST, K_LIST, R_LIST, S2_STATIC))

    # B) Relative-to-s:
    # grid = []
    # for s, k, res in itertools.product(S_LIST, K_LIST, R_LIST):
    #     for s2 in s2_candidates_for(s):
    #         grid.append((s, k, res, s2))

    with mp.Pool(mp.cpu_count()) as pool:
        rows = list(pool.map(run_row, grid))

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s","k","res","s2","acc"])
        w.writerows(rows)

    rows_valid = [r for r in rows if r[4] == r[4]]
    rows_valid.sort(key=lambda x: x[4], reverse=True)

    print("Best (top 15):")
    for s,k,res,s2,acc in rows_valid[:15]:
        print(f"s={s:>4}  k={k:>3}  res={res:>3}  s2={s2:>4}  acc={acc:.4f}")
    print(f"\nSaved CSV -> {OUT_CSV}")
