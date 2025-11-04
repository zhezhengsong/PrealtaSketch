import itertools, os, re, subprocess, csv, statistics, multiprocessing as mp

PROG = os.environ.get("PROG", "./CountSketch_v2")

S_LIST = [64, 128, 256, 512, 1024]
K_LIST = [5, 10, 15, 20, 30]
R_LIST = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

REPEATS = 3  
OUT_CSV = f"results_{os.getpid()}.csv"
PAT = re.compile(r"Clustering Accuracy\s*=\s*([0-9]*\.?[0-9]+)")

def run_once(s, k, r):
    try:
        out = subprocess.check_output([PROG, str(s), str(k), str(r)], text=True)
        m = PAT.search(out)
        return float(m.group(1)) if m else float("nan")
    except Exception:
        return float("nan")

def run_avg(args):
    s, k, r = args
    accs = [run_once(s, k, r) for _ in range(REPEATS)]
    accs = [a for a in accs if a == a]
    mean = statistics.mean(accs) if accs else float("nan")
    std  = statistics.pstdev(accs) if len(accs) > 1 else 0.0
    return (s, k, r, mean, std, len(accs))

if __name__ == "__main__":
    grid = list(itertools.product(S_LIST, K_LIST, R_LIST))
    with mp.Pool(mp.cpu_count()) as pool:
        rows = list(pool.map(run_avg, grid))

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s","k","res","acc_mean","acc_std","n_runs"])
        w.writerows(rows)

    rows_valid = [r for r in rows if r[3] == r[3]]
    rows_valid.sort(key=lambda x: x[3], reverse=True)

    print("Best (top 10):")
    for r in rows_valid[:10]:
        s,k,res,mean,std,n = r
        print(f"s={s:>4}  k={k:>3}  res={res:>4}  acc_mean={mean:.4f}  (std={std:.4f}, n={n})")
    print(f"\nSaved CSV -> {OUT_CSV}")
