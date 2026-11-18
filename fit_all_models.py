"""Optional: pre-fit every Bayesian model used in the paper, several at a time.

The notebooks fit any model that is not cached yet, one after another (each fit
takes ~2-5 min, ~46 fits in total). Running this script first fills the same
cache in parallel, which is much faster on a multi-core machine:

    uv run python fit_all_models.py --jobs 4

Each job samples its 4 chains sequentially (1 core); PyMC seeds every chain from
random_seed = 42, so the draws are identical to a run with parallel chains. Every
job holds one ~360 MB log-likelihood in memory while it computes PSIS-LOO, so
allow ~2 GB of RAM per job.
"""
from __future__ import annotations

import argparse
import time

from joblib import Parallel, delayed

import frogcall as fc

K_VALUES = [2, 3, 4, 5]                       # Table S4 / Fig. S4
LADDER = [m for m, _ in fc.LADDER_BLOCKS]     # model_0 ... model_4_ord


def jobs():
    out = []
    for sp in fc.SPECIES:
        # Table 1, Table 2, Figs 2, S1-S3: the K = 2 ladder plus the two other
        # heterospecific-calling encodings.
        for m in LADDER + ["model_4_bin", "model_4_lag1"]:
            out.append((sp, m, 2, 1.0))
        # Table S3: model_4_ord under the tighter N(0, 0.5) coefficient prior.
        out.append((sp, "model_4_ord", 2, 0.5))
        # Table S4 / Fig. S4: the ladder at K = 3, 4, 5 harmonics.
        for K in K_VALUES:
            if K != 2:
                for m in LADDER:
                    out.append((sp, m, K, 1.0))
    return out


def run_one(sp, m, K, sigma):
    df = fc.load_data()
    d = fc.prepare_species_df(df, sp)
    res = fc.fit_cached(fc.ladder_formula(m, K), d, sp, f"{m}_K{K}", prior_sigma=sigma, cores=1)
    return sp, m, K, sigma, res["loo"]["elpd_loo"], res["seconds"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4, help="models fitted in parallel")
    args = ap.parse_args()
    todo = jobs()
    print(f"{len(todo)} models; fitting with {args.jobs} parallel jobs")
    t0 = time.time()
    out = Parallel(n_jobs=args.jobs, verbose=5)(delayed(run_one)(*j) for j in todo)
    for sp, m, K, sigma, elpd, secs in out:
        print(f"{fc.SP_SHORT[sp]:<17} {m:<13} K={K} N(0,{sigma:g}): elpd_loo = {elpd:9.1f}"
              + ("" if secs != secs else f"  ({secs / 60:.1f} min)"))
    print(f"total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
