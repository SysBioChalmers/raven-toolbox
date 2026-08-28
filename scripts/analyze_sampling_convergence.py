#!/usr/bin/env python3
"""Between-chain convergence for ACHR flux sampling (Gelman-Rubin R-hat).

``sampling.md``'s existing ``thinning`` result is a *single-chain* diagnostic
(lag-1 autocorrelation / effective sample size): it asks how independent
consecutive samples are within one Markov chain. It cannot detect a chain that
mixes fine locally but never reaches part of the flux polytope -- for that you
need multiple independent chains started from different points and a check
that they agree on the same distribution. That's what this script adds.

For each reaction, run ``n_chains`` independent ``random_sampling`` calls
(different seeds) and compute the Gelman-Rubin statistic:

    W     = mean within-chain variance
    B     = n/(m-1) * sum_over_chains((chain_mean - grand_mean)^2)
    var+  = (n-1)/n * W + B/n
    R-hat = sqrt(var+ / W)

R-hat close to 1.0 means the chains agree; R-hat > 1.1 (the common threshold)
or > 1.01 (the stricter one used for published MCMC work) flags a reaction
whose sampled distribution still depends on where the chain started -- i.e.
not yet trustworthy for that reaction. Constant reactions (W == 0, e.g. a
reaction pinned by tight bounds) are excluded rather than producing a
divide-by-zero.

Usage
-----
    python scripts/analyze_sampling_convergence.py \
        --model /path/to/yeast-GEM.xml --out work/ --n-chains 4 \
        --n-samples 300 --thinning 100 --warmup 1000

Chains run in parallel (``--workers``, default = ``n_chains``) via
``ProcessPoolExecutor`` -- each chain gets its own model copy in its own
process, since cobrapy models aren't guaranteed thread-safe for concurrent
solves. Results are cached per (model, n_chains, n_samples, thinning, warmup)
config under ``--out`` so a re-run with the same config is instant.
"""
from __future__ import annotations

import argparse
import pickle
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cobra
import numpy as np
import pandas as pd


def _run_one_chain(model_path: str, seed: int, n_samples: int, thinning: int, warmup: int, method: str) -> pd.DataFrame:
    # Re-imported in the child process (ProcessPoolExecutor pickles by reference,
    # not by closure), and the model is reloaded from disk rather than pickling a
    # cobra.Model across the process boundary (its solver interface doesn't pickle
    # cleanly on every backend).
    from raven_toolbox.analysis.sampling import random_sampling

    model = cobra.io.read_sbml_model(model_path)
    result = random_sampling(
        model, n_samples, method=method, seed=seed, thinning=thinning, warmup=warmup,
    )
    return result.samples


def _gelman_rubin(chains: list[pd.DataFrame]) -> pd.Series:
    """Per-reaction R-hat across chains of equal length. NaN where all chains are constant."""
    n = chains[0].shape[0]
    stacked = np.stack([c.to_numpy() for c in chains], axis=0)  # (m, n, r)
    chain_means = stacked.mean(axis=1)  # (m, r)
    grand_mean = chain_means.mean(axis=0)  # (r,)
    m = stacked.shape[0]

    within_var = stacked.var(axis=1, ddof=1)  # (m, r)
    W = within_var.mean(axis=0)  # (r,)
    B_over_n = ((chain_means - grand_mean) ** 2).sum(axis=0) / (m - 1)  # (r,)

    var_plus = (n - 1) / n * W + B_over_n
    with np.errstate(divide="ignore", invalid="ignore"):
        r_hat = np.sqrt(var_plus / W)
    return pd.Series(r_hat, index=chains[0].columns, name="r_hat")


def _cache_path(out: Path, model_name: str, n_chains: int, n_samples: int, thinning: int, warmup: int, method: str) -> Path:
    return out / f"chains_{model_name}_{method}_m{n_chains}_n{n_samples}_t{thinning}_w{warmup}.pkl"


def run(model_path: Path, out: Path, n_chains: int, n_samples: int, thinning: int, warmup: int, workers: int | None, method: str = "achr") -> pd.Series:
    out.mkdir(parents=True, exist_ok=True)
    model_name = model_path.stem
    cache = _cache_path(out, model_name, n_chains, n_samples, thinning, warmup, method)

    if cache.exists():
        print(f"[cache hit] {cache}")
        with open(cache, "rb") as fh:
            chains = pickle.load(fh)
    else:
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=workers or n_chains) as pool:
            futures = [
                pool.submit(_run_one_chain, str(model_path), seed, n_samples, thinning, warmup, method)
                for seed in range(n_chains)
            ]
            chains = [f.result() for f in futures]
        print(f"{n_chains} chains x {n_samples} samples ({method}) on {model_name}: {time.time() - t0:.1f}s wall")
        with open(cache, "wb") as fh:
            pickle.dump(chains, fh)

    r_hat = _gelman_rubin(chains)
    n_total = len(r_hat)
    n_const = r_hat.isna().sum()
    scored = r_hat.dropna()
    print(f"\n{model_name}: {n_total} reactions, {n_const} constant (excluded), {len(scored)} scored")
    print(f"  R-hat median={scored.median():.4f}  p90={scored.quantile(0.9):.4f}  max={scored.max():.4f}")
    print(f"  R-hat > 1.01: {(scored > 1.01).sum()} ({100 * (scored > 1.01).mean():.1f}%)")
    print(f"  R-hat > 1.1:  {(scored > 1.1).sum()} ({100 * (scored > 1.1).mean():.1f}%)")
    worst = scored.sort_values(ascending=False).head(10)
    print("  worst 10:")
    for rxn_id, val in worst.items():
        print(f"    {rxn_id}: {val:.4f}")
    return r_hat


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("work"))
    p.add_argument("--n-chains", type=int, default=4)
    p.add_argument("--n-samples", type=int, default=300)
    p.add_argument("--thinning", type=int, default=100)
    p.add_argument("--warmup", type=int, default=1000)
    p.add_argument("--workers", type=int, default=None, help="default: one per chain")
    p.add_argument("--method", choices=["achr", "chrr"], default="achr")
    args = p.parse_args()

    run(args.model, args.out, args.n_chains, args.n_samples, args.thinning, args.warmup, args.workers, args.method)
