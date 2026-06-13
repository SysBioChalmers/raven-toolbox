"""Random-objective flux sampling — RAVEN's ``randomSampling`` (port + improvements).

Samples the flux solution space by the **random-objective** method of Bordel et al.
(2010, PLoS Comput Biol, doi:10.1371/journal.pcbi.1000859), as ported from RAVEN's
``randomSampling``: each sample maximises a small random linear combination of
reactions, so every sample is an *extreme point* (vertex) of the flux polytope.

This is a different statistical object from cobrapy's ``cobra.sampling`` (OptGP /
ACHR), which draw a (near-)uniform Markov-chain sample of the polytope *interior*.
Use cobra's samplers when you need the uniform flux distribution; use this when you
want a fast, robust spread of diverse optimal states — the workflow RAVEN uses to
compare conditions, and one that stays well-behaved on large or tightly-constrained
models where MCMC mixing is poor. cobrapy has no equivalent, so this is a genuine
addition, not a wrapper.

Improvements over RAVEN (see IMPROVEMENTS SAMP1):

* **`good_reactions` via one FVA pass**, not a hand-rolled per-reaction ``parfor``
  loop. A reaction is usable as a random objective if it can carry flux and is not
  stuck in a stoichiometrically-infeasible loop (its range blows past the arbitrary
  large bound). ``cobra``'s FVA computes exactly that, faster and in far less code,
  and can optionally be made ``loopless``.
* **Reproducible** via ``seed`` (RAVEN has no seed control).
* **`n_objectives` is a parameter** (RAVEN hard-codes 2, though its docstring claims
  3).
* **Tidy output**: a ``samples`` DataFrame shaped samples × reactions (matching
  ``cobra.sampling``), plus the reusable ``good_reactions`` list — instead of a
  reactions × samples matrix and a parallel index vector.
"""
from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass

import cobra
import numpy as np
import pandas as pd
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis, pfba
from optlang.symbolics import add

logger = logging.getLogger(__name__)


@dataclass
class RandomSamplingResult:
    """Output of :func:`random_sampling`.

    ``samples`` is a DataFrame of flux vectors shaped *n_samples × n_reactions*
    (one sample per row, reaction ids as columns — the ``cobra.sampling`` layout).
    ``good_reactions`` is the list of reaction ids that were eligible as random
    objectives; pass it back in to skip the (one-off) FVA on a repeat run.
    """

    samples: pd.DataFrame
    good_reactions: list[str]


def find_good_reactions(
    model: cobra.Model,
    *,
    flux_tol: float = 1e-9,
    loopless: bool = True,
    exclude_reactions: Iterable[str] | None = None,
) -> list[str]:
    """Reactions usable as random objectives: carry real (non-loop) flux.

    A reaction is kept if its FVA range spans more than ``flux_tol``. With
    ``loopless`` (default) the FVA is loopless (``cycleFreeFlux``), so reactions
    that can carry flux *only* through a thermodynamically-infeasible cycle have a
    ~0 loopless range and are dropped — the right test for "loopy", unlike a fixed
    bound threshold which wrongly drops legitimate reactions that simply reach the
    model's default (e.g. 1000) bound. Pass ``loopless=False`` for a faster, looser
    pass that keeps any flux-carrying reaction (loops included).
    """
    fva = flux_variability_analysis(
        model, fraction_of_optimum=0.0,
        loopless="cycleFreeFlux" if loopless else None,
    )
    excluded = set(exclude_reactions or ())
    return [
        rxn_id
        for rxn_id, lo, hi in zip(fva.index, fva["minimum"], fva["maximum"], strict=True)
        if rxn_id not in excluded and max(abs(lo), abs(hi)) > flux_tol
    ]


def random_sampling(
    model: cobra.Model,
    n_samples: int = 1000,
    *,
    n_objectives: int = 2,
    good_reactions: Iterable[str] | None = None,
    replace_max_bound: bool = False,
    min_flux: bool = False,
    loopless_good_reactions: bool = True,
    exclude_reactions: Iterable[str] | None = None,
    max_attempts: int = 100,
    suppress_errors: bool = False,
    seed: int | None = None,
) -> RandomSamplingResult:
    """Random-objective sampling of ``model``'s flux space (Bordel et al. 2010).

    Each sample maximises ``sum(w_i * v_i)`` over ``n_objectives`` reactions drawn at
    random from ``good_reactions``, with weights ``w_i = U(0,1) * (±1)`` (a random
    sign per reaction, as in RAVEN). The resulting flux vector is one sample.

    Parameters
    ----------
    n_samples
        Number of flux vectors to return.
    n_objectives
        Reactions combined into each random objective (RAVEN's fixed 2).
    good_reactions
        Reaction ids eligible as objectives. If ``None`` they are computed once with
        :func:`find_good_reactions` and returned for reuse.
    replace_max_bound
        RAVEN's ``replaceBoundsWithInf``: replace the largest upper bound with
        ``+inf`` (and the smallest negative lower bound with ``-inf``) before
        sampling, so a reaction whose biological maximum exceeds the model's
        arbitrary cap is not pinned at it. **Off by default** — unlike RAVEN. It
        applies only to the sampling phase (``good_reactions`` is always found on
        the finite bounds), and it can open unbounded directions through loops
        that show up as large fluxes in non-objective reactions; pair it with
        ``min_flux`` if you enable it.
    min_flux
        After maximising the random objective, re-solve parsimoniously
        (:func:`cobra.flux_analysis.pfba`) to minimise total flux at that optimum —
        squeezes residual loops out of each individual sample.
    loopless_good_reactions, exclude_reactions
        Forwarded to :func:`find_good_reactions` when it is invoked (loopless loop
        detection is on by default).
    max_attempts, suppress_errors
        A sample is retried if the random objective is degenerate (zero flux). After
        ``max_attempts`` failures this raises, unless ``suppress_errors`` (then the
        degenerate solution is kept with a warning).
    seed
        Seed for reproducible objective draws.

    Returns
    -------
    RandomSamplingResult
    """
    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")
    rng = np.random.default_rng(seed)
    model = model.copy()

    if model.slim_optimize(error_value=None) is None:
        raise ValueError(
            "The model has no feasible solution, likely due to incompatible constraints."
        )

    # good_reactions must be found on the finite bounds (FVA cannot handle inf),
    # before any bound replacement.
    if good_reactions is None:
        good_reactions = find_good_reactions(
            model, loopless=loopless_good_reactions,
            exclude_reactions=exclude_reactions,
        )
    good_reactions = list(good_reactions)

    if replace_max_bound:
        max_ub = max(r.upper_bound for r in model.reactions)
        min_lb = min(r.lower_bound for r in model.reactions)
        for r in model.reactions:
            if r.upper_bound == max_ub:
                r.upper_bound = float("inf")
            if min_lb < 0 and r.lower_bound == min_lb:
                r.lower_bound = float("-inf")

    if len(good_reactions) < n_objectives:
        raise ValueError(
            f"Only {len(good_reactions)} usable reactions found, need at least "
            f"n_objectives={n_objectives}. Check the model's constraints."
        )

    good_rxn_objs = [model.reactions.get_by_id(r) for r in good_reactions]
    reaction_ids = [r.id for r in model.reactions]
    samples = np.zeros((n_samples, len(reaction_ids)))

    for i in range(n_samples):
        for attempt in range(1, max_attempts + 1):
            chosen = rng.choice(len(good_rxn_objs), size=n_objectives, replace=False)
            signs = rng.choice((-1.0, 1.0), size=n_objectives)
            weights = rng.random(n_objectives) * signs
            terms = [w * good_rxn_objs[j].flux_expression
                     for j, w in zip(chosen, weights, strict=True)]
            # add() (not sum()) builds the symbolic objective in one pass; sum()
            # re-canonicalises the optlang expression on every term (O(n^2)).
            model.objective = model.problem.Objective(add(terms), direction="max")
            sol = model.optimize()
            if sol.status == "optimal" and abs(sol.objective_value) > 1e-8:
                fluxes = (pfba(model) if min_flux else sol).fluxes.reindex(reaction_ids)
                if fluxes.isna().any():  # solver returned an unexpected reaction set
                    missing = fluxes.index[fluxes.isna()].tolist()
                    raise OptimizationError(
                        "solver returned fluxes missing reaction(s) "
                        f"{missing[:5]}; cannot assemble a NaN-free sample matrix."
                    )
                samples[i, :] = fluxes.to_numpy()
                break
            if attempt == max_attempts:
                if not suppress_errors:
                    raise OptimizationError(
                        "Could not find a non-zero, loop-free solution after "
                        f"{max_attempts} attempts for sample {i}. Review the model's "
                        "constraints, or set suppress_errors=True."
                    )
                logger.warning("Sample %d: kept a degenerate solution after %d attempts.",
                               i, max_attempts)
                samples[i, :] = sol.fluxes.reindex(reaction_ids).to_numpy()

    return RandomSamplingResult(
        samples=pd.DataFrame(samples, columns=reaction_ids),
        good_reactions=good_reactions,
    )
