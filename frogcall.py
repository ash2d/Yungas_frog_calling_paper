"""Shared helpers for the notebooks: paths, data preparation and cached Bayesian fits.

Every Bayesian model in the paper is a bambi cumulative-logit model fitted with the
same sampler settings (4 chains x 1000 draws after 1000 tuning steps,
target_accept = 0.9, random_seed = 42). Fits are cached in ``outputs/`` so each
model is sampled once and then shared by every notebook:

* ``outputs/models/<name>.nc``   posterior + sample_stats (no log-likelihood, a few MB)
* ``outputs/elpd/<name>.npz``    pointwise PSIS-LOO elpd and the LOO summary

The cache file name contains a hash of the formula, the prior SD and the sampler
settings, so a changed model is never confused with a cached one. Delete
``outputs/`` to refit everything from scratch.
"""
from __future__ import annotations

import hashlib
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "data" / "df_helechos_with2020.csv"
AUDIO_DIR = ROOT / "data" / "audio"
FIG_DIR = ROOT / "figures"
TABLE_DIR = ROOT / "tables"
OUT_DIR = ROOT / "outputs"
MODEL_DIR = OUT_DIR / "models"
ELPD_DIR = OUT_DIR / "elpd"
CV_DIR = OUT_DIR / "cv"
for _d in (FIG_DIR, TABLE_DIR, MODEL_DIR, ELPD_DIR, CV_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ── Study constants ──────────────────────────────────────────────────────────
SPECIES = ["Gastrotheca chrysosticta", "Oreobates berdemenos"]
SP_SHORT = {
    "Gastrotheca chrysosticta": "G. chrysosticta",
    "Oreobates berdemenos": "O. berdemenos",
}
TEMP_COL, RH_COL = "Temp", "RH%"

# ── Sampler settings (identical for every model in the paper) ────────────────
BAYES_DRAWS = 1000
BAYES_TUNE = 1000
BAYES_CHAINS = 4
BAYES_TARGET_ACCEPT = 0.9
BAYES_SEED = 42
# Chains run in parallel when cores > 1. PyMC derives one seed per chain from
# BAYES_SEED, so the draws do not depend on the number of cores.
BAYES_CORES = int(os.environ.get("FROGCALL_CORES", min(4, os.cpu_count() or 1)))

# ── The nested model ladder ──────────────────────────────────────────────────
LADDER_BLOCKS = [
    ("model_0", []),
    ("model_1", ["temp_z", "rh_z"]),
    ("model_2", ["I(temp_z ** 2)", "I(rh_z ** 2)"]),
    ("model_3", ["temp_z:rh_z"]),
    ("model_4_ord", ["other_calling_ord"]),
]
M4_VARIANT_TERM = {"model_4_bin": "other_calling_bin", "model_4_lag1": "other_calling_lag1"}


def fourier_terms(K: int) -> List[str]:
    """sin/cos pairs for harmonics 1..K of the diel and then the seasonal cycle."""
    diel = [t for j in range(1, K + 1) for t in (f"sin_h{j}", f"cos_h{j}")]
    seas = [t for j in range(1, K + 1) for t in (f"sin_d{j}", f"cos_d{j}")]
    return diel + seas


def ladder_terms(model: str, K: int = 2) -> List[str]:
    if model in M4_VARIANT_TERM:
        return ladder_terms("model_3", K) + [M4_VARIANT_TERM[model]]
    terms = fourier_terms(K)
    for name, block in LADDER_BLOCKS:
        terms = terms + block
        if name == model:
            return terms
    raise KeyError(model)


def ladder_formula(model: str, K: int = 2, response: str = "response_cat") -> str:
    return f"{response} ~ " + " + ".join(ladder_terms(model, K))


# ── Data ─────────────────────────────────────────────────────────────────────
def load_data() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH)
    df["DateTime"] = pd.to_datetime(df["DateTime"], errors="coerce")
    return df.set_index("DateTime").sort_index()


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=0)


def prepare_species_df(df_in: pd.DataFrame, response_col: str, k_max: int = 5) -> pd.DataFrame:
    """Model frame for one focal species (same rows for every harmonic count).

    Rows with a missing value in any model column are dropped; this includes the
    lag-1 heterospecific term, so every model is fitted to the same 11,302 hours.
    """
    d = df_in.copy()
    d["temp_z"] = zscore(pd.to_numeric(d[TEMP_COL], errors="coerce"))
    d["rh_z"] = zscore(pd.to_numeric(d[RH_COL], errors="coerce"))
    d["hour"] = d.index.hour
    d["doy"] = d.index.dayofyear
    two_pi = 2 * np.pi
    for j in range(1, k_max + 1):
        d[f"sin_h{j}"] = np.sin(j * two_pi * d["hour"] / 24.0)
        d[f"cos_h{j}"] = np.cos(j * two_pi * d["hour"] / 24.0)
        d[f"sin_d{j}"] = np.sin(j * two_pi * d["doy"] / 366.0)
        d[f"cos_d{j}"] = np.cos(j * two_pi * d["doy"] / 366.0)
    other_col = [c for c in SPECIES if c != response_col][0]
    d["other_calling_ord"] = pd.to_numeric(d[other_col], errors="coerce")
    d["other_calling_bin"] = (d["other_calling_ord"] > 0).astype(int)
    d["other_calling_lag1"] = d["other_calling_bin"].shift(1)
    d["response_raw"] = pd.to_numeric(d[response_col], errors="coerce")
    harmonic_cols = [f"{f}_{c}{j}" for c in ("h", "d") for j in range(1, k_max + 1)
                     for f in ("sin", "cos")]
    cols = ["response_raw", "temp_z", "rh_z", "other_calling_ord", "other_calling_bin",
            "other_calling_lag1", "hour", "doy"] + harmonic_cols
    dm = d[cols].dropna()
    cat = pd.Categorical(dm["response_raw"], ordered=True)
    dm["response_code"] = cat.codes
    dm["response_cat"] = cat
    return dm


# ── Cached Bayesian fitting ──────────────────────────────────────────────────
def cache_name(formula: str, species: str, label: str, prior_sigma: float = 1.0) -> str:
    key = "|".join([formula, f"{prior_sigma:g}", str(BAYES_DRAWS), str(BAYES_TUNE),
                    str(BAYES_CHAINS), str(BAYES_TARGET_ACCEPT), str(BAYES_SEED)])
    h = hashlib.sha1(key.encode()).hexdigest()[:8]
    sd = f"{prior_sigma:g}".replace(".", "p")
    return f"{label}_sd{sd}_{species.replace(' ', '_')}_{h}"


def strip_heavy(idata):
    """Drop the ~360 MB log_likelihood (and any posterior_predictive) group."""
    slim = idata.copy()
    for grp in ("log_likelihood", "posterior_predictive"):
        if grp in slim.children:
            del slim[grp]
    return slim


def psis_loo(idata):
    import arviz as az
    res = az.loo(idata, pointwise=True)
    elpd_i = np.asarray(res.elpd_i, dtype=float).ravel()
    k = np.asarray(res.pareto_k, dtype=float).ravel()
    good_k = float(getattr(res, "good_k", 0.7) or 0.7)
    stats = {"elpd_loo": float(res.elpd), "loo_se": float(res.se), "p_loo": float(res.p),
             "n_pareto_bad": int((k > good_k).sum())}
    return elpd_i, stats


def build_model(formula: str, data: pd.DataFrame, prior_sigma: float = 1.0):
    import bambi as bmb
    return bmb.Model(formula=formula, data=data, family="cumulative",
                     priors={"common": bmb.Prior("Normal", mu=0, sigma=prior_sigma)})


def fit_cached(formula: str, data: pd.DataFrame, species: str, label: str,
               prior_sigma: float = 1.0, refit: bool = False,
               cores: Optional[int] = None, verbose: bool = True) -> Dict:
    """Fit (or load) one cumulative-logit model.

    Returns a dict with ``idata`` (posterior + sample_stats + observed_data),
    ``elpd_i`` (pointwise PSIS-LOO elpd), ``loo`` (summary), ``name`` and
    ``seconds`` (sampling time, NaN when loaded from the cache).
    """
    import arviz as az

    name = cache_name(formula, species, label, prior_sigma)
    nc, npz = MODEL_DIR / f"{name}.nc", ELPD_DIR / f"{name}.npz"
    tag = f"{SP_SHORT.get(species, species):<17} {label:<16} N(0, {prior_sigma:g})"

    if nc.exists() and npz.exists() and not refit:
        idata = az.from_netcdf(str(nc))
        z = np.load(npz)
        loo = {k: float(z[k]) for k in ("elpd_loo", "loo_se", "p_loo", "n_pareto_bad")}
        if verbose:
            print(f"{tag}  loaded from cache ({name})")
        return {"idata": idata, "elpd_i": np.asarray(z["elpd_i"], float), "loo": loo,
                "name": name, "seconds": np.nan}

    if verbose:
        print(f"{tag}  sampling ...", flush=True)
    t0 = time.time()
    model = build_model(formula, data, prior_sigma)
    idata = model.fit(draws=BAYES_DRAWS, tune=BAYES_TUNE, chains=BAYES_CHAINS,
                      cores=BAYES_CORES if cores is None else cores,
                      target_accept=BAYES_TARGET_ACCEPT, random_seed=BAYES_SEED,
                      idata_kwargs={"log_likelihood": True}, progressbar=False)
    secs = time.time() - t0
    elpd_i, loo = psis_loo(idata)
    np.savez_compressed(npz, elpd_i=elpd_i, **loo)
    slim = strip_heavy(idata)
    del idata
    slim.to_netcdf(str(nc))
    if verbose:
        print(f"{tag}  done in {secs / 60:.1f} min, elpd_loo = {loo['elpd_loo']:.1f}")
    return {"idata": slim, "elpd_i": elpd_i, "loo": loo, "name": name, "seconds": secs}


def paired_diff(elpd_a: np.ndarray, elpd_b: np.ndarray):
    """elpd(A) - elpd(B) and its paired SE, sqrt(n) * sd(elpd_i^A - elpd_i^B)."""
    d = np.asarray(elpd_a) - np.asarray(elpd_b)
    return float(d.sum()), float(np.sqrt(d.size) * d.std(ddof=1))


# ── Figures ──────────────────────────────────────────────────────────────────
def save_figure(fig, name: str) -> None:
    """Save a figure into figures/ as PDF."""
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    print(f"  saved figures/{name}.pdf")
