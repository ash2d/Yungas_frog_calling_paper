# Microclimate and calling activity of two Yungas frogs — analysis code

Code and data to reproduce every figure and table in

> Ash, D., Boullhesen, M., Bonnaffé, W. & Akmentins, M. S. *Non-linear microclimatic constraints on
> calling activity in two terrestrial breeding frogs of the Yungas Andean forests.*

The study models hourly calling intensity (amphibian calling index, 0–3) of *Gastrotheca chrysosticta*
and *Oreobates berdemenos*, recorded by passive acoustic monitoring in Baritú National Park, Argentina,
with Bayesian cumulative-logit (ordinal) regression on temperature, relative humidity, heterospecific
calling and cyclic (Fourier) time-of-day and day-of-year terms.

## Contents

```
01_main_analysis.ipynb                  main analysis: Figs 1, 2, S1, S2, S3, S5; Tables 1, 2, S1, S3; text statistics
02_fourier_harmonics_sensitivity.ipynb  Table S4, Fig. S4 (2–5 Fourier harmonics)
03_blocked_cross_validation.ipynb       Table S2 (6 h / 12 h blocked cross-validation)
frogcall.py                             shared paths, data preparation and cached model fitting
fit_all_models.py                       optional: fit all 46 Bayesian models in parallel
data/df_helechos_with2020.csv           hourly calling indices, temperature and relative humidity
data/audio/                             the ten recordings used for the spectrograms (Fig. S5)
figures/                                every figure, as produced by the notebooks (PDF)
tables/                                 every table (CSV; LaTeX where the notebooks write it)
tables/supporting/                      further tables written along the way
pyproject.toml, uv.lock                 pinned Python environment
requirements.txt                        the same pins, for pip
```

The notebooks in this repository were executed top to bottom; their saved outputs, `figures/` and
`tables/` are exactly what that run produced.

## Data

`data/df_helechos_with2020.csv` has one row per annotated hour (13,257 rows, 2 Sep 2018 – 18 Jul 2020):

| Column | Meaning |
|---|---|
| `DateTime` | hour of the 1-min recording (local clock time) |
| `Date`, `Month`, `Day`, `hour` | the same timestamp split up (`Day` = day of year) |
| `Temp` | air temperature, °C (HOBO MX2301A) |
| `RH%` | relative humidity, % (missing after 16 Dec 2019) |
| `Gastrotheca chrysosticta`, `Oreobates berdemenos` | calling index: 0 none, 1 one male, 2 several distinguishable males, 3 full chorus |

Models use the 11,302 hours that have both temperature and humidity. The annotations and microclimate
data are also archived at <http://hdl.handle.net/11336/301442>.

`data/audio/` holds ten 60-s recordings from the 2018 deployment (8 kHz, 16-bit stereo WAV): the five
candidate calls per species shown in the spectrogram section, including the two used for Fig. S5.

## Installation

With [uv](https://docs.astral.sh/uv/) (recommended; installs Python 3.12 and the exact locked versions):

```bash
git clone <this repository>
cd frogs_bayes_weather_paper
uv sync
```

Or with pip in a Python 3.12 virtual environment:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Main versions: PyMC 6.0.1, Bambi 0.18.0, ArviZ 1.1.0, statsmodels 0.14.5, NumPy 2.3.2, pandas 2.3.2.

## Running

Open the notebooks in Jupyter (`uv run jupyter lab`) or with your IDE and run them in order:

```bash
uv run jupyter nbconvert --to notebook --execute --inplace 01_main_analysis.ipynb
uv run jupyter nbconvert --to notebook --execute --inplace 02_fourier_harmonics_sensitivity.ipynb
uv run jupyter nbconvert --to notebook --execute --inplace 03_blocked_cross_validation.ipynb
```

**Run time.** The paper rests on 46 Bayesian model fits (4 chains × 1,000 draws after 1,000 tuning
steps, `target_accept = 0.9`, `random_seed = 42`), each taking about 2–6 min. Every fit is cached in
`outputs/` the first time it runs and loaded afterwards, and the cache is shared by the three notebooks
(`frogcall.fit_cached`). Run sequentially from the notebooks, a first run takes a few hours. To fill
the cache in parallel first (about 1 h with 4 jobs on an 8-core laptop; allow ~2 GB RAM per job):

```bash
uv run python fit_all_models.py --jobs 4
```

then run the notebooks as above (about 15 min in total once the cache exists). Delete `outputs/` to
refit everything.

Chains are seeded individually from `random_seed = 42`, so the draws do not depend on the number of
cores used. Set the environment variable `FROGCALL_CORES` to change how many chains run in parallel.

## Where each paper item comes from

| Paper item | Notebook | Output |
|---|---|---|
| Fig. 1 (seasonal and diel calling) | 01 §9 | `figures/Fig1_top_seasonal_calling_environment`, `figures/Fig1_bottom_diel_rose_calling` |
| Fig. 2 (response to temperature and humidity) | 01 §7 | `figures/Fig2_top_response_curves_temperature`, `figures/Fig2_bottom_response_curves_humidity` |
| Table 1 (model comparison, PSIS-LOO) | 01 §4.1 | `tables/Table1_model_comparison.csv` |
| Table 2 (posterior summaries) | 01 §4.2 | `tables/Table2_posterior_summaries.csv` |
| Results 3.1 (data summary) | 01 §2.1 | `tables/Results3.1_data_summary.csv` |
| Results 3.4 (temperature optima) | 01 §9.2 | `tables/Results3.4_temperature_optima.csv` |
| Other statistics in the text | 01 §10 | `tables/Text_statistics.csv` |
| Table S1 (convergence diagnostics) | 01 §5.6 | `tables/TableS1_convergence_diagnostics.csv` |
| Fig. S1 (coefficients across models) | 01 §5.1 | `figures/FigS1_top_coefficients_across_models`, `figures/FigS1_bottom_fourier_coefficients_across_models` |
| Table S2 (blocked cross-validation) | 03 | `tables/TableS2_blocked_cross_validation.csv` |
| Fig. S2 (posterior predictive checks) | 01 §6 | `figures/FigS2_posterior_predictive_checks` |
| Fig. S3 (temperature × humidity surface) | 01 §7 | `figures/FigS3_response_surface_temp_rh` |
| Table S3 (prior sensitivity) | 01 §5.5 | `tables/TableS3_prior_sensitivity.{csv,tex,md,html}` |
| Table S4 (Fourier harmonics) | 02 | `tables/TableS4_fourier_harmonics_hypothesis_ladder.{csv,tex}` |
| Fig. S4 (coefficients by harmonic count) | 02 | `figures/FigS4_coefficients_by_harmonics` |
| Fig. S5 (call spectrograms) | 01 §0.1 | `figures/FigS5_call_spectrogram_*`, `figures/FigS5_candidate_*` |

Table S5 (literature review) is not computed. Every figure in `figures/` is in the paper.

## Notes on reproduction

### TODO: confirm that fresh fits match the original fits

The saved outputs in this repository were produced from the **original project's model fits**:
the posterior draws and pointwise PSIS-LOO values of all 46 models were copied into `outputs/` and
loaded by the notebooks; no model was re-sampled for this run. For the two N(0, 0.5) fits, PSIS-LOO
was recomputed from the saved posteriors and matched the original values (−2557.75 and −2232.98).
All tables listed below match the manuscript.

Still to do: run the notebooks from an empty `outputs/` folder (or `python fit_all_models.py` first)
and confirm the fresh fits give the same results. So far:

- The four fresh fits finished before the full run was stopped (*G. chrysosticta* `model_0`–`model_3`,
  K = 2) gave elpd_LOO −3028.2, −2645.4, −2540.5, −2538.5. Table 1 gives −3028, −2645, −2540, −2538.
- The one fresh fit compared draw by draw (*G. chrysosticta* `model_1`, K = 2) was **not**
  bit-identical to the original fit; posterior means agreed only to Monte Carlo error (`temp_z`
  1.3729 vs 1.3737). The cause is not yet known. The original fits were made with the same PyMC
  (6.0.1), Bambi (0.18.0), PyTensor, NumPy and SciPy versions as this environment, and refitting
  with the original notebook gives identical draws every time. The only differences in the sampling
  call are `progressbar=False` and extra (unused) columns in the data frame.
- In this environment the draws are deterministic: the same seed gives identical draws with 1 or more cores.

### Checked against the manuscript (from the original fits)

| Item | Status |
|---|---|
| Table 1, Table 2, Table S3 (incl. ΔELPD −18.7 ± 3.5 / −8.7 ± 2.2), Table S4 | identical to the manuscript |
| Table S2 (1 h, 6 h, 12 h) | identical; the manuscript's 12 h column is the blocks starting at midnight (`12 h blocks`); all 300 refits converged |
| Table S1 | identical (R̂, bulk and tail ESS) |
| Temperature optima 18.3 °C [17.5, 19.3] and 19.8 °C [19.1, 20.5] | identical |
| Fig. S2 KL values and stratum sizes (769, 1,015, 3,819) | identical |
| Residual lag-1 r = 0.068 / 0.075; diel / seasonal PPC r = 0.90–0.98 | identical |
| Figs 1, 2, S1–S5 | regenerated |

### What differs from the original notebooks

The notebooks are assembled from the project's original notebooks (`Bayes_GLM_cum_ord_logit_hypoth`,
`fourier_harmonics_sensitivity`, `fourier_harmonics_hour_vs_day_loo`, `blockcv_residual_diagnostics`,
and `reports/manuscript_checks.py`). The analysis code is unchanged except for the following:

- **Paths and caching.** Paths are relative to the repository. All Bayesian fitting goes through
  `frogcall.fit_cached`, which samples with the original settings and caches the result.
- **Table cells.** New cells produce Table 1, Table 2, Table S1, Table S2, the Table S4 CSV and the
  text statistics.
- **Removed.** Sections that produce nothing in the paper: ALE plots, Bayes factors, the prior checker, K = 1 and K = 6 fits, frequentist AIC/BIC plots and the AIC/BIC scan over K,
  and every figure or table not in the paper (prior-sensitivity plots, the 2 × 2 PPC composite, the
  temperature-optimum density, the ladder and block-length plots, the ΔELPD autocorrelation plot, the
  frequentist coefficient table and the alternative Fourier-ladder table).
