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
data/audio/                             one call clip per species (Fig. S5)
figures/                                every figure, produced by the notebooks (PDF)
tables/                                 every table (CSV)
tables/supporting/                      tables not in the paper
pyproject.toml, uv.lock                 pinned Python environment
requirements.txt                        for pip install environment
```

The notebooks in this repository were executed top to bottom; their saved outputs, `figures/` and
`tables/` are what that run produced.

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
data are archived at <http://hdl.handle.net/11336/301442>.

`data/audio/` holds the two clips used for Fig. S5, one isolated call per species (44.1 kHz, mono,
16-bit WAV, about 2 s each). Each file name gives the species, the recording date and start time, the
channel and the call's position in the original recording, e.g. `gastrotheca_20191013_200000_ch1_49.49s.wav`.

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

**Run time.** The paper fits 46 Bayesian models (main analysis and sensitivities) (4 chains × 1,000 draws after 1,000 tuning
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
| Fig. S5 (call spectrograms) | 01 §0.1 | `figures/FigS5_call_spectrogram_*` |

Table S5 (literature review) is not computed here. Every figure in `figures/` is in the paper.