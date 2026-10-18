# Bio-economic simulation of site-specific weed management under imperfect sensing

Code accompanying:

> Economic evaluation of site-specific weed management under imperfect sensing: a multi-year bio-economic simulation study.

This repository contains the four scripts used to produce every result in the paper: grid construction and weed-count aggregation, the bio-economic simulation, the analysis of the simulation output, and the marginal shadow cost figure.

## Contents

| File | Supplementary Code | Language | Purpose |
|---|---|---|---|
| `01_grid_creation_weed_counts.py` | 3 | Python | Builds a row-aligned grid over a field boundary and aggregates UAV weed detections into per-cell counts. |
| `02_bioeconomic_simulation.py` | 1 | Python (QGIS) | Runs the bio-economic simulation for one field, one resolution and one decision strategy, including the marginal shadow costs. |
| `03_main_results_analysis.R` | 2 | R | Produces all tables and figures from the simulation output. |
| `04_marginal_shadow_costs.R` | 4 | R | Produces the marginal shadow cost figure from the per-cell output of the simulation. |

## Pipeline

```
field boundary + tramline + weed detection points
        |
        |  01_grid_creation_weed_counts.py      (once per field per resolution)
        v
   grid_<size>m_weed_counts.gpkg
        |
        |  02_bioeconomic_simulation.py         (once per field x resolution x strategy)
        v
   results_thesis.csv                    --+
   results_thesis_recall_curve.csv         +--> 03_main_results_analysis.R ---> tables + figures
   results_thesis_dd_kill_sweep.csv      --+
   results_thesis_yearly.csv                    (yearly costs; for other discount rates or horizons)
   cells/<project>_SCENARIO_1_ECON.csv   -----> 04_marginal_shadow_costs.R  ---> shadow cost figure
```

The published results cover 19 fields, 3 spraying resolutions (1 x 1 m, 50 x 50 cm, 25 x 25 cm) and 4 decision strategies, i.e. 228 runs of `02_bioeconomic_simulation.py`.

## Requirements

**Python (script 01):** Python >= 3.8 with `geopandas >= 0.12`, `numpy`, `pandas`, `shapely`.

**Python (script 02):** QGIS 3.x with its bundled Python (tested with QGIS 3.36). The script uses the `qgis.core` API and is run from the QGIS Python console. No additional packages are required.

**R (scripts 03 and 04):** R >= 4.2 (tested with R 4.3.1) with `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `multcomp`, `multcompView` and `pbkrtest`. A cairo-enabled R build is recommended so that the euro sign renders correctly in the PDF output.

## Running the pipeline

Every script has a configuration section at the top. No paths are hard-coded elsewhere. By default, all scripts read and write in `C:/SSWM_simulation/results`.

**1. Build the grids.** Set the input paths and `CELL_SIZE` in `01_grid_creation_weed_counts.py` and run it once per field and resolution. Load the resulting GeoPackage into a QGIS project; the simulation expects the grid layer to be named `grid_rotated` with the detected weed counts in the fields `DicotCount` and `MonoCount`.

The **QGIS project file name** identifies the field and the resolution in all outputs. Name the projects `<field>` for 1 x 1 m, `<field>-50` for 50 x 50 cm and `<field>-25` for 25 x 25 cm (e.g. `R101_econ_model`, `R101_econ_model-50`, `R101_econ_model-25`); `03_main_results_analysis.R` reads the resolution from this suffix.

**2. Run the simulation.** Open the QGIS project for one field at one resolution, set `RESULTS_CSV_PATH` and `SCENARIO_ID` in `02_bioeconomic_simulation.py`, and run it from the QGIS Python console. Repeat for each field, resolution and strategy.

`SCENARIO_ID` selects the decision strategy: `1` economic threshold (ET), `2` presence-based (PB), `3` dual-dose (DD), `4` forward-looking economic threshold.

Per-cell results are written to the grid layer. For the ET strategy at 1 x 1 m, they are also written to `cells/<project>_SCENARIO_1_ECON.csv` next to the results file (`ENABLE_CELL_EXPORT`, `CELL_EXPORT_SCENARIOS`, `CELL_EXPORT_CELL_SIZE_M`); these files are the input of script 04.

The four `results_thesis*.csv` files are opened in **append** mode, so results accumulate across runs. If you re-run a strategy, archive or delete the existing CSVs first, or remove the superseded rows afterwards. Every row of the main file records the configuration that produced it, so mixed files can also be filtered after the fact. The per-cell files are overwritten.

**3. Analyse.** Set `csv_folder` in `03_main_results_analysis.R` to the folder holding the CSVs and source the script. It writes all tables and a multi-page PDF of the figures.

Before producing any output the script checks that all runs share the same model configuration, that every field x resolution x strategy combination appears exactly once, and that the dual-dose rows and the recall curve reproduce the baseline. It stops with a diagnostic message if a check fails. These checks exist because the results are assembled from many separate runs, where a missing or duplicated row would otherwise silently change a mean rather than raise an error.

**4. Marginal shadow costs.** Set `INPUT_FOLDER` in `04_marginal_shadow_costs.R` to the `cells` folder written in step 2 and source the script. It reads the files of the ET strategy (`FILE_PATTERN`); the published figure uses the ET strategy at 1 x 1 m for all 19 fields.

## Reproducing the published results

The default settings in `02_bioeconomic_simulation.py` are those used for the reported results. The model retains alternative settings for validation; they are described in the script header and are not used for any reported result. Two of them are worth stating explicitly:

- `DD_FUTURE_COST_MODE = 'dose_scaled'` charges herbicide in proportion to the applied dose in every simulated year, so a dual-dose background cell is charged 50% of the full rate. The alternative, `'v13'`, charges the full cost in background-dose years.
- `DD_BASE_KILL = 0.85` is the kill fraction of the reduced background dose; full-dose cells always use `full_kill_rate = 0.95`.

Setting `DD_FUTURE_COST_MODE = 'v13'` together with `DD_BASE_KILL = 0.95` makes the dual-dose strategy algebraically identical to broadcast application, so its differential is zero everywhere. This is retained as a validation identity and must not be used to produce reported results; `03_main_results_analysis.R` rejects input produced with `'v13'`.

Similarly, `FUTURE_DETECTION = 'probabilistic'` is the reported behaviour. The `'expected'` setting, in which every cell containing emerged weeds is detected from year 1 onward, is retained only for comparison and removes the compounding of repeated detection failures that drives the resolution effect reported in the paper.

The model is deterministic: future-year detection is a probabilistic process evaluated analytically over weighted seedbank trajectories rather than by sampling, so repeated runs of the same configuration give identical results.

## Data availability

The UAV weed map data used in the analysis are available from the corresponding author on reasonable request.

## License

MIT License, see `LICENSE`.
