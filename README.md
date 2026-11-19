# RFR-LME

This repository contains MATLAB code to create regional ocean-acidification estimates for eleven U.S. Large Marine Ecosystems (LMEs). It uses surface ocean CO<sub>2</sub> observations from the Surface Ocean CO<sub>2</sub> Atlas ([SOCAT](https://socat.info)) with environmental predictors, segments each LME into overlapping objective clusters, fits random-forest regression (RFR) models, and produces gridded indicators and regional statistics.

## Workflow

The main workflow is `create_RFR_LME.m`. It:

1. Loads and filters the SOCAT observations according to the user-provided version tag.
2. Aggregates observations into monthly 0.25-degree latitude/longitude cells.
3. Imports and aligns environmental predictors, including salinity, sea-surface temperature, ice, chlorophyll, wind, pressure, atmospheric CO2, and bathymetry. The default data sources for environmental predictors are set in `create_RFR_LME.m`.
4. Extracts data for each LME, prepares model inputs, clusters the data, and trains regional RFR models.
5. Predicts *p*CO<sub>2</sub>, estimates total alkalinity using [ESPER](https://github.com/BRCScienceProducts/ESPER), and calculates other ocean-acidification indicators using [CO2SYS](https://github.com/jonathansharp/CO2-System-Extd), then calculates error and regional/global statistics.
6. Saves figures, model files, tables, and gridded NetCDF products.

## Setup

Run the workflow in MATLAB. Downloaded predictor data may require network access and credentials for the selected providers.

1. Download the latest [SOCAT](https://socat.info) `.tsv` file and place it in the directory configured by `paths.socat`. The file name must match the version tag, such as `SOCATv2026.tsv`.
2. Set `vrs` near the top of `create_RFR_LME.m` to that version tag.
3. Set the input and output directories in `rfr_lme_paths.m`:
   - `paths.data`, `paths.figures`, and `paths.models` are the workflow's data, figure, and model directories.
   - `paths.socat` is the directory containing the SOCAT `.tsv` file.
   - `paths.dpath` is the downloaded-predictor data directory. It is also where the workflow expects the 60-arcsecond [ETOPO 2022](https://doi.org/10.5194/essd-17-1835-2025) bathymetry file.
4. Configure provider access in `create_RFR_LME.m` as needed: the Copernicus Marine toolbox path and credentials, NASA Earthdata token, and Climate Data Store token. Do not commit personal credentials or tokens.
5. In MATLAB, set the Current Folder to the repository root, then run:

   ```matlab
   startup
   create_RFR_LME
   ```

`startup.m` adds the project code to the MATLAB path and removes `old_code` from the path. The main workflow initializes the `paths` structure once and passes it to functions that read or write the configured directories.

## Outputs and runtime

- Intermediate data and NetCDF products are written under `paths.data`.
- Models are written under `paths.models`.
- Figures are written under `paths.figures`.
- Regional and global summary tables are written under `IndsAndStats/` in the current working directory.

The end of `create_RFR_LME.m` also retrains models for cluster counts from 1 to 10. This sensitivity run can take additional time.
