# RFR-LME

This repository contains code to produce RFR-LME indicators of ocean acidification in eleven U.S. large marine ecosystems (LMEs).

Instructions:
1. Download the latest SOCAT `.tsv` file and place it in the configured SOCAT directory.
2. Set `vrs` in `create_RFR_LME.m` to the matching dataset version (for example, `'SOCATv2026'`).
3. Update the directory locations in `rfr_lme_paths.m` if your storage layout differs. This includes the SOCAT input, `Data`, `Figures`, and `Models` directories, plus `dpath` for downloaded predictor datasets. The workflow initializes this path struct once and passes it into the functions that need it.
4. Run `startup.m`, then run `create_RFR_LME.m` from the repository.
