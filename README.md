# Glaucoma-related vision loss

Analysis code for GBD 2021 estimates of glaucoma-related vision loss among adults aged 55 years and older, 1990–2021. Historical analyses cover the global population, five SDI groups, 21 regions and 204 countries and territories. Conditional global projections extend to 2050.

The analysis includes truncated age-standardized rates, EAPCs with ordinary and HAC confidence intervals, SDI associations, country-cluster bootstrap frontiers, HFPG-attributable burden, BAPC forecasts, temporal holdout validation and prior/period sensitivity analyses.

## Requirements

Python 3.11 with the packages in `requirements.txt`; R with BAPC, INLA and data.table. Maps additionally require sf, ggplot2, patchwork, scales, Cairo support and the `pdfunite` executable from Poppler.

The recorded model environment used R 4.2.2, BAPC 0.0.36, INLA 23.04.24 and data.table 1.17.0. The Python dependency versions reflect the recorded supplementary-analysis environment. BAPC and INLA require their respective distribution sources; they are not installed by pip. Map generation was exercised with sf 1.0-16, ggplot2 4.0.3, patchwork 1.3.2 and scales 1.4.0.

```bash
python3 -m pip install -r requirements.txt
python3 run.py --data-dir /path/to/inputs --output-dir /path/to/results
```

Inputs are not distributed with the source. Prepare the directory described in [data/README.md](data/README.md). Local defaults are `data/local/` for inputs and `results/` for outputs. Both are excluded from version control. Input and output directories must be separate, with neither containing the other.

The default command calculates historical results, fits BAPC models, calculates HFPG results, exports supplementary tables and draws the final figures. Fitting BAPC can take substantial time. Individual modules can be run with `--stage historical`, `bapc`, `hfpg`, `tables` or `figures`; run historical calculations before the dependent modules and BAPC before projection figures.

```bash
python3 run.py --stage historical
python3 run.py --stage bapc
python3 run.py --stage hfpg
python3 run.py --stage tables
python3 run.py --stage figures
```

The R model script saves fitted objects under the output directory and reuses existing model files there. Use a new output directory for changed inputs or settings. `--model-cache-dir /path/to/models` copies existing RDS models into the output directory for reuse; these must match the input data and model configuration. Reusing models regenerates their summaries and does not constitute a new fit.

See [METHODS.md](METHODS.md) for implementation assumptions and [OUTPUTS.md](OUTPUTS.md) for figure and table mapping. The code does not download data, assemble manuscripts or submit files.
