# Input data

Place inputs in `data/local/`, or supply an external directory with `--data-dir`. Raw GBD exports, support snapshots and fitted model objects are not included in the repository or source archive. Obtain and use them under their applicable terms. The local working copy can contain these files without making them part of the source distribution.

```text
inputs/
  cause/*.csv
  hfpg/metadata/hfpg_download_manifest.json
  hfpg/<paths listed in the manifest>
  support/standard_weights_source.csv
  support/SDI_support_1990_2021.csv
  support/population_support_global_2011_2050.csv
  support/countries.csv
  support/world_map_support.geojson
```

## GBD disease exports

Obtain GBD 2021 exports from the IHME GBD Results Tool, https://vizhub.healthdata.org/gbd-results/. This code uses cause 670 (glaucoma-related vision loss), measures 5 (Prevalence) and 2 (DALYs), metrics 1 (Number) and 3 (Rate per 100,000), sexes 1 (Male), 2 (Female) and 3 (Both), years 1990–2021, and ages 16–20, 30–32, 235 and 228. The nine individual age IDs correspond to 55–59, 60–64, 65–69, 70–74, 75–79, 80–84, 85–89, 90–94 and 95+; 228 is the 55+ aggregate.

Place the CSV exports directly in `cause/`. The existing analysis export layout comprises eight CSVs, including aggregate HFPG files. Filenames are not used as identifiers. Required columns are `measure_id`, `location_id`, `location_name`, `sex_id`, `sex_name`, `age_id`, `cause_id`, `metric_id`, `year`, `val`, `lower`, `upper`; aggregate risk exports additionally have `rei_id`. Other descriptive columns may remain.

Disease data must cover 231 locations: global (1), five SDI groups (44635, 44634, 44639, 44636, 44637), 21 GBD regions and 204 countries/territories. Regional IDs are 5, 9, 21, 32, 42, 56, 65, 70, 73, 96, 100, 104, 120, 124, 134, 138, 159, 167, 174, 192, 199. The code identifies the country set from an aggregate disease export containing exactly the global/five-SDI/204-country locations (210 locations) and only age 228. Preserve that export grouping.

The disease key is `(measure_id, location_id, sex_id, age_id, cause_id, metric_id, year)`. Equal overlaps between exports are accepted; conflicting values fail. The aggregate risk CSVs use age 228, cause 670, risk 105, measure 2, and metrics 1/3. They provide the original aggregate count UIs and crude rates, and must cover global/five-SDI/21-region locations for all three sexes and years.

`python3 run.py --stage historical` converts these raw CSVs directly to the prepared historical age-specific input in the output directory. No manually prepared Parquet or hidden conversion step is required.

## Age-specific HFPG exports

Use GBD 2021 risk export version 8016, cause 670, risk 105, measure 2, metrics 1/3, the nine age IDs above, all three sexes, all 231 locations and years 1990–2021. Each annual CSV contains 12,474 rows (399,168 total). Unlike disease files, the source column names are `location`, `sex`, `age`, `cause`, `rei`, `measure`, `metric`, `year`, `val`, `lower`, `upper`. The code explicitly maps these ID fields to `*_id`.

The unique key is `(location, sex, age, metric, year)` under the fixed cause/risk/measure selection. Preserve original lower limits, including negative values.

The manifest identifies the version and each annual CSV relative to `hfpg/`:

```json
{
  "study": "GBD 2021",
  "version": 8016,
  "status": "complete",
  "files": [
    {"year": 1990, "path": "annual/1990.csv", "sha256": "SHA256_OF_THE_ACTUAL_FILE"}
  ]
}
```

The `files` array must contain all 32 years, not just the illustrated record. Compute SHA256 over each original file's bytes with `shasum -a 256`; the script verifies these hashes before reading the annual exports. Manifest paths must remain within `hfpg/`.

## Analytical support inputs

| File | Required fields and key | Coverage / units |
|---|---|---|
| standard_weights_source.csv | `age_id`, `Population_weight`; unique age ID | Nine selected ages; weights normalized over 55+ by the analysis |
| SDI_support_1990_2021.csv | `location_id`, `year`, `sdi`; unique location/year | Global, regions and countries, 1990–2021; SDI dimensionless |
| population_support_global_2011_2050.csv | `age_id`, `year`, `sex_name`, `val`; unique age/year/sex | Global, nine ages, Male/Female/Both, 2011–2050; persons |
| countries.csv | `location_id`, `location_name`; unique location ID | The 204 GBD countries and territories |
| world_map_support.geojson | GeoJSON with `location_id` property and valid geometries | Country polygons matching GBD location IDs |

The local weights, SDI, population and map support objects were exported from existing analytical support files. Their complete upstream release metadata, acquisition dates and redistribution permissions were not established. No authenticated dataset/version or forecast scenario is assigned to the future-population snapshot. Obtain those source records from the data provider before interpreting it as a specific published projection release. A similarly shaped substitute is not guaranteed to reproduce these results.

The actual population model uses only the global support file. Historical denominators come from the GBD disease Number/Rate pairs; future values are rescaled to the 2021 denominator for each age/sex. The unused SDI population object is not a required input. A data hash establishes byte identity, not dataset provenance.
