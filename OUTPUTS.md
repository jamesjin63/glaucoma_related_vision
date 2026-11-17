# Outputs

All output paths below are relative to the configured output directory.

| Figure | File prefix in `figures/` | Content |
|---|---|---|
| 1 | Fig1_age_sex_dual_axis | Global 2021 age/sex numbers and age-specific rates |
| 2 | Fig2_SDI_age_specific_trends | Nine age-specific series within five SDI groups |
| 3 | Fig3_regional_sex_standardized_rates | Regional male/female truncated ASRs |
| 4 | Fig4_regional_SDI_trajectories | Regional SDI trajectories, 1990/2021 markers and separate cross-sections |
| 5 | Fig5_country_maps | A cases, B ASPR, C DALYs, D ASDR, E ASPR EAPC, F ASDR EAPC |
| 6 | Fig6_country_SDI_associations | 2021 country SDI associations |
| 7 | Fig7_DALY_reference_frontier_and_gaps | Historical DALY frontier and 2021 gaps; reversed vertical axes |
| 8 | Fig8_HFPG_standardized_trends_and_age_specific | HFPG truncated ASDR trends and 2021 age-specific rates |
| 9 | Fig9_age_specific_BAPC_projections | Nine age-specific projections for each outcome |
| S1 | FigS1_1990_2021_summary | Global/five-SDI endpoint summary |
| S2 | FigS2_global_projected_rates_and_counts | Global ASRs and projected totals by sex |
| S3 | FigS3_temporal_validation | Holdout predictions and constant-rate comparator |

Figures are exported as PNG and PDF. Figure 5 PDF contains an overview followed by six full-size panels. Some small locations lack polygons and remain in the country tables.

| Numerical output in `tables/` | Content |
|---|---|
| Table1_Prevalence_Global_SDI_21regions.csv / Table2_DALYs_Global_SDI_21regions.csv | Main historical estimates, endpoint changes and trends |
| eTable1_Country_Prevalence.csv / eTable2_Country_DALYs.csv | 612 country/sex records each, 26 columns |
| TableS4_HFPG_Age_specific_2021.csv | 54 records, 17 columns |
| TableS5_HFPG_Annual_1990_2021.csv | 192 records, 20 columns |
| TableS_BAPC_prior_and_period_sensitivity_2050.csv | Baseline, period RW1 and prior sensitivities |
| Table8_temporal_validation_fit1990_2011_predict2012_2021.csv | Holdout performance by outcome and sex |
| Table7_BAPC_global_annual_1990_2050.csv | Annual global fitted and projected ASRs/counts |
| TableS_BAPC_age_specific_1990_2050.csv | Age-specific fitted and projected results |
| Figure5_six_panel_source.csv | Country-level numerical values for all map panels |
| Figure4_annual_source.csv / Figure4_1990_2021_regional_correlations.csv | Trajectory points and cross-section statistics |
| Figure1_age_sex_source.csv | Numbers, age-specific rates and original uncertainty limits |

Other machine-readable tables provide annual rates, trend/HAC results, frontier curves/gaps, SDI effects, HFPG source values, standard weights, demographic scenario populations and scientific consistency diagnostics. Prepared historical age-specific data and SDI associations are exported to `data/` within the output directory. No manuscript text or publication workbook layout is generated.
