# Statistical implementation

## Historical estimates

Disease inputs are GBD 2021 cause 670, measures 5 (prevalence) and 2 (DALYs), metrics 1 (Number) and 3 (Rate), sexes 1–3, years 1990–2021. Nine age groups cover 55–59 through 90–94 years and 95+; age 228 is the aggregate 55+ export. Overlapping exports are deduplicated only after confirming equal point estimates and bounds for each key. Conflicting overlaps stop execution.

Population is inferred independently from disease Number / Rate × 100,000 and checked between outcomes. Nine age-group populations are summed for the crude denominator. Original aggregate number and crude-rate uncertainty limits come from age 228; age-group marginal limits are not summed to reconstruct an aggregate UI.

Truncated ASRs are sums of the nine age-specific rates multiplied by fixed standard-population weights normalized over ages 55+. Weighted lower and upper age-specific limits are approximate marginal bounds, without verified joint 95% coverage. Posterior GBD draws are unavailable.

EAPC is 100 × (exp(beta) − 1), from log(rate) regressed on year using all 32 annual observations. Ordinary regression confidence intervals and a three-lag HAC sensitivity analysis are exported separately. SDI Spearman correlations and log-linear percentage associations use independent cross-sections in 1990 and 2021. Historical trajectories do not make repeated annual observations independent.

Frontiers use 1,000 whole-country bootstrap samples and a non-increasing empirical reference curve on a fixed SDI grid (seed 20261113). Exported gaps are non-negative descriptive rate differences and do not identify avoidable burden or healthcare efficiency.

## BAPC

Global outcomes are modeled separately by sex. Both prevalence and DALY responses are rounded with R `round()` before being passed to BAPC's count model. Applying a Poisson count model to rounded DALYs is an approximation, not a literal model of integer DALY events. Age groups are five years wide and periods are annual (`gf=5`); the terminal 95+ group is an approximation.

Age, period and cohort components use RW2 in the baseline model, each with `prior="loggamma", param=c(1, 0.00005)`. An iid overdispersion term uses `param=c(1, 0.005)`. These are the parameter vectors supplied to BAPC/INLA for log-precision priors. Other component options follow the recorded BAPC 0.0.36 defaults; they are not asserted to be explicitly scaled or constrained beyond the implementation in that version.

Sensitivity scenarios replace the period RW2 with RW1 while retaining age/cohort RW2, or change the age/period/cohort second prior parameter to 0.000025 or 0.0001. Overdispersion parameters remain unchanged. Both-sex 2050 results are compared with the baseline. The labels `half_prior_scale` and `double_prior_scale` refer to that supplied second parameter, not a demonstrated proportional change in posterior variance.

Temporal validation fits 1990–2011 and predicts 2012–2021 using the known holdout-period populations. The comparator holds 2011 age-specific rates constant. MAPE, bias and interval coverage are exported. Prediction intervals for totals are sums of marginal age-specific bounds and are approximate. Age-standardized intervals and age-specific intervals are taken from BAPC outputs; no GBD input draws or future population uncertainty are propagated.

Historical populations are inferred from disease exports. Future population is a local support scenario, scaled separately for each sex and age to the 2021 inferred population. Its upstream dataset release and forecast scenario have not been authenticated. These are conditional projections. The inconsistent SDI population support object is not used in fitting or distributed as a required input.

## HFPG

GBD 2021 cause 670, risk 105, measure 2 exports provide nine age groups, 231 locations, three sexes and Number/Rate for each year. Original negative uncertainty limits are retained. Attributable fractions divide attributable DALYs by matched total glaucoma DALYs; the aggregate fraction is a ratio of counts, not an average of age-specific fractions. HFPG-attributable results do not establish that glucose-lowering treatment prevents glaucoma.
