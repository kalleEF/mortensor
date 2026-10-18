# Methods and sources

mortensor implements published methods. This page lists, for every
method of the package, its source and how the implementation relates to
it:

- **published**: the method as published (numerical safeguards aside);
- **published, deviations**: the published method with stated changes;
- **package-specific**: no published source is known to us; the help
  page says so. Package-specific methods are never the default unless
  stated here, and deprecated ones warn when used.

The help page of each function gives the details and the full
references. Corrections are welcome: a published source for a method
listed as package-specific, or a deviation that is not listed.

## Data and harmonisation

| Function, option | Method | Source | Status |
|----|----|----|----|
| [`mt_as_exposure()`](../reference/mt_as_exposure.md) | person-years from population counts: the mid-year population, or the mean of the populations at the start and the end of the year (1 January or 31 December) | Preston, Heuveline and Guillot (2001) | published |
| [`mt_data()`](../reference/mt_data.md) | deaths of unknown age distributed over the age groups of the same stratum and period in proportion to the known deaths | practice of the Human Mortality Database | published |
| `mt_fill_suppressed(hierarchy = )` | iterative proportional fitting of the missing cells to published totals, from the expected deaths at the pooled rate | Deming and Stephan (1940); Bishop, Fienberg and Holland (1975) | published, deviations: bounds of the suppressed counts enforced by clipping within each step (package-specific; no convergence guarantee once a bound binds) |
| `mt_fill_suppressed(lower = , upper = )` | cells that no published total covers: mean of the Poisson distribution with the expected deaths at the pooled rate, truncated to the bounds | none | package-specific |
| `mt_ungroup(method = "pclm")` (default) | penalised composite link model with the exposures as offset, 1D and 2D; the estimates are the fitted values (`preserve_totals = FALSE`) | Eilers (2007); Rizzi, Gampe and Eilers (2015); Rizzi et al. (2019); defaults of the ungroup package | published, deviations: in 2D both smoothing parameters are chosen by BIC (ungroup fixes them); an open age interval of the grid is split with stationary-population weights (package-specific, default) |
| `mt_ungroup(preserve_totals = TRUE)` | estimates rescaled to the observed group totals after the fit | none | package-specific, opt-in |
| `mt_ungroup(method = "monotone")`, `exposure_method = "monotone"` | monotone cubic interpolation of cumulative counts | Hyman (1983) | published |
| `mt_ungroup(method = "hyman_kannisto", "natural_gompertz", "gam")` | spline heuristics with parametric tails | none | package-specific, deprecated |
| `mt_disaggregate(method = "shape")` (default), [`mt_shape_rates()`](../reference/mt_shape_rates.md) | forecast draws of age groups split by a fixed age shape within each group, from the rates of the last observed years | none | package-specific (default) |
| `mt_disaggregate(method = "hyman_kannisto", "natural_gompertz")`, `mt_shape_rates(method = "gam")` | splines through the group rates of each draw; age shapes from an extrapolated GAM surface | none | package-specific, deprecated |
| `mt_smooth(method = "pspline")` | Poisson P-splines over age, or age and time, smoothing parameters by BIC | Currie, Durban and Eilers (2004) | published |
| `mt_smooth(method = "whittaker")` | weighted Whittaker smoother of the log rates with the deaths as weights, over age or time, smoothing parameter by GCV for every fit | Whittaker (1923); Eilers (2003) | published |
| `mt_smooth(method = "gam")` | quasi-Poisson GAM with the exposures as offset, `s(age)` or `te(age, year)`, REML | Wood (2017) | published |
| [`mt_smooth()`](../reference/mt_smooth.md) conventions | the first year of life is left out by default; an open interval is placed at its lower bound plus half the width of the interval below it (Whittaker smoother over age, GAMs) | none | package-specific conventions |
| [`mt_smooth_check()`](../reference/mt_smooth_check.md) | Pearson estimate of the dispersion of the raw around the smoothed deaths | McCullagh and Nelder (1989) | published |
| [`step_shocks()`](../reference/mt_prep.md) | shock years removed or interpolated log-linearly before fitting | none | package-specific (interpolated years enter as data) |

## Models

| Function, option | Method | Source | Status |
|----|----|----|----|
| [`lee_carter()`](../reference/lee_carter.md) | Poisson log-bilinear model (default), or SVD with `k` re-estimated to the total deaths (`adjust = "dt"`), the deaths by age (`"dxt"`) or life expectancy (`"e0"`) | Lee and Carter (1992); Brouhns, Denuit and Vermunt (2002); Booth, Maindonald and Smith (2002); Lee and Miller (2001) | published; age groups without deaths get a constant rate, and the ridge and SVD fallbacks for sparse counts are package-specific safeguards that warn |
| `lee_carter(choose_period = TRUE)` | fitting period by the ratio of the mean deviances of a linear and of the Lee-Carter index | Booth, Maindonald and Smith (2002) | published |
| `lee_carter(jump_off = "observed")`, `li_lee(jump_off = "observed")` | forecasts start from the observed rates of the last year | Lee and Miller (2001) | published |
| `lee_carter(rotate = TRUE)` | rotation of the age pattern of decline as life expectancy rises | Li, Lee and Gerland (2013); Sevcikova et al. (2016) | published, deviations: the path of life expectancy is that of the unrotated draws (MortCast takes it from a separate model) |
| `lee_carter(bootstrap = TRUE)` | parameters re-estimated on Poisson resamples of the fitted deaths | Renshaw and Haberman (2008) | published |
| `bootstrapped(type = "grid")` (default) | Poisson resamples of the fitted deaths on the grid, any model refitted | Renshaw and Haberman (2008) | published; models without fitted rates resample the deaths on the grid |
| `bootstrapped(type = "data")` | Poisson resamples of the observed deaths | Brouhns, Denuit and Van Keilegom (2005) | published, deviations: refitting the preprocessing in each replicate is package-specific |
| [`rwd()`](../reference/rwd.md) | random walk with drift on the log rate of every age; innovations resampled as age vectors of observed changes (default) or normal | standard benchmark, no reference given | published |
| `rwd(observation_noise = TRUE)` | moment estimator of a local level model with drift | none (the published counterpart is the model fitted by maximum likelihood; Harvey 1989) | package-specific, deprecated: use `arima_series(order = c(0, 1, 1), drift = TRUE)` |
| [`arima_series()`](../reference/arima_series.md), [`ets_series()`](../reference/arima_series.md) | ARIMA or exponential smoothing model per age, chosen automatically or given by `order` | Hyndman and Khandakar (2008) | published |
| `arima_series(order = c(0, 1, 1), drift = TRUE)` | local level model with drift in its reduced form, maximum likelihood | Harvey (1989) | published |
| [`hyndman_ullah()`](../reference/hyndman_ullah.md) | functional demographic model | Hyndman and Ullah (2007), via demography | published; smoothing without the monotonicity constraint where the monotonic spline cannot be fitted is a package-specific safeguard |
| [`coherent_fdm()`](../reference/coherent_fdm.md) | product-ratio method with functional data models | Hyndman, Booth and Yasmeen (2013), via demography | published |
| [`product_ratio()`](../reference/product_ratio.md) | product-ratio method with Lee-Carter components | none | package-specific, deprecated: use [`coherent_fdm()`](../reference/coherent_fdm.md) |
| [`li_lee()`](../reference/li_lee.md) | augmented common factor model | Li and Lee (2005); Lee and Carter (1992) | published, deviations: AR(1) processes without intercept; their coefficient bounded at 0.99 (package-specific safeguard, warns) |
| [`coda()`](../reference/coda.md) | compositional Lee-Carter model of life-table deaths, and its coherent version | Oeppen (2008); Bergeron-Boucher et al. (2017); MortalityForecast | published, deviations: the log rate of the open age interval is regressed on the indexes |
| [`coda_causes()`](../reference/coda_causes.md) | multiple-decrement CoDa model for causes of death | Oeppen (2008), abstract only | package reading of a published idea; experimental |
| [`tucker3()`](../reference/tucker3.md) | three-way Lee-Carter model: Tucker3 decomposition of the centred log rates | Russolillo, Giordano and Haberman (2011); De Lathauwer, De Moor and Vandewalle (2000) | published, deviations: HOOI estimation, the forecasts of the period factors and the ranks are the package’s; experimental |
| `tucker3(residual_noise = TRUE)`, `tucker3(transform = "coda")` | in-sample residual variation added to the draws; Tucker3 decomposition of transformed life-table deaths | none | package-specific: `residual_noise` opt-in (counts the observation noise twice), `"coda"` deprecated: use [`coda()`](../reference/coda.md) |
| [`cause_shares()`](../reference/cause_shares.md) | Tucker model of log-ratio coordinates of cause shares | none | package-specific, deprecated: use [`reconciled()`](../reference/reconciled.md) or [`coda_causes()`](../reference/coda_causes.md) |
| [`stmomo()`](../reference/stmomo.md) | GAPC models LC, CBD, APC, RH, M6, M7, M8 | Villegas, Kaishev and Millossovich (2018); Cairns et al. (2009) | published (StMoMo); `"rh"` is StMoMo’s version with a cohort effect without age modulation; initial exposures of the logit-link models approximated as central exposures plus half the deaths (package-specific) |
| [`stanmomo()`](../reference/stanmomo.md) | Bayesian GAPC models with a Poisson or negative binomial likelihood | Barigou et al. (2023), via StanMoMo; Delwarde, Denuit and Partrat (2007) | published |
| [`apc_bayes()`](../reference/apc_bayes.md) | Bayesian APC model with random-walk priors, fitted where the deaths were observed | Schmid and Held (2007), via bamp; Knorr-Held and Rainer (2001) | published; the mapping of draws of age groups to the grid by shape factors is package-specific |
| [`apc_select()`](../reference/apc_select.md) | choice of the random-walk orders by DIC | Spiegelhalter et al. (2002) | published |
| zero deaths | 0.5 deaths in cells without deaths in models of log rates | none | package-specific convention |
| [`mt_order_strata()`](../reference/mt_order_strata.md) | keeps an ordering between strata path by path | none | package-specific |
| predictive draws (`mt_forecast(predictive = TRUE)`, backtests) | Poisson deaths given the model’s rates; posterior negative binomial for `stanmomo(family = "nb")` | Brouhns, Denuit and Vermunt (2002) | published |

## Combining models and coherence

| Function, option | Method | Source | Status |
|----|----|----|----|
| [`ensemble()`](../reference/ensemble.md) | mixture (linear pool) of predictive distributions, equal weights by default | Gneiting and Ranjan (2013) | published |
| `ensemble(weights = "crps")` | weights inverse to the mean CRPS of an inner backtest | analogous to Bates and Granger (1969) | package-specific |
| `ensemble(weights = "stack")` | stacking of point forecasts, used as mixture weights | Breiman (1996); Kessy et al. (2022) | published, deviations: the weights optimise the medians and are used for the mixture |
| `ensemble(top_k = )` | trimming to the components with the lowest inner CRPS | none | package-specific |
| [`stratified()`](../reference/stratified.md) | different models for groups of strata | none | package-specific plumbing |
| [`reconciled()`](../reference/reconciled.md), [`mt_reconcile()`](../reference/mt_reconcile.md): bottom-up, OLS, structural, WLS, MinT | projection of the draws onto coherent forecasts | Hyndman et al. (2011); Athanasopoulos et al. (2017); Hyndman, Lee and Wang (2016); Wickramasuriya, Athanasopoulos and Hyndman (2019); Schäfer and Strimmer (2005) | published, deviations: WLS and MinT weights from the variances and covariances of the base draws (the papers use in-sample residuals) |
| non-negativity | negatives held at zero, the rest projected again | related to Wickramasuriya, Turlach and Hyndman (2020) | package-specific heuristic |
| `fixed =` | nodes kept unchanged (immutable forecasts) | Zhang et al. (2023) | published, approximated by a near-zero variance |
| `method = "proportional"` | top-down by forecast proportions | Gross and Sohl (1990); Athanasopoulos, Ahmed and Hyndman (2009) | published, deviations: one level, the proportions of the bottom forecasts below each fixed node |
| `method = "energy_score"` | score-optimal reconciliation | Panagiotelis et al. (2023); ProbReco | published, deviations: estimated per age; observation noise in the training forecasts; draws of each step taken without replacement; rescaling before the optimisation |
| `dependence = "empirical"` | reordering of separately simulated strata | related to the Schaake shuffle (Clark et al. 2004) | package-specific, opt-in (default `"none"`) |
| [`recalibrated()`](../reference/recalibrated.md) | spread of the draws rescaled by one factor for all horizons that minimises the CRPS of inner rolling-origin forecasts | member-by-member post-processing (Van Schaeybroeck and Vannitsem 2015); minimum-CRPS estimation (Gneiting et al. 2005) | published, deviations: the draws are rescaled around their median on the log scale, path by path |
| `recalibrated(scale = "horizon", bias = TRUE, ages = , quantity = )` | power law of the factor over the horizon, drift correction, age bands, life expectancy and ASR as targets | none | package-specific, opt-in |
| `recalibrated(target = "interval")` | interval score of one central interval | Gneiting and Raftery (2007) | published |
| [`anchored()`](../reference/anchored.md) | forecasts start from the observed rates of the last training year, with the offset kept in every forecast year | Lee and Miller (2001) | published (default) |
| `anchored(years = , decay = , shrink = TRUE)` | offset of several years, linear fade, empirical Bayes shrinkage of the offsets | DerSimonian and Laird (1986) for the variance estimator of the shrinkage | package-specific, opt-in |
| [`mt_jump_off()`](../reference/mt_jump_off.md) | step between the last observed or fitted year and the first forecast year | none | package-specific diagnostic |
| [`tuned()`](../reference/mt_tune.md), [`mt_tune()`](../reference/mt_tune.md) | choice of settings by rolling-origin or hold-out accuracy; one-standard-error rule | Tashman (2000); Hastie, Tibshirani and Friedman (2009) | published, deviations: the standard error of the one-SE rule is approximate; `strategy = "two_stage"` and the minimum number of training years are package heuristics |

## Evaluation

| Function, option | Method | Source | Status |
|----|----|----|----|
| [`mt_backtest()`](../reference/mt_backtest.md) | rolling-origin evaluation at the resolution of the observed deaths | Tashman (2000) | published; a cell without deaths counts as half a death (package convention) |
| CRPS from draws | empirical-CDF estimator | Gneiting and Raftery (2007); Grimit et al. (2006) | published |
| WIS | weighted interval score with the 23 quantiles of the COVID-19 forecast hubs | Bracher et al. (2021) | published |
| PIT | randomised where draws equal the observation | Czado, Gneiting and Held (2009) | published |
| MASE ([`mt_score()`](../reference/mt_score.md)) | absolute error scaled by the mean absolute year-to-year change of the training data | Hyndman and Koehler (2006) | published |
| [`mt_dm_test()`](../reference/mt_dm_test.md) | Diebold-Mariano test of one loss per origin: equal-weighted cosine variance with fixed-smoothing t critical values (default); the test of [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html) with the correction of Harvey, Leybourne and Newbold for `varestimator = "acf"` (or Bartlett weights, `"bartlett"`) | Diebold and Mariano (1995); Lazarus et al. (2018); Coroneo and Iacone (2020); Harvey, Leybourne and Newbold (1997) | published |

## Results and derived quantities

| Function, option | Method | Source | Status |
|----|----|----|----|
| [`life_table()`](../reference/life_table.md) | `q` from `m` with `a_x`, open interval `L = l / m` | Preston, Heuveline and Guillot (2001) | published |
| [`life_table()`](../reference/life_table.md), `a_x` | `n / 2` for every closed interval (default) or given by `ax`; the age-0 values of Andreev and Kingkade (2015) used by the Human Mortality Database are not implemented | Preston, Heuveline and Guillot (2001); Andreev and Kingkade (2015) | published |
| [`life_table()`](../reference/life_table.md) safeguards | constant-hazard `q = 1 - exp(-n m)` where the formula exceeds 1; the last positive rate for an open interval with a rate of 0 | none | package-specific safeguards |
| `life_table(open = "gompertz", "kannisto")` | closure of the open interval by a Gompertz or Kannisto curve | Thatcher, Kannisto and Vaupel (1998) | published, deviations: fitted by least squares to the log rates (Gompertz) or their logits (Kannisto) of the last `k` ages; the source fits by maximum likelihood |
| [`life_expectancy()`](../reference/life_expectancy.md), [`prob_dying()`](../reference/derived.md) | period life table | Preston, Heuveline and Guillot (2001) | published |
| [`asr()`](../reference/derived.md) | direct standardisation on the age groups of the standard population | Eurostat (2013) | published |
| `asr(age_groups = "grid")` | standard weights shared evenly among the grid cells of each standard age group | none | package-specific, opt-in |
| [`standard_population()`](../reference/standard_population.md) | European Standard Populations 2013 and 1976, WHO World Standard Population | Eurostat (2013); Ahmad et al. (2001) | published |
| [`yll()`](../reference/derived.md) | period expected years of life lost with the remaining life expectancy at the average age at death; standard expected years of life lost with `reference_ex` | Murray (1994) | published, deviations: life expectancy at `x + a_x` interpolated linearly between the bounds of the interval |
| [`mt_observed()`](../reference/mt_observed.md) | Poisson bootstrap of observed indicators | none | package-specific construction (parametric bootstrap) |
| [`mt_contrast()`](../reference/mt_contrast.md) | differences and ratios between strata, path by path, with the share of positive differences | none | package-specific summary of the draws |
| [`mt_inequality()`](../reference/mt_inequality.md) | SII and RII from a weighted linear regression on the relative rank | Mackenbach and Kunst (1997); Moreno-Betancur et al. (2015) | published (identity link) |
| [`lifespan_variation()`](../reference/lifespan_variation.md) | e-dagger, entropy, Gini coefficient, variance, standard deviation, coefficient of variation, interquartile range | Vaupel (1986); Keyfitz (1977); Leser (1955); Goldman and Lord (1986); Hanada (1983); Shkolnikov, Andreev and Begun (2003); LifeIneq | published, deviations: variance and standard deviation measured from the conditioning age `at` (LifeIneq: from the first age of abridged tables); interquartile range by linear interpolation of the survival curve |
| `mt_decompose(method = "arriaga")` | symmetric Arriaga decomposition; causes by their implied sensitivity | Arriaga (1984); Riffe, Tursun-zade and Trias Llimós (2024) | published |
| `mt_decompose(method = "horiuchi")` | line integral of continuous change | Horiuchi, Wilmoth and Pletcher (2008) | published |
| `mt_decompose(method = "stepwise")` | stepwise replacement from the youngest age, averaged over both directions | Andreev, Shkolnikov and Begun (2002) | published, deviations: all causes of an age are replaced together |
| `mt_decompose(method = "kitagawa")` | rate and composition effects of the crude death rate | Kitagawa (1955) | published |
| `mt_decompose(composition = "radix")` | subgroup life tables blended at the radix | Riffe, Tursun-zade and Trias Llimós (2024); Vaupel (2002) | published |
| `mt_decompose(composition = "pooled")` | composition of pooled rates by the line integral | Horiuchi, Wilmoth and Pletcher (2008) for the line integral | package-specific parametrisation |
| [`mt_excess()`](../reference/mt_excess.md) | observed minus predictive deaths of a baseline model | Msemburi et al. (2023) and others compare with the expected deaths of a baseline | published idea; estimand (prediction interval of the realised excess), P-score quantiles and `p_positive` as documented |
