#' Lee-Carter model
#'
#' `log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
#' `sum(k) = 0`, fitted per stratum:
#'
#' * `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
#'   Vermunt 2002), maximum likelihood by Newton updates;
#' * `method = "svd"`: original method (Lee and Carter 1992): singular value
#'   decomposition of centred log rates, with `k[t]` re-estimated in a
#'   second stage (`adjust`, as in `demography::lca()`): so that the fitted
#'   deaths match the observed total deaths of each year (`"dt"`, Lee and
#'   Carter 1992), by Poisson maximum likelihood of the deaths by age given
#'   `a` and `b` (`"dxt"`, Booth, Maindonald and Smith 2002), so that the
#'   fitted life expectancy at the first age matches the observed one
#'   (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).
#'
#' `choose_period = TRUE` fits the model to the last years only, starting in
#' the year that minimises the ratio of the mean deviance of a fit with a
#' linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
#' Maindonald and Smith 2002): the longest recent period over which the
#' decline was close to linear (`demography::lca(chooseperiod = TRUE,
#' breakmethod = "bms")` gives the same period). Fitted rates of earlier
#' years are missing.
#'
#' Compared with `demography::lca()`, `k[t]` is centred again after the
#' second stage (a constant shift of `k` and `a`, with identical fitted
#' rates), and life tables use `a_x = n/2` at age 0 (as [life_table()]), so
#' `adjust = "e0"` matches a slightly different life expectancy than
#' demography's life table where age 0 is in the data.
#'
#' The period index follows a random walk with drift. Draws include the
#' innovations and the uncertainty of the drift. With `bootstrap = TRUE` the
#' parameters are re-estimated on `n_boot` Poisson resamples of the fitted
#' deaths of the fitting period (semi-parametric bootstrap with fitted deaths
#' of Renshaw and Haberman 2008), which adds the uncertainty of `a` and `b`.
#' The bootstrap of Brouhns, Denuit and Van Keilegom (2005) resamples the
#' observed deaths instead: `bootstrapped(lee_carter(), type = "data")`
#' ([bootstrapped()]).
#'
#' Age groups without deaths in the fitting period get a constant rate (0.5
#' deaths over the period) and no trend; with `adjust = "e0"` they enter the
#' life tables at that rate.
#'
#' **Safeguards of this package for sparse counts.** Where the Poisson
#' likelihood has no finite maximum (`method = "poisson"` with very few
#' deaths), the fit falls back to two package-specific estimators (reference
#' unknown), each with a warning (class `mt_warning_model_fit`; within
#' [mt_fit()] it is reported with the stratum): first to the Poisson fit with
#' a ridge penalty that shrinks `b` towards equal shares `1 / n_ages`, and
#' if that also fails to the SVD estimates (`method = "svd"` with `adjust =
#' "dt"`). Their forecasts are not those of the Poisson Lee-Carter model;
#' wider age groups ([mt_regroup()]) or pooling ([li_lee()]) avoid them.
#'
#' Draws describe the expected rates. Predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#' given the rates and the exposures, the observation model of the Poisson
#' Lee-Carter model (Brouhns, Denuit and Vermunt 2002); no other observation
#' noise is added. If intervals are too narrow on real data, that is a
#' property of the model: see [recalibrated()], [ensemble()], or
#' overdispersed models such as `family = "nb"` below or `stanmomo("lc",
#' family = "nb")` ([stanmomo()]).
#'
#' **Jump-off.** The jump-off rates are the rates of the last year `T` of
#' the fitting period, from which the forecast starts. Stoeldraijer et al.
#' (2018) compare six options for the Lee-Carter model; their names are in
#' quotes:
#'
#' * `jump_off = "fitted"`: "model values", `a[x] + b[x] k[T]` (Lee and
#'   Carter 1992);
#' * `jump_off = "observed"`: "last observed values", the observed log rates
#'   of year `T` (Lee and Miller 2001): `a[x]` is replaced by them and `k` is
#'   0 in year `T`;
#' * `jump_off = "observed", jump_off_years = n` (`n` = 2 to 5 in the
#'   paper): "average n years observed". `a[x]` is the mean of the observed
#'   log rates of the last `n` years and `k` is 0 at the midpoint of these
#'   years; the random walk with drift carries the mean from the midpoint to
#'   year `T`, so the jump-off rates are that mean plus `b[x] d (n - 1) / 2`,
#'   with `d` the drift of `k`. The paper does not say how the index at the
#'   midpoint relates to the fitted index of year `T`; this reading (the
#'   forecast does not depend on `k[T]`) reproduces its results for French
#'   men (below), where the mean of the fitted index over the `n` years does
#'   not.
#'
#' Stoeldraijer et al. call a forecast accurate if its out-of-sample errors
#' in historical data are small (mean absolute error of life expectancy at
#' 65 one and five years ahead), and robust if it changes little when one
#' more year of data is added (standard deviation of the changes of the
#' forecast life expectancy at 65 in 2020 between successive fitting
#' periods). In eight Western European countries the last observed values
#' were the most accurate for most countries and an average of observed
#' years the most robust; more years averaged gave more robust and less
#' accurate forecasts. For French men (Human Mortality Database in
#' `StanMoMo::FRMaleData`, a later release than the paper's; ages 0-94 and
#' 95+; fitting periods from 1960 to 2005, ..., 2014; `method = "svd"`), the
#' mean absolute errors and the standard deviations of the last observed
#' values and of the averages are within 0.03 years of Table 8a of the
#' paper. Those of the model values differ by up to 0.08 years: the mean
#' error of the fitted log rates at ages 65 and over in 2005-2014 is -0.002
#' with these data and 0.005 in the paper (Table 9a), and a shift of the log
#' rates at these ages by 0.007 moves life expectancy at 65 by 0.05 years.
#' [anchored()] moves the jump-off of any model.
#'
#' **Filtered jump-off** (`jump_off = "filtered"`; De Jong and Tickle 2006).
#' The paper writes the Lee-Carter model in state-space form (eq. 6): `y[t]
#' = a + b alpha[t] + eps[t]` for the vector of log rates of year `t`, and
#' `alpha[t + 1] = delta + alpha[t] + lambda eta[t]`, with independent normal
#' `eps[t]` (variance `sigma^2` at every age) and `eta[t]` (variance
#' `sigma^2`), `b'b = 1` and the signal-to-noise ratio `lambda`. The forecast
#' starts from `a + b E(alpha[T] | y[1], ..., y[T])`, the Kalman filter
#' estimate of the index in the jump-off year (eq. 15). The paper estimates
#' `b` and `lambda` by maximum likelihood with the Kalman filter, then `a`,
#' `delta` and `sigma` given them (Section 5), for a model with B-spline age
#' effects (eq. 2; eq. 6 is its case without splines, the model fitted
#' here). Here `a` and `b` are those of the Lee-Carter fit (`method`), in
#' the scale `sum(b) = 1` of `k`; given them, `sigma`, `lambda`, `delta`
#' and the level of the index in the first year are estimated by maximum
#' likelihood of the paper's model for the observed log rates. The filtered
#' index is computed in closed form (generalised least squares for the
#' level and the drift and the best linear prediction of the index), which
#' equals the Kalman filter with a diffuse initial state. Without state
#' noise (`lambda = 0`, which is the estimate in some short or noisy series)
#' it is the value in year `T` of the linear trend fitted to the per-year
#' indices `z[t] = b'(y[t] - a) / b'b`; with little
#' observation noise (large `lambda`) it approaches `z[T]`, the projection
#' of the observed log rates of year `T` on the age pattern, which is `k[T]`
#' for `method = "svd", adjust = "none"`. Only the start changes: the
#' forecast keeps the drift and the innovations of `k` (the paper's forecast
#' uses the maximum likelihood drift and adds the variance of the filtered
#' index). Package-specific (reference unknown): cells without deaths enter
#' with 0.5 deaths, as for `jump_off = "observed"`; age groups without deaths
#' in the fitting period are left out; with `bootstrap = TRUE` each
#' bootstrap fit filters the observed rates with its own `a` and `b`.
#' [mt_parameters()] reports the filtered index (`k_filtered`), its standard
#' error (`se_k_filtered`, which includes the uncertainty of the level and
#' the drift) and the estimates `filter_sigma`, `filter_lambda` and
#' `filter_drift` (per year, in the units of `k`). The paper's data
#' (Australian women, 1921-2000) are not in an R package; tests check the
#' limiting cases, a Kalman filter and recovery on simulated data.
#'
#' **Negative binomial version** (`family = "nb"`; Delwarde, Denuit and
#' Partrat 2007). Following the paper:
#'
#' * the deaths are negative binomial with mean `mu = E exp(a[x] + b[x]
#'   k[t])` and variance `mu + phi mu^2`, a Poisson distribution mixed over
#'   a gamma distributed random effect of each cell (Section 2.2). The
#'   dispersion `phi` (the paper's `k`) is one parameter for all ages and
#'   years;
#' * `a`, `b`, `k` and `phi` are estimated by maximum likelihood (eq. 5):
#'   uni-dimensional Newton updates of `a`, `b` and `k` (eq. 6-8), each given
#'   the other parameters, alternated with the estimate of `phi` given the
#'   fitted deaths;
#' * `k[t]` follows a random walk with drift (Section 3.5);
#' * predictive draws (`mt_forecast(predictive = TRUE)`) are negative
#'   binomial deaths with the estimated dispersion (size `1 / phi`), the
#'   observation model of Section 2.2.
#'
#' Package-specific (reference unknown): the Newton updates start from the
#' Poisson estimates (the paper: `a = 0`, `b = 1`, `k = 0`, or random
#' values); `phi` is maximised by Newton steps on `log(phi)` with step
#' halving (the paper: Newton steps on `1 / phi`, eq. 9); the estimate is 0
#' (the Poisson model) where the derivative of the log-likelihood at `phi =
#' 0` is not positive; the log-likelihood uses gamma functions, so deaths
#' need not be whole numbers. The paper states no identification
#' constraints; `sum(b) = 1` and `sum(k) = 0` of the Poisson model
#' (Brouhns, Denuit and Vermunt 2002) agree with its Figures 2 and 3. As for
#' the Poisson model, draws include the innovations of `k` and the
#' uncertainty of its drift; the paper reports point forecasts only. With
#' `bootstrap = TRUE` the resampled deaths are negative binomial with the
#' estimated dispersion, and each draw has the dispersion of its bootstrap
#' fit (package-specific). The ridge and SVD fallbacks are not used: a fit
#' that does not converge fails. Like other maximum likelihood estimates of
#' a variance, `phi` is biased downwards when many parameters are estimated
#' relative to the number of cells. For France, ages 0-99, the fit to the
#' Human Mortality Database series in `demography::fr.mort` reproduces the
#' paper's dispersions (Section 3.3), AIC (Section 3.4) and drifts and
#' standard deviations of `k` (Sections 3.5 and 3.6) within 1%; the paper's
#' column labelled `sigma^2` holds the standard deviation (divisor `T - 1`).
#'
#' **Dispersion by age** (`family = "nb"` with `dispersion = "age"`; Li,
#' Hardy and Tan 2009). Following the paper:
#'
#' * each cell consists of clusters whose deaths are Poisson given a gamma
#'   distributed frailty with mean 1 and a variance that depends on age. The
#'   deaths of the cell are then negative binomial with mean `mu = E
#'   exp(a[x] + b[x] k[t])` and variance `mu + alpha[x] mu^2` (Section 3, eq.
#'   5-12), with one dispersion `alpha[x]` per age group of the data;
#'   `alpha[x]` near 0 gives the Poisson model;
#' * `a`, `b`, `k` and `alpha` are estimated by maximum likelihood (eq. 13):
#'   the Newton updates of `a`, `b` and `k` of Appendix 1, alternated with
#'   the estimate of each `alpha[x]` given the fitted deaths of its age
#'   group;
#' * `jump_off = "constrained"` applies the forecast origin correction
#'   (Section 3 and Appendix 2): `a[x]` is the observed log rate of the last
#'   year `T` (the paper writes `a[x] = m[x, T]`; `a` is on the log scale of
#'   eq. 1), `k[T] = 0`, and `b` and `k` are estimated by maximum
#'   likelihood under these constraints. The fitted rates of the last year
#'   are the observed ones, and the log-likelihood is lower than without the
#'   constraints (footnote 3). With `family = "poisson"` it gives the
#'   original model with the correction, which the paper fits for comparison
#'   (Section 4.2);
#' * interval forecasts come from the parametric bootstrap of Section 3:
#'   `bootstrap = TRUE` simulates negative binomial deaths from the fitted
#'   model, re-estimates `a`, `b` and `k` (with `jump_off = "constrained"`,
#'   `a` comes from the simulated deaths of the last year, as the correction
#'   is part of the estimation in Section 3), and simulates
#'   `k` from the time series model of each refit. The dispersions are those
#'   of the fit: step 2(a) re-estimates `a`, `b` and `k`. The paper's method
#'   is `lee_carter(family = "nb", dispersion = "age", jump_off =
#'   "constrained", bootstrap = TRUE)`;
#' * predictive draws (`mt_forecast(predictive = TRUE)`) are negative
#'   binomial deaths with the dispersion of their age group (size `1 /
#'   alpha[x]`), the observation model of eq. (10). The paper's intervals
#'   describe the central death rates, which are the draws without
#'   `predictive`.
#'
#' Package-specific (reference unknown): as for one dispersion, the updates
#' start from the Poisson estimates (Appendix 1 recommends them), and each
#' `alpha[x]` is maximised by Newton steps on `log(alpha[x])` with step
#' halving (the paper: Newton steps on `alpha[x]`). Its estimate is 0 where
#' the derivative of the log-likelihood at 0 is not positive. The paper fits
#' an ARIMA model to `k` without naming the model it chose; here `k` follows
#' a random walk with drift, and the draws of each bootstrap fit include the
#' uncertainty of its drift, as all draws of `lee_carter()`. With `jump_off =
#' "constrained"`, an age group without deaths in the last year starts from
#' the rate of 0.5 deaths. Age groups without deaths in all years keep their
#' constant rate and have Poisson deaths (dispersion 0). As with one
#' dispersion, a fit
#' that does not converge fails; `jump_off = "constrained"` has no
#' fallbacks either. Each `alpha[x]` is estimated from cells that also
#' estimate `a[x]`, `b[x]` and `k`. Like other maximum likelihood estimates
#' of a variance it is biased downwards: by about `2 / T` with `T` years,
#' and more at the ages that carry most of the information on `k` (many
#' deaths and a small dispersion). The paper uses data of the Human
#' Mortality Database for the United States and Canada, ages 0-99,
#' 1950-2004 (Section 4.1). The packages demography and StMoMo do not hold
#' these series, so its estimates are not reproduced.
#'
#' `rotate = TRUE` rotates the age pattern of decline as life expectancy
#' rises (Li, Lee and Gerland 2013; status: published, with the deviations
#' listed below). In long-lived populations the decline of mortality slows
#' down at young ages and speeds up at old ages, which a fixed `b[x]` does
#' not capture over long horizons. For each path and year `t`, the life
#' expectancy at birth `e0` of the unrotated forecast gives the weight `w =
#' (0.5 (1 + sin(pi / 2 (2 u - 1))))^power` with `u = (e0 - e0_start) /
#' (e0_end - e0_start)`; `w` is 0 below `e0_start` and 1 from `e0_end` (eq.
#' 2 of the paper). The age pattern of that year is `B[x, t] = (1 - w) b[x] +
#' w bu[x]` (eq. 3). The defaults are those of the paper: the rotation starts
#' at 80 (about the median life expectancy of the 20 populations with the
#' highest life expectancy in 2005-2010), is complete at 102 (the value with
#' which the ratio of the death rates at age 0 and at ages 15-19 approaches
#' its estimated long-run level of 7.7), and the exponent is 0.5. The
#' ultimate pattern `bu` has the mean of `b` over ages 15-64 at all ages
#' below 65, as in the paper. From age 65 it keeps the shape of `b`, scaled
#' so that the first age group from 65 equals that mean, and it sums to 1
#' like `b`: the paper only notes that lowering `b` at young ages raises it at
#' older ages because `b` sums to 1, and this scaling is that of MortCast
#' (Sevcikova et al. 2016). The index is then solved again, path by path and
#' year by year, so that the rotated rates have the same life expectancy at
#' birth as the unrotated ones (eq. 4): rotation changes the age pattern of
#' the forecast and keeps its life expectancy at birth. Decline thus moves
#' from young ages to the ages where `bu` exceeds `b`. Where `b` at 65 is
#' already above its mean over ages 15-64 (as for French women), it moves
#' mostly to ages 15-64. Fitted rates, and paths with `e0` below `e0_start`,
#' are unchanged. Rotation needs ages from 0 (life expectancy at birth) with
#' age groups starting at 15 and 65. With a closed last age group, a warning
#' notes that life expectancy is then truncated at the highest age. Rotated
#' forecasts compute life tables for every path and year and are much slower
#' (about 20 times for single years of age).
#'
#' Deviations of the rotation from Li, Lee and Gerland (2013):
#'
#' * the paper forecasts `b` and `k` for both sexes combined (Li and Lee
#'   2005) and takes the weights from the life expectancy of both sexes
#'   combined; `lee_carter()` fits each stratum separately, so the weights
#'   come from the life expectancy of the stratum (MortCast takes the path
#'   of life expectancy from a separate model);
#' * the paper illustrates the method with the median path of life
#'   expectancy; here every draw is rotated with its own path of life
#'   expectancy;
#' * rates change from the jump-off year `T`, `log m[x, t] = log m[x, T] +
#'   B[x, t] (k*[t] - k[T])`. Eq. (4) writes `a[x] + B[x, t] K[t]` with `a`
#'   the last observed log rates (as `jump_off = "observed"`) or their mean
#'   over the years; with `jump_off = "fitted"` (or `"filtered"`) the fitted
#'   (or filtered) rates of year `T` are the starting point
#'   (package-specific), which keeps the jump-off rates when `B` departs
#'   from `b`;
#' * the paper extends the open age group to age 130 with a logistic model
#'   before fitting (footnote 3); here the age groups of the data are used,
#'   and life expectancy closes the open interval with `L = l / m`
#'   ([life_table()]).
#'
#' @param method `"poisson"` (maximum likelihood with the distribution of
#'   `family`) or `"svd"`.
#' @param adjust Re-estimation of `k[t]` for `method = "svd"`: `"dt"`,
#'   `"dxt"`, `"e0"` or `"none"`.
#' @param choose_period Choose the fitting period (Booth, Maindonald and
#'   Smith 2002).
#' @param min_years `choose_period` considers fitting periods of more than
#'   `min_years` years (`minperiod` of `demography::lca()`).
#' @param jump_off `"fitted"` (forecasts start from the fitted rates of the
#'   last year), `"observed"` (from the observed rates of the last year; Lee
#'   and Miller 2001; with `jump_off_years > 1` from their mean over the last
#'   years, Stoeldraijer et al. 2018), `"filtered"` (from the rates at the
#'   Kalman filter estimate of the index in the last year; De Jong and
#'   Tickle 2006) or `"constrained"` (the maximum likelihood fit is
#'   constrained to the observed rates of the last year, the forecast origin
#'   correction of Li, Hardy and Tan 2009; `method = "poisson"` only). See
#'   "Jump-off" above.
#' @param jump_off_years For `jump_off = "observed"`: the number of last
#'   observed years whose mean log rates give the jump-off rates (default 1,
#'   the last year; 2 to 5 in Stoeldraijer et al. 2018).
#' @param bootstrap Semi-parametric bootstrap of the parameters.
#' @param n_boot Number of bootstrap fits.
#' @param rotate Rotate the age pattern of decline `b` towards an ultimate
#'   pattern as life expectancy at birth rises (Li, Lee and Gerland 2013).
#' @param family Distribution of the deaths for `method = "poisson"`:
#'   `"poisson"` (Brouhns, Denuit and Vermunt 2002) or `"nb"` (negative
#'   binomial; Delwarde, Denuit and Partrat 2007).
#' @param dispersion Dispersion of the negative binomial deaths (`family =
#'   "nb"`): `"global"` (one for all ages; Delwarde, Denuit and Partrat
#'   2007) or `"age"` (one per age group of the data; Li, Hardy and Tan
#'   2009).
#' @param rotation Settings of the rotation, a named vector: life expectancy
#'   at birth at which the rotation starts (`e0_start`, default 80) and is
#'   complete (`e0_end`, default 102), and the exponent of the weights
#'   (`power`, default 0.5). The defaults are those of Li, Lee and Gerland
#'   (2013, eq. 2) and MortCast; elements left out keep their defaults.
#' @return A model specification.
#' @references
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear regression
#' approach to the construction of projected lifetables. *Insurance:
#' Mathematics and Economics* 31(3): 373-393.
#' \doi{10.1016/S0167-6687(02)00185-3}
#'
#' Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
#' conditions of variable mortality decline. *Population Studies* 56(3):
#' 325-336. \doi{10.1080/00324720215935}
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#' \doi{10.1353/dem.2001.0036}
#'
#' Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
#' strategies for matching mortality forecasts to the most recently observed
#' data: exploring the trade-off between accuracy and robustness. *Genus*
#' 74: 16. \doi{10.1186/s41118-018-0041-y}
#'
#' De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting.
#' *Mathematical Population Studies* 13(1): 1-18.
#' \doi{10.1080/08898480500452109}
#'
#' Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the
#' Lee-Carter model for mortality forecasting. *Applied Stochastic Models in
#' Business and Industry* 23(5): 385-401. \doi{10.1002/asmb.679}
#'
#' Li JS-H, Hardy MR, Tan KS (2009). Uncertainty in mortality forecasting:
#' an extension to the classical Lee-Carter approach. *ASTIN Bulletin* 39(1):
#' 137-164. \doi{10.2143/AST.39.1.2038060}
#'
#' Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
#' log-bilinear model for mortality forecasting. *Scandinavian Actuarial
#' Journal* 2005(3): 212-224.
#'
#' Renshaw AE, Haberman S (2008). On simulation-based approaches to risk
#' measurement in mortality with specific reference to Poisson Lee-Carter
#' modelling. *Insurance: Mathematics and Economics* 42(2): 797-816.
#'
#' Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
#' the rotation of age patterns of mortality decline for long-term
#' projections. *Demography* 50(6): 2037-2051.
#' \doi{10.1007/s13524-013-0232-2}
#'
#' Li N, Lee R (2005). Coherent mortality forecasts for a group of
#' populations: an extension of the Lee-Carter method. *Demography* 42(3):
#' 575-594. \doi{10.1353/dem.2005.0021}
#'
#' Sevcikova H, Li N, Kantorova V, Gerland P, Raftery AE (2016). Age-specific
#' mortality and fertility rates for probabilistic population projections.
#' In Schoen R (ed.) *Dynamic Demographic Analysis*. Springer, 285-310.
#' \doi{10.1007/978-3-319-26603-9_15}
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
#' \donttest{
#' # rotation needs ages from 0: French women (HMD, via demography)
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   fr <- as_mt_data(demography::fr.mort, series = "female", years = 1970:2006)
#'   fr <- mt_regroup(fr, age = c(0, 1, seq(5, 100, 5)))
#'   fc_rot <- mt_forecast(mt_fit(lee_carter(rotate = TRUE), fr), h = 50, n_draws = 100, seed = 1)
#'   life_expectancy(fc_rot)
#' }
#' }
#' # negative binomial deaths (Delwarde, Denuit and Partrat 2007)
#' fit_nb <- mt_fit(lee_carter(family = "nb"), x)
#' mt_parameters(fit_nb)[parameter == "nb_dispersion"]
#' # dispersion by age with the forecast origin correction and the parametric
#' # bootstrap (Li, Hardy and Tan 2009)
#' fit_lht <- mt_fit(lee_carter(family = "nb", dispersion = "age", jump_off = "constrained",
#'                              bootstrap = TRUE, n_boot = 20), mt_subset(x, sex == "female"))
#' mt_parameters(fit_lht)[parameter == "nb_dispersion"]
#' # jump-off: mean of the last three observed years (Stoeldraijer et al.
#' # 2018) and the Kalman filter estimate (De Jong and Tickle 2006)
#' fit_av <- mt_fit(lee_carter(jump_off = "observed", jump_off_years = 3), x)
#' fit_kf <- mt_fit(lee_carter(jump_off = "filtered"), x)
#' mt_parameters(fit_kf)[grepl("filter", parameter)]
#' @export
lee_carter <- function(method = c("poisson", "svd"), jump_off = c("fitted", "observed", "filtered", "constrained"),
                       bootstrap = FALSE, n_boot = 100, adjust = c("dt", "dxt", "e0", "none"), choose_period = FALSE,
                       min_years = 20, rotate = FALSE, rotation = c(e0_start = 80, e0_end = 102, power = 0.5),
                       family = c("poisson", "nb"), dispersion = c("global", "age"), jump_off_years = 1) {
  method <- match.arg(method)
  jump_off <- match.arg(jump_off)
  adjust <- match.arg(adjust)
  family <- match.arg(family)
  dispersion <- match.arg(dispersion)
  check_flag(choose_period, "choose_period")
  check_number(min_years, "min_years", lower = 5, integerish = TRUE)
  check_number(jump_off_years, "jump_off_years", lower = 1, integerish = TRUE)
  if (jump_off_years > 1 && jump_off != "observed") {
    mt_abort("`jump_off_years` averages observed log rates; use it with `jump_off = \"observed\"`.",
             class = "argument")
  }
  check_flag(rotate, "rotate")
  rotation <- rotation_settings(rotation)
  if (family == "nb" && method != "poisson") {
    mt_abort(paste("`family = \"nb\"` is estimated by maximum likelihood (Delwarde, Denuit and Partrat 2007);",
                   "use it with `method = \"poisson\"`."), class = "argument")
  }
  if (family == "nb" && choose_period) {
    mt_abort(paste("`choose_period` compares Poisson deviances (Booth, Maindonald and Smith 2002) and is not",
                   "available with `family = \"nb\"`."), class = "argument")
  }
  if (!identical(dispersion, "global") && family != "nb") {
    mt_abort("`dispersion` sets the dispersion of the negative binomial deaths; use it with `family = \"nb\"`.",
             class = "argument")
  }
  if (jump_off == "constrained" && (method != "poisson" || choose_period)) {
    mt_abort(paste("`jump_off = \"constrained\"` constrains the maximum likelihood fit (Li, Hardy and Tan 2009);",
                   "use it with `method = \"poisson\"` and without `choose_period`."), class = "argument")
  }
  new_mt_spec("lee_carter", list(method = method, jump_off = jump_off, bootstrap = bootstrap, n_boot = n_boot,
                                 adjust = adjust, choose_period = choose_period, min_years = min_years,
                                 rotate = rotate, rotation = rotation, family = family, dispersion = dispersion,
                                 jump_off_years = jump_off_years),
              min_periods = max(if (choose_period) min_years else 5, jump_off_years),
              label = paste0("Lee-Carter", if (method == "svd") sprintf(" (SVD, %s)", adjust) else "",
                             if (family == "nb") " (negative binomial" else "",
                             if (family == "nb" && dispersion == "age") ", dispersion by age" else "",
                             if (family == "nb") ")" else "",
                             if (rotate) rotation_label(rotation) else "",
                             if (choose_period) ", chosen period" else "",
                             switch(jump_off, fitted = "", filtered = ", filtered jump-off",
                                    constrained = ", constrained jump-off",
                                    observed = if (jump_off_years > 1) {
                                      sprintf(", observed jump-off (mean of %d years)", jump_off_years)
                                    } else {
                                      ", observed jump-off"
                                    }),
                             if (bootstrap) ", bootstrap" else ""))
}


# Settings of the rotation: defaults completed by the elements given
rotation_settings <- function(rotation) {
  out <- c(e0_start = 80, e0_end = 102, power = 0.5)
  if (is.list(rotation)) rotation <- unlist(rotation)
  nm <- names(rotation)
  if (!is.numeric(rotation) || (length(rotation) && (is.null(nm) || !all(nm %in% names(out)) || anyDuplicated(nm)))) {
    mt_abort("`rotation` must be a named numeric vector with elements `e0_start`, `e0_end` and `power`.",
             class = "argument", call = NULL)
  }
  out[nm] <- rotation
  if (anyNA(out) || any(!is.finite(out)) || out[["e0_end"]] <= out[["e0_start"]] || out[["power"]] <= 0) {
    mt_abort("`rotation` needs finite values with `e0_start` < `e0_end` and `power` > 0.", class = "argument",
             call = NULL)
  }
  out
}

rotation_label <- function(rotation) {
  default <- c(e0_start = 80, e0_end = 102, power = 0.5)
  if (isTRUE(all.equal(rotation, default))) return(", rotated")
  sprintf(", rotated (e0 %g-%g, power %g)", rotation[["e0_start"]], rotation[["e0_end"]], rotation[["power"]])
}

#' @export
fit_model.mt_spec_lee_carter <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  if (isTRUE(p$rotate)) rotation_check(ages)
  zero <- rowSums(D) == 0
  if (sum(!zero) < 2) mt_abort("Fewer than two age groups with deaths.", class = "model_fit")
  fallback <- character()
  # dispersion group of each age group, one per age group with a dispersion
  # by age (NULL: one dispersion), and the forecast origin correction (Li,
  # Hardy and Tan 2009)
  groups <- if (identical(p$dispersion, "age")) seq_along(ages$lower)
  constrained <- identical(p$jump_off, "constrained")
  # Poisson estimates; where the Poisson likelihood has no finite maximum
  # (very sparse counts) the package's safeguards: a ridge penalty on b, then
  # SVD estimates (`fallback` records which). `fixed`: constant rates of the
  # age groups left out of D (NA for the rows of D), for life tables over all
  # ages. `disp`: dispersions kept fixed (one per row of D), `grp`: the
  # dispersion groups of the rows of D
  est_one <- function(D, E, fixed = NULL, disp = NULL, grp = NULL) {
    if (p$method == "svd") return(lc_svd(D, E, adjust = p$adjust %||% "dt", ages = ages, fixed = fixed))
    # negative binomial maximum likelihood, without the safeguards of the
    # Poisson fit
    if (identical(p$family, "nb")) return(lc_nb(D, E, disp = disp, groups = grp, origin = constrained))
    if (constrained) {
      # the Poisson model with the forecast origin correction: the limit of
      # the negative binomial model with dispersion 0
      out <- lc_nb(D, E, disp = 0, origin = TRUE)
      out$method <- "poisson"
      out$dispersion <- out$nb_size <- NULL
      return(out)
    }
    tryCatch(lc_poisson(D, E), mt_error_model_fit = function(e) {
      tryCatch({
        out <- lc_poisson(D, E, ridge = 1)
        fallback <<- c(fallback, "ridge")
        out
      }, mt_error_model_fit = function(e) {
        fallback <<- c(fallback, "svd")
        lc_svd(D, E)
      })
    })
  }
  # `disp`: dispersions of all age groups kept fixed (bootstrap refits)
  fit_lc <- function(D, E, disp = NULL) {
    if (!any(zero)) return(est_one(D, E, disp = disp, grp = groups))
    # age groups without any death: constant rate of 0.5 deaths over the
    # training years, no trend (the Poisson estimate would be minus infinity)
    fixed <- rep(NA_real_, nrow(D))
    fixed[zero] <- 0.5 / rowSums(E[zero, , drop = FALSE])
    grp <- if (!is.null(groups)) match(groups[!zero], unique(groups[!zero]))
    e <- est_one(D[!zero, , drop = FALSE], E[!zero, , drop = FALSE], fixed, disp = disp[!zero], grp = grp)
    a <- b <- numeric(nrow(D))
    a[!zero] <- e$a
    b[!zero] <- e$b
    a[zero] <- log(fixed[zero])
    e$a <- a
    e$b <- b
    if (!is.null(e$age_nb_size)) {
      # dispersion 0 (Poisson) for the age groups without deaths
      d <- numeric(nrow(D))
      d[!zero] <- e$dispersion
      e$dispersion <- d
      e$age_nb_size <- 1 / d
    }
    e
  }
  first <- 1L
  if (isTRUE(p$choose_period)) {
    first <- bms_start(D, E, est_one, p$min_years)
    D <- D[, first:ncol(D), drop = FALSE]
    E <- E[, first:ncol(E), drop = FALSE]
    zero <- rowSums(D) == 0
  }
  # fallbacks of the final fit only (not of the candidate periods)
  fallback <- character()
  est <- fit_lc(D, E)
  if (any(zero)) {
    mt_warn(sprintf(paste("No deaths in age group(s) %s in the training years: constant rate (0.5 deaths in all",
                          "years), no trend."), paste(ages$labels[zero], collapse = ", ")), class = "sparse")
  }
  if ("ridge" %in% fallback) {
    mt_warn(paste("The Poisson likelihood of the Lee-Carter model has no finite maximum (sparse counts); fitted",
                  "with a ridge penalty that shrinks the age pattern of change b towards equal shares (a",
                  "safeguard of this package without a published source)."), class = "model_fit")
  }
  if ("svd" %in% fallback) {
    mt_warn(paste("The Poisson likelihood of the Lee-Carter model has no finite maximum (sparse counts), also with",
                  "a ridge penalty; the SVD estimates (method = \"svd\", adjust = \"dt\") are used instead (a",
                  "safeguard of this package)."), class = "model_fit")
  }
  # observed jump-off: the log rates of the last year, or their mean over the
  # last `jump_off_years` years (Stoeldraijer et al. 2018)
  n_jo <- p$jump_off_years %||% 1
  if (n_jo > ncol(D)) {
    mt_abort(sprintf("`jump_off_years = %d` exceeds the %d years of the fitting period.", n_jo, ncol(D)),
             class = "model_data")
  }
  last <- seq(ncol(D) - n_jo + 1, ncol(D))
  est$last_log_rate <- if (n_jo == 1) {
    log_rates(D[, ncol(D)], E[, ncol(E)])
  } else {
    rowMeans(log_rates(D[, last, drop = FALSE], E[, last, drop = FALSE]))
  }
  est$jump_off_years <- n_jo
  est$fitted <- exp(est$a + outer(est$b, est$k))
  est$jump_off <- p$jump_off
  # filtered jump-off (De Jong and Tickle 2006): the observed log rates of the
  # age groups with deaths, filtered with the age effects of each fit
  filtered <- identical(p$jump_off, "filtered")
  if (filtered) {
    y_obs <- log_rates(D[!zero, , drop = FALSE], E[!zero, , drop = FALSE])
    add_filter <- function(s) {
      s$filter <- lc_filter(y_obs, s$a[!zero], s$b[!zero])
      s
    }
    est <- add_filter(est)
  }
  if (isTRUE(p$bootstrap)) {
    # resamples of the (chosen) fitting period: Poisson deaths, or negative
    # binomial deaths with the estimated dispersion (family = "nb"); with a
    # dispersion by age the refits keep it (Li, Hardy and Tan 2009, Section 3,
    # step 2a, re-estimates a, b and k)
    fitted <- E * est$fitted
    keep <- if (!is.null(est$age_nb_size)) est$dispersion
    est$boot <- lapply(seq_len(p$n_boot), function(i) {
      Db <- matrix(r_deaths(as.vector(fitted), est$age_nb_size %||% est$nb_size), nrow(fitted))
      tryCatch(fit_lc(Db, E, disp = keep), error = function(e) NULL)
    })
    est$boot <- Filter(Negate(is.null), est$boot)
    if (filtered) est$boot <- lapply(est$boot, add_filter)
  }
  if (isTRUE(p$rotate)) {
    est$rotation <- c(list(age = ages), as.list(p$rotation))
    est$b_ultimate <- lc_ultimate_b(est$b, ages, zero)
    if (is.null(est$b_ultimate)) {
      mt_abort(paste("Rotation needs a positive age pattern of decline b at ages 15-64 (mean) and in the first age",
                     "group from 65."), class = "model_fit")
    }
    # bootstrap fits without a valid ultimate pattern are dropped
    if (length(est$boot)) {
      est$boot <- Filter(Negate(is.null), lapply(est$boot, function(s) {
        s$b_ultimate <- lc_ultimate_b(s$b, ages, zero)
        if (is.null(s$b_ultimate)) NULL else s
      }))
    }
  }
  if (first > 1L) {
    # fitted rates only for the chosen period
    est$fitted <- cbind(matrix(NA_real_, nrow(D), first - 1L), est$fitted)
    est$first_year <- periods[first]
  }
  est
}

#' @export
simulate_model.mt_spec_lee_carter <- function(spec, fit, h, n_draws, ...) {
  sets <- if (length(fit$boot)) fit$boot else list(fit)
  which_set <- rep_len(seq_along(sets), n_draws)
  out <- array(NA_real_, c(length(fit$a), h, n_draws))
  for (i in seq_along(sets)) {
    js <- which(which_set == i)
    if (!length(js)) next
    s <- sets[[i]]
    kp <- rw_paths(s$k[length(s$k)], s$drift, matrix(s$dk - s$drift, nrow = 1), h, length(js), "gaussian",
                   sigma = s$sigma_k, exp_out = FALSE)
    # paths of k start from k[T]; the jump-off rates of the set are
    # base + b k[T]
    base <- if (fit$jump_off == "fitted") s$a else lc_jump_log_rate(fit, s) - s$b * s$k[length(s$k)]
    for (m in seq_along(js)) {
      kk <- kp[1, , m]
      out[, , js[m]] <- exp(base + outer(s$b, kk))
    }
    if (!is.null(fit$rotation)) out[, , js] <- lc_rotate(out[, , js, drop = FALSE], matrix(kp[1, , ], h), s, fit)
  }
  # negative binomial family: the size 1 / dispersion of the parameter set of
  # each draw, for predictive draws of the deaths (see predictive_rates());
  # with a dispersion by age one size per age group [age, draw], the same in
  # all draws (bootstrap refits keep the dispersions)
  if (!is.null(fit$age_nb_size)) {
    attr(out, "nb_size") <- matrix(fit$age_nb_size, length(fit$a), n_draws)
  } else if (!is.null(fit$nb_size)) {
    attr(out, "nb_size") <- vapply(sets, function(s) s$nb_size, 0)[which_set]
  }
  out
}

# Rotation of the age pattern of mortality decline (Li, Lee and Gerland
# 2013). rotation_weights() and lc_ultimate_b() are adapted from MortCast
# 2.8-0 (rotate.leecarter() and ultimate.bx() in R/LC.R), GPL (>= 2), by Hana
# Sevcikova, Nan Li and Patrick Gerland. Changes: the ultimate pattern is
# defined on any age axis with groups starting at 0, 15 and 65 (the mean
# over ages 15-64 is weighted by the widths of the age groups, which gives
# MortCast's mean for 5-year groups and single years), age groups without
# deaths keep no trend, and the sum is that of b. In MortCast the index is
# solved by bisection (to 0.01 years) for an external path of life
# expectancy; here rotate_rates() solves it for the life expectancy of each
# unrotated Lee-Carter path, and rates change from the jump-off year.

# Ages from 0 (life expectancy at birth) with age groups starting at 15 and
# 65 (the ultimate pattern)
rotation_check <- function(age) {
  lo <- age$lower
  if (lo[1] != 0 || !all(c(15, 65) %in% lo)) {
    mt_abort(sprintf(paste("lee_carter(rotate = TRUE) needs ages from 0 with age groups starting at 15 and 65",
                           "(life expectancy at birth, ultimate pattern of decline); the grid has ages %s to %s.",
                           "Extend or regroup the ages with mt_regroup()."), age$labels[1], age$labels[length(age)]),
             class = "model_data", call = NULL)
  }
  if (is.finite(age$upper[length(age)])) {
    mt_warn(sprintf(paste("The last age group (%s) is closed: the rotation uses life expectancy up to age %s, below",
                          "life expectancy at birth."), age$labels[length(age)], format(age$upper[length(age)])),
            class = "model_data", call = NULL)
  }
  invisible(age)
}

# Weights of the ultimate pattern for life expectancies e0 (MortCast's
# rotate.leecarter(): 0 below e0_start, 1 from e0_end)
rotation_weights <- function(e0, e0_start = 80, e0_end = 102, power = 0.5) {
  u <- (e0 - e0_start) / (e0_end - e0_start)
  w <- (0.5 * (1 + sin(pi / 2 * (2 * u - 1))))^power
  w[!is.na(e0) & e0 < e0_start] <- 0
  w[!is.na(e0) & e0 >= e0_end] <- 1
  w
}

# Ultimate age pattern of decline (MortCast's ultimate.bx()): the mean of b
# over ages 15-64 at ages below 65; from 65 the shape of b, scaled so that
# the first age group from 65 equals that mean; the same sum as b. Age groups
# without deaths (`fixed`) keep b = 0. NULL where the mean or b at 65 is not
# positive.
lc_ultimate_b <- function(b, age, fixed = NULL) {
  fixed <- fixed %||% rep(FALSE, length(b))
  lo <- age$lower
  n <- age$upper - lo
  young <- lo >= 15 & age$upper <= 65 & !fixed
  i65 <- which(lo == 65)
  if (!any(young) || length(i65) != 1L) return(NULL)
  m <- sum(b[young] * n[young]) / sum(n[young])
  if (!isTRUE(m > 0) || !isTRUE(b[i65] > 0)) return(NULL)
  bu <- ifelse(lo >= 65, b * (m / b[i65]), m)
  bu[fixed] <- 0
  bu / sum(bu) * sum(b)
}

# Rotated rates for the unrotated draws R [age, h, draw] of parameter set s,
# with index paths K [h, draw]. Rates change from the jump-off year: log m =
# J + B (k* - k[T]), J the jump-off log rates (lc_jump_log_rate())
lc_rotate <- function(R, K, s, fit) {
  kj <- s$k[length(s$k)]
  J <- lc_jump_log_rate(fit, s)
  d <- dim(R)
  array(rotate_rates(matrix(R, d[1]), J, s$b, s$b_ultimate, as.vector(K) - kj, fit$rotation), d)
}

# Rotation of unrotated rates M [age, path-year] = exp(J + b kap0): per
# column, the weight of the ultimate pattern bu from the life expectancy at
# birth of M, and the index solved so that exp(J + B kap) keeps that life
# expectancy. Columns below e0_start are returned as they are. `rot` holds
# the age axis and the settings (fit$rotation).
rotate_rates <- function(M, J, b, bu, kap0, rot) {
  e0 <- life_expectancy(M, rot$age)
  w <- rotation_weights(e0, rot$e0_start, rot$e0_end, rot$power)
  act <- which(!is.na(w) & w > 0)
  if (!length(act)) return(M)
  B <- outer(b, 1 - w[act]) + outer(bu, w[act])
  kap <- lc_solve_e0(J, B, e0[act], kap0[act], rot$age)
  ok <- !is.na(kap)
  if (!all(ok)) {
    mt_warn(sprintf(paste("Rotation: no index matches the life expectancy of %d of %d path-years; their rates are",
                          "not rotated."), sum(!ok), length(ok)), class = "model", call = NULL)
  }
  M[, act[ok]] <- exp(J + B[, ok, drop = FALSE] * rep(kap[ok], each = nrow(M)))
  M
}

# Index kap (one per column of B) such that the life expectancy at birth of
# exp(J + B[, j] kap[j]) equals e0[j]: secant steps for all columns at once
# from `start`, uniroot() for the columns that do not converge; NA where no
# root is found
lc_solve_e0 <- function(J, B, e0, start, age, tol = 1e-10, maxit = 30) {
  na <- length(J)
  f <- function(kap, cols) {
    life_expectancy(exp(J + B[, cols, drop = FALSE] * rep(kap, each = na)), age) - e0[cols]
  }
  out <- rep(NA_real_, length(start))
  x0 <- start
  f0 <- f(x0, seq_along(x0))
  conv <- is.finite(f0) & abs(f0) < tol
  out[conv] <- x0[conv]
  x1 <- x0 - 0.01 * (1 + abs(x0))
  act <- which(!conv & is.finite(f0))
  for (it in seq_len(maxit)) {
    if (!length(act)) break
    xa <- x1[act]
    f1 <- f(xa, act)
    conv <- is.finite(f1) & abs(f1) < tol
    out[act[conv]] <- xa[conv]
    step <- f1 * (xa - x0[act]) / (f1 - f0[act])
    x0[act] <- xa
    f0[act] <- f1
    x1[act] <- xa - step
    act <- act[!conv & is.finite(step)]
  }
  for (j in which(is.na(out))) {
    out[j] <- tryCatch(stats::uniroot(function(kap) f(kap, j), start[j] + c(-1, 1) * (1 + abs(start[j])),
                                      extendInt = "yes", tol = 1e-12)$root,
                       error = function(e) NA_real_)
  }
  out
}

# Jump-off log rates (year T) of parameter set s: the fitted rates (also for
# the fit constrained to the last year); the observed log rates of the last
# year, or the mean of the last n years moved from their midpoint to T by the
# drift (Stoeldraijer et al. 2018); or the rates at the filtered index (De
# Jong and Tickle 2006)
lc_jump_log_rate <- function(fit, s) {
  switch(fit$jump_off,
         observed = {
           n <- fit$jump_off_years %||% 1
           if (n > 1) fit$last_log_rate + s$b * s$drift * (n - 1) / 2 else fit$last_log_rate
         },
         filtered = s$a + s$b * s$filter$level,
         s$a + s$b * s$k[length(s$k)])
}

# Identification sum(b) = 1 and sum(k) = 0 (rates unchanged); `centre =
# FALSE` keeps a and the level of k (forecast origin correction: k[T] = 0)
lc_finish <- function(a, b, k, centre = TRUE) {
  s <- sum(b)
  kbar <- if (centre) mean(k) else 0
  a <- a + b * kbar
  k <- (k - kbar) * s
  b <- b / s
  dk <- diff(k)
  drift <- mean(dk)
  list(a = a, b = b, k = k, dk = dk, drift = drift, sigma_k = sqrt(sum((dk - drift)^2) / max(length(dk) - 1, 1)))
}

# `fixed`: for adjust = "e0" with rows of `ages` left out of D, their constant
# rates (NA for the rows of D), so that life tables cover the whole age range
lc_svd <- function(D, E, adjust = "dt", ages = NULL, fixed = NULL) {
  y <- log_rates(D, E)
  a <- rowMeans(y)
  s <- svd(y - a, nu = 1, nv = 1)
  b <- s$u[, 1]
  k <- s$d[1] * s$v[, 1]
  if (sum(b) < 0) {
    b <- -b
    k <- -k
  }
  # second stage: re-estimate k (Lee and Carter: total deaths; Booth et al.:
  # deaths by age; Lee and Miller: life expectancy)
  if (adjust == "e0" && is.null(ages)) adjust <- "dt"
  all_ages <- function(m) {
    if (is.null(fixed)) return(array(m, c(nrow(m), ncol(m), 1)))
    out <- matrix(fixed, length(fixed), ncol(m))
    out[is.na(fixed), ] <- m
    array(out, c(nrow(out), ncol(out), 1))
  }
  if (adjust == "e0") e_obs <- life_expectancy(all_ages(pmax(D, 0.5) / E), ages)[, 1]
  for (t in seq_len(ncol(D))) {
    if (adjust == "none") break
    if (adjust == "dxt") {
      kt <- k[t]
      for (it in 1:50) {
        mu <- E[, t] * exp(a + b * kt)
        step <- sum((D[, t] - mu) * b) / sum(mu * b^2)
        kt <- kt + step
        if (abs(step) < 1e-10) break
      }
      if (is.finite(kt)) k[t] <- kt
      next
    }
    f <- if (adjust == "dt") {
      function(kt) sum(E[, t] * exp(a + b * kt)) - sum(D[, t])
    } else {
      function(kt) life_expectancy(all_ages(matrix(exp(a + b * kt))), ages)[1, 1] - e_obs[t]
    }
    rng <- k[t] + c(-1, 1) * (10 + 10 * abs(k[t]))
    # without a root in reach, total deaths keep the SVD index (with `b` of
    # both signs the nearest root of total deaths can lie far off and makes
    # the index erratic); life expectancy searches for the nearest root
    k[t] <- tryCatch(stats::uniroot(f, rng, extendInt = "yes", tol = 1e-10)$root, error = function(e) {
      if (adjust == "e0") nearest_root(f, k[t], 0.1 * (1 + abs(k[t]))) else k[t]
    })
  }
  out <- lc_finish(a, b, k)
  out$method <- "svd"
  out
}

# Root of f nearest to x0, for f that is not monotone (life expectancy as a
# function of k when b has both signs, with sparse counts): steps of growing
# size in both directions until the sign changes; x0 if there is none
nearest_root <- function(f, x0, step) {
  f0 <- f(x0)
  if (!is.finite(f0) || f0 == 0) return(x0)
  prev <- c(x0, x0)
  for (d in step * 2^(-4:12)) {
    for (i in 1:2) {
      x <- x0 + c(-1, 1)[i] * d
      fx <- f(x)
      if (is.finite(fx) && sign(fx) != sign(f0)) {
        return(tryCatch(stats::uniroot(f, sort(c(prev[i], x)), tol = 1e-10)$root, error = function(e) x0))
      }
      prev[i] <- x
    }
  }
  x0
}

lc_poisson <- function(D, E, init = NULL, maxit = 500, tol = 1e-8, ridge = 0) {
  if (is.null(init)) init <- lc_svd(D, E)
  a <- init$a
  b <- if (ridge > 0) rep(1 / nrow(D), nrow(D)) else init$b
  k <- init$k
  b0 <- 1 / nrow(D)
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    Dh <- E * exp(a + outer(b, k))
    a <- a + rowSums(D - Dh) / rowSums(Dh)
    Dh <- E * exp(a + outer(b, k))
    k <- k + colSums((D - Dh) * b) / colSums(Dh * b^2)
    a <- a + b * mean(k)
    k <- k - mean(k)
    Dh <- E * exp(a + outer(b, k))
    kk <- matrix(k, nrow(D), ncol(D), byrow = TRUE)
    # a ridge penalty shrinks b towards equal shares where the counts are
    # too sparse for a finite maximum likelihood estimate
    b <- b + (rowSums((D - Dh) * kk) - ridge * (b - b0)) / (rowSums(Dh * kk^2) + ridge)
    Dh <- E * exp(a + outer(b, k))
    ll <- sum(D * log(Dh) - Dh)
    if (!is.finite(ll)) mt_abort("Lee-Carter fit diverged.", class = "model_fit")
    if (abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  out <- lc_finish(a, b, k)
  out$method <- "poisson"
  out$loglik <- ll
  out$iterations <- it
  out
}

# Negative binomial Lee-Carter model (Delwarde, Denuit and Partrat 2007,
# Section 2): deaths with mean mu = E exp(a + b k) and variance mu + phi
# mu^2, one dispersion phi (the paper's k) for all ages and years. Maximum
# likelihood (eq. 5): the uni-dimensional Newton updates of a, b and k of
# the paper (eq. 6-8), each given the others and phi, and phi
# maximised given the fitted deaths (nb_dispersion()), alternated until the
# log-likelihood converges. The updates start from the Poisson estimates
# (the paper: a = 0, b = 1, k = 0, "random values can also be used"). After
# each round, sum(b) = 1 and sum(k) = 0 (rates unchanged). `disp`: a fixed
# dispersion instead of its estimate (phi = 0 gives the Poisson model).
#
# Dispersion by age (Li, Hardy and Tan 2009, Section 3): `groups` gives the
# dispersion group (1, 2, ...) of each row; the deaths of the rows of group
# g have variance mu + alpha[g] mu^2. The updates of a, b and k are those of
# Appendix 1 (the same terms per cell as above, with the dispersion of the
# row), alternated with the estimate of each alpha[g] given the fitted deaths
# of its rows (eq. 13). `disp` is then one fixed dispersion per row (or one
# for all). The result has the dispersion of each row (`dispersion`) and its
# inverse, the size of each row (`age_nb_size`).
#
# `origin = TRUE`: forecast origin correction (Appendix 2): a fixed at the
# observed log rates of the last year (0.5 deaths where there are none) and
# k[T] = 0, b and k estimated under these constraints (sum(b) = 1).
lc_nb <- function(D, E, disp = NULL, groups = NULL, origin = FALSE, maxit = 5000, tol = 1e-12) {
  init <- tryCatch(lc_poisson(D, E), mt_error_model_fit = function(e) lc_svd(D, E))
  a <- init$a
  b <- init$b
  k <- init$k
  nt <- ncol(D)
  if (origin) {
    a <- log_rates(D[, nt], E[, nt])
    k <- k - k[nt]
  }
  fit_mu <- function() E * exp(a + outer(b, k))
  # dispersion: one (Delwarde et al. 2007) or one per group of rows (Li et
  # al. 2009), and its value per row (a vector recycled along the rows)
  one <- is.null(groups)
  est_disp <- function(mu, start = NULL) {
    if (one) return(nb_dispersion(D, mu, start))
    vapply(seq_len(max(groups)), function(g) {
      nb_dispersion(D[groups == g, , drop = FALSE], mu[groups == g, , drop = FALSE], start[g])
    }, 0)
  }
  phi_g <- if (is.null(disp)) est_disp(fit_mu())
  phi <- disp %||% if (one) phi_g else phi_g[groups]
  if (!one) phi <- rep_len(phi, nrow(D))
  # score and observed information of the log-likelihood in the linear
  # predictor of each cell (the terms of eq. 6-8)
  score <- function(mu) {
    list(s = (D - mu) / (1 + phi * mu), w = mu * (1 + phi * D) / (1 + phi * mu)^2)
  }
  ll_old <- -Inf
  converged <- FALSE
  for (it in seq_len(maxit)) {
    if (!origin) {
      g <- score(fit_mu())
      a <- a + rowSums(g$s) / rowSums(g$w)
    }
    g <- score(fit_mu())
    kk <- matrix(k, nrow(D), nt, byrow = TRUE)
    b <- b + rowSums(g$s * kk) / rowSums(g$w * kk^2)
    g <- score(fit_mu())
    step <- colSums(g$s * b) / colSums(g$w * b^2)
    s <- sum(b)
    if (origin) {
      # k[T] stays 0 (Appendix 2); the scale of b and k is free
      step[nt] <- 0
      k <- (k + step) * s
    } else {
      k <- k + step
      a <- a + b * mean(k)
      k <- (k - mean(k)) * s
    }
    b <- b / s
    mu <- fit_mu()
    if (is.null(disp)) {
      phi_g <- est_disp(mu, phi_g)
      phi <- if (one) phi_g else phi_g[groups]
    }
    ll <- nb_loglik(D, mu, phi)
    if (!is.finite(ll)) mt_abort("The negative binomial Lee-Carter fit diverged.", class = "model_fit")
    if (abs(ll - ll_old) <= tol * abs(ll)) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }
  if (!converged) {
    mt_abort("The negative binomial Lee-Carter fit did not converge.", class = "model_fit")
  }
  out <- lc_finish(a, b, k, centre = !origin)
  out$method <- "nb"
  out$dispersion <- phi
  if (one) out$nb_size <- 1 / phi else out$age_nb_size <- 1 / phi
  out$loglik <- ll
  out$iterations <- it
  out
}

# Negative binomial log-likelihood (Delwarde, Denuit and Partrat 2007, eq. 5;
# Li, Hardy and Tan 2009, eq. 13) of deaths D with means mu and dispersion
# phi (variance mu + phi mu^2), in a form that is accurate for small phi and
# accepts non-integer deaths; phi = 0 gives the Poisson log-likelihood.
# lgamma(D + z) - lgamma(z) - lgamma(D + 1) = -lbeta(z, D) - log(D) for D >
# 0, with z = 1 / phi. `phi`: one dispersion, or one per cell (a vector of
# one per row of the matrix D is recycled along its rows).
nb_loglik <- function(D, mu, phi) {
  if (length(phi) > 1L) {
    phi <- rep_len(phi, length(D))
    pois <- phi <= 0
    out <- 0
    if (any(pois)) out <- nb_loglik(D[pois], mu[pois], 0)
    if (any(!pois)) out <- out + nb_loglik_cells(D[!pois], mu[!pois], phi[!pois])
    return(out)
  }
  if (phi <= 0) return(sum(ifelse(D > 0, D * log(mu), 0) - mu - lgamma(D + 1)))
  nb_loglik_cells(D, mu, phi)
}

# nb_loglik() for phi > 0 (one value, or one per cell)
nb_loglik_cells <- function(D, mu, phi) {
  phi <- rep_len(phi, length(D))
  z <- 1 / phi
  out <- -z * log1p(phi * mu)
  pos <- D > 0
  out[pos] <- out[pos] - lbeta(z[pos], D[pos]) - log(D[pos]) +
    D[pos] * (log(phi[pos] * mu[pos]) - log1p(phi[pos] * mu[pos]))
  sum(out)
}

# Maximum likelihood estimate of the dispersion phi >= 0 given the expected
# deaths mu. The derivative of the log-likelihood (eq. 5) at phi = 0 is
# sum((D - mu)^2 - D) / 2; where it is not positive, the estimate is 0
# (Poisson).
# Otherwise Newton steps on log(phi) with step halving (the paper: Newton
# steps on z = 1 / phi, eq. 9), from `start` or the moment estimate
# sum((D - mu)^2 - D) / sum(mu^2).
nb_dispersion <- function(D, mu, start = NULL, tol = 1e-10, maxit = 200) {
  s0 <- sum((D - mu)^2 - D)
  if (!is.finite(s0)) mt_abort("The dispersion cannot be estimated.", class = "model_fit")
  if (s0 <= 0) return(0)
  u <- log(if (isTRUE(start > 0)) start else max(s0 / sum(mu^2), 1e-10))
  ll <- nb_loglik(D, mu, exp(u))
  for (it in seq_len(maxit)) {
    z <- exp(-u)
    # derivatives of the log-likelihood in z (eq. 9) and in u = log(phi) = -log(z)
    g <- sum(digamma(D + z) - digamma(z) - log1p(mu / z) + (mu - D) / (z + mu))
    g1 <- sum(trigamma(D + z) - trigamma(z) + 1 / z - 1 / (z + mu) - (mu - D) / (z + mu)^2)
    d1 <- -z * g
    d2 <- z * g + z^2 * g1
    step <- if (is.finite(d2) && d2 < 0) -d1 / d2 else sign(d1)
    step <- max(min(step, 5), -5)
    for (half in 0:40) {
      ll_new <- nb_loglik(D, mu, exp(u + step))
      if (is.finite(ll_new) && ll_new >= ll) break
      step <- step / 2
    }
    if (!is.finite(ll_new) || ll_new < ll) break
    u <- u + step
    ll <- ll_new
    if (abs(step) < tol) break
  }
  exp(u)
}

# First year of the fitting period with the smallest ratio of the mean
# deviance of a linear-k fit to the mean deviance of the Lee-Carter fit
# (Booth, Maindonald and Smith 2002)
bms_start <- function(D, E, est_one, min_years) {
  nt <- ncol(D)
  na <- nrow(D)
  if (nt <= min_years + 1) return(1L)
  starts <- seq_len(nt - min_years)
  ratio <- vapply(starts, function(s) {
    cols <- s:nt
    Ds <- D[, cols, drop = FALSE]
    Es <- E[, cols, drop = FALSE]
    est <- tryCatch(est_one(Ds, Es), error = function(e) NULL)
    if (is.null(est)) return(Inf)
    n <- length(cols)
    # linear k through its mean with the mean annual change as slope
    klin <- mean(est$k) + mean(diff(est$k)) * (seq_len(n) - (n + 1) / 2)
    dev_base <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, est$k)))
    dev_lin <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, klin)))
    (dev_lin / (na * (n - 2))) / (dev_base / ((na - 1) * (n - 2)))
  }, 0)
  starts[which.min(ratio)]
}
