# Decomposition as post-processing: plan

Status: implemented in 0.2 (work package Q1); revised after applied reviews 3 and 4 (work package Q2: totals refused, `composition` required for subgroups, `observed = "fitted"`, jump-off, change of a gap; see `dev/progress.md`).

## Aim

Explain differences in a mortality summary (life expectancy, probability of
dying, age-standardised or crude rate) by the ages, causes and population
subgroups that produce them, for

* observed data (`mt_data`, optionally with Poisson intervals),
* forecasts (`mt_draws`, draw by draw, so that every contribution has an
  interval), and
* both at once (an observed year against a forecast year of the same draws,
  or observed data against a forecast object).

Decomposition is a post-processing step: it takes rates that the package
already produces and does not change models or draws.

## Methods

| Method | Quantity | Components | Exact | Reference |
|---|---|---|---|---|
| `"arriaga"` | life expectancy | age (cause by implied sensitivity) | yes | Arriaga 1984; symmetric form and cause split as in Riffe et al. |
| `"horiuchi"` | any function of rates | age, cause, subgroup rates and composition | up to the number of steps | Horiuchi, Wilmoth and Pletcher 2008 |
| `"stepwise"` | any function of rates | as above | yes (order-dependent; symmetric average) | Andreev, Shkolnikov and Begun 2002 |
| `"kitagawa"` | crude death rate | rate and age (or age-by-subgroup) composition | yes | Kitagawa 1955 |
| subgroups, `composition = "radix"` | life expectancy, probability of dying | subgroup rates by age and cause, one composition term | yes with `"arriaga"`/`"stepwise"` | Riffe, Tursun-zade and Trias Llimós 2024 |
| subgroups, `composition = "pooled"` | any function of the pooled rates | subgroup rates by age and cause, composition by age | as the age method | Horiuchi et al. 2008 applied to `m_a = sum_s pi_a^s m_a^s` (Vaupel 2002, "current rates") |

### Arriaga meets Kitagawa (Riffe, Tursun-zade and Trias Llimós 2024)

Total life expectancy at the radix age `a0` is the mixture of subgroup life
tables blended at that age, `e = sum_s pi^s e^s` (eq. 2), with `pi^s` the
population share of subgroup `s` at `a0`. For two populations 1 and 2:

* composition: `Delta_pi = sum_s (pi^{s,2} - pi^{s,1}) * mean(e^s)` (eq. 3);
  its subgroup terms are not identified and are reported as one total;
* rate effect of subgroup `s`: `Delta_e^s = (e^{s,2} - e^{s,1}) * mean(pi^s)`
  (eq. 5);
* the age (and cause) decomposition of `e^{s,2} - e^{s,1}` by any exact
  method, rescaled by `mean(pi^s)` (eq. 14, which reduces to this weight).

Symmetric Arriaga (eqs. 8-11): the average of the forward decomposition
(1 to 2) and the negated backward one (2 to 1). Causes (eqs. 12-13): the
implied all-cause sensitivity `s_a = Delta_a / (m_a^2 - m_a^1)` times the
cause-specific rate differences; where the all-cause difference is close to
zero, the sensitivity is taken from a numerical derivative of life
expectancy at the average rates (the robust route the paper recommends).

The pooled alternative (`composition = "pooled"`) decomposes the life
expectancy of the pooled rates that national offices publish (and that the
aggregate strata of a hierarchy hold) into subgroup rate effects and
age-specific composition effects with the line integral.

## Interface

```r
mt_decompose(x, y = NULL, var = NULL, ref = NULL, periods = NULL,
             quantity = "life_expectancy", method = "arriaga",
             cause = NULL, groups = NULL, hierarchy = NULL,
             composition = c("radix", "pooled"), at = NULL, to = NULL,
             standard = "esp2013", exposure = NULL, n_steps = 20,
             n_draws = 1, probs = c(0.025, 0.5, 0.975), seed = 1)
summary(dec, by = c("age", "cause"))   # sums path by path, then quantiles
mt_plot_decomposition(dec, by = "age", fill = "cause")
```

Comparisons:

1. `var` and `ref`: each level of a stratum variable against a reference
   level (for example quintile Q5 against Q1), matched on the other stratum
   variables, in each of `periods` (like `mt_contrast()`);
2. `y`: `y` minus `x` for matching strata and periods (two models, or
   observed data against a forecast);
3. neither: change from the first of `periods` to each later one, per
   stratum; with draws, `periods` may mix observed and forecast years.

`cause` and `groups` name stratum variables whose levels are the causes and
the subgroups; their total levels (from `hierarchy`, if given) are dropped.
`quantity` may be a function `f(m, age)` of a rate matrix `[age, n]`
returning `n` values, for `"horiuchi"` and `"stepwise"`.

Uncertainty: every comparison is computed draw by draw. Observed data give a
single path (`n_draws = 1`) or Poisson draws (`n_draws > 1`, as
`mt_observed()`). Observed years inside forecast draws are fixed across
draws. Sums over ages, causes or subgroups are formed path by path before
quantiles are taken.

## Validation

* Exactness: contributions add up to the difference (`arriaga`, `stepwise`,
  `kitagawa` to machine precision, `horiuchi` to the step error).
* Identity: forward Arriaga equals stepwise replacement from the oldest to
  the youngest age, so the symmetric versions agree exactly for life
  expectancy.
* Reference implementations (tests, `skip_if_not_installed()`): DemoDecomp
  (`horiuchi()`, `stepwise_replacement()`) and LEdecomp (`arriaga_sym()`) on
  the same rates and life table.
* The Kitagawa step reproduces Table 1 of Riffe et al. (2024) from its
  rounded inputs.
* Validation check `decomposition`: exactness and agreement of the methods
  on ONS and HMD data.

## Documentation

* Article "Decomposing differences in life expectancy" on `ons_data`: the
  deprivation gap by age and cause (observed), the change 2001-2019 of
  England with quintile composition (Arriaga-Kitagawa), the forecast change
  2019-2035 by age with intervals, and two models compared.
* Function map, README table, FAQ entry, pkgdown reference section.

## Later

* Lifespan variation (e-dagger) and other inequality measures as built-in
  quantities (now possible through a function).
* Continuous change over time (Vaupel and Canudas-Romo 2003), Das Gupta's
  multi-factor standardisation, cohort decompositions.

## References (checked)

* Andreev EM, Shkolnikov VM, Begun AZ (2002). Algorithm for decomposition of
  differences between aggregate demographic measures and its application to
  life expectancies, healthy life expectancies, parity-progression ratios
  and total fertility rates. *Demographic Research* 7: 499-522.
* Arriaga EE (1984). Measuring and explaining the change in life
  expectancies. *Demography* 21(1): 83-96. doi:10.2307/2061029
* Horiuchi S, Wilmoth JR, Pletcher SD (2008). A decomposition method based
  on a model of continuous change. *Demography* 45(4): 785-801.
  doi:10.1353/dem.0.0033
* Kitagawa EM (1955). Components of a difference between two rates.
  *Journal of the American Statistical Association* 50(272): 1168-1194.
  doi:10.1080/01621459.1955.10501299
* Riffe T, Tursun-zade R, Trias Llimós S (2024). Arriaga meets Kitagawa.
  Life expectancy decomposition with population subgroups. MPIDR Working
  Paper WP-2024-029. doi:10.4054/MPIDR-WP-2024-029. Published in
  *Geographical Analysis* (2026), doi:10.1111/gean.70033.
* Vaupel JW (2002). Life expectancy at current rates vs. current
  conditions. *Demographic Research* 7: 365-378.
