# Function map

This page maps the functions of mortensor onto the steps of an analysis
and shows how their results feed into each other. Each box names a
function and its main purpose. Paths through the maps are workflows; the
section “What is not covered” lists steps that need other tools.

## The main path

Every analysis passes the same objects from left to right. Functions in
the boxes create the objects; the labels on the arrows are the functions
that connect them.

``` mermaid
flowchart TD
  T["`**your tables**
  deaths, populations or exposures (CSV, Excel, official exports)`"]
  D["`**mt_data object**
  observations at any age and period resolution, strata`"]
  S["`**model specification**
  lee_carter(), li_lee(), rwd(), ... and wrappers ensemble(), reconciled()`"]
  F["`**mt_fit**
  estimates per stratum or joint`"]
  P["`**mt_draws**
  simulated paths of future rates`"]
  B["`**mt_backtest**
  scores of past forecasts`"]
  R["`**tables and plots**
  life expectancy, rates, deaths, gaps, inequality indices`"]
  T -- "mt_data()" --> D
  D -- "prepare and check (map 1)" --> D
  D -- "mt_fit(spec, data, prep)" --> F
  S --> F
  F -- "mt_forecast()" --> P
  P -- "map 3" --> R
  D -- "mt_backtest(spec, data, prep)" --> B
  B -- "mt_score()" --> R
  D -- "mt_observed(): the past" --> R
```

## Map 1: from tables to analysis-ready data

Each step is an analysis decision. Most data sets need only some of
them.

``` mermaid
flowchart TD
  subgraph READ["Read"]
    A1["`**mt_data()**
    deaths, exposures, populations in long
    or wide tables; age labels such as 85+,
    '85 years and over', '85 Jahre und mehr';
    suppression symbols; totals and unknown
    ages handled with warnings`"]
    A2["`**mt_read_genesis(), mt_read_gbe()**
    German official exports`"]
    A3["`**mt_from_hmd(), as_mt_data()**
    Human Mortality Database, StMoMo,
    demography objects`"]
  end
  subgraph PREP["Prepare (decisions)"]
    B1["`**mt_as_exposure(type)**
    person-years from mid-year,
    1 January or 31 December populations`"]
    B2["`**mt_subset()**
    ages, years, strata to analyse`"]
    B3["`**mt_group_levels(), mt_relabel()**
    aggregation layers: cause groups,
    quintiles, regions; recode labels`"]
    B4["`**mt_fill_suppressed()**
    estimate suppressed counts from
    published totals or known bounds`"]
  end
  subgraph CHECK["Check"]
    C1["`**mt_validate()**
    impossible values, missing cells,
    totals among levels, sums of a hierarchy`"]
  end
  subgraph HARM["Harmonise resolution"]
    H1["`**mt_regroup()**
    common age groups or years
    ('deaths', 'common')`"]
    H2["`**mt_ungroup()**
    single years from age groups
    (penalised composite link)`"]
    H3["`**mt_smooth()**
    smooth rates over age, time or both
    (P-splines, Whittaker, GAM);
    mt_smooth_check(), mt_plot_smooth()`"]
    H4["`**mt_prep() with step_*()**
    the same steps inside fits
    and backtests, without leakage`"]
  end
  A1 --> B1
  A2 --> B1
  A3 --> B2
  B1 --> B2 --> B3 --> B4 --> C1
  C1 --> H1
  C1 --> H2
  H2 --> H3
  H1 --> M["model-ready mt_data"]
  H3 --> M
  H4 -.-> M
```

## Map 2: models, coherence and evaluation

``` mermaid
flowchart TD
  M["model-ready mt_data"]
  subgraph ONE["Models for one population per stratum"]
    direction LR
    S1["`**lee_carter()**
    age pattern of a common trend`"]
    S2["`**rwd()**
    random walk per age; benchmark;
    observation_noise for small counts`"]
    S3["`**hyndman_ullah(), stmomo(),
    arima_series()**
    demography, StMoMo, forecast`"]
    S4["`**apc_bayes(), apc_select()**
    Bayesian age-period-cohort`"]
  end
  subgraph MANY["Models for several strata jointly"]
    direction LR
    J1["`**li_lee(common_by)**
    common trend, bounded deviations`"]
    J2["`**product_ratio(), tucker3(),
    cause_shares()**
    other joint and cause models`"]
  end
  subgraph WRAP["Wrappers around any model"]
    direction LR
    W1["`**ensemble()**
    mixture of models,
    better calibrated`"]
    W2["`**reconciled(model, hierarchy)**
    totals and parts add up
    draw by draw`"]
    W3["`**bootstrapped(), tuned()**
    parameter uncertainty,
    inner tuning`"]
    W4["`**calibrated(model, years = 10)**
    spread of the draws fitted to
    the errors of the last years`"]
    W5["`**anchored(model, years = 3)**
    start of the forecast moved to
    the observed level`"]
  end
  H["`**mt_hierarchy(), mt_aggregate_strata()**
  totals and groups of strata; causes share exposures;
  fixed nodes (a trusted total) in reconciliation`"]
  FIT["`**mt_fit()**`"]
  INS["`**summary(), mt_parameters(), residuals(), mt_plot_residuals()**
  inspect the model, its assumptions and the fit`"]
  FC["`**mt_forecast()**
  draws; warns about implausible rates;
  mt_jump_off() reports the first-year step`"]
  M --> FIT
  H --> W2
  ONE --> WRAP
  MANY --> WRAP
  ONE --> FIT
  MANY --> FIT
  WRAP --> FIT
  FIT --> INS
  FIT --> FC
```

Choosing between specifications:

``` mermaid
flowchart TD
  M["model-ready mt_data"]
  SP["`**candidate specifications**
  with their preprocessing (mt_prep)`"]
  BT["`**mt_backtest()**
  rolling origins; preprocessing and model refitted
  on each training window; suppressed cells not scored`"]
  SC["`**mt_score()**
  CRPS, interval score, coverage, PIT, MASE;
  for rates, life expectancy or ASR`"]
  PL["`**mt_plot_scores(), mt_plot_coverage()**
  by horizon, age or model`"]
  CH["chosen specification"]
  M --> BT
  SP --> BT
  BT --> SC --> CH
  BT --> PL --> CH
```

## Map 3: from draws to results

``` mermaid
flowchart TD
  P["`**mt_draws**
  forecast paths`"]
  O["`**mt_observed(data)**
  observed rates with Poisson draws`"]
  subgraph IND["Indicators, path by path"]
    direction LR
    I1["`**life_expectancy()**`"]
    I2["`**asr()**
    standards esp2013, esp1976,
    who2000 or your own`"]
    I3["`**prob_dying()**`"]
    I4["`**mt_deaths(), yll()**
    deaths with projected
    populations`"]
  end
  subgraph COMP["Comparisons, path by path"]
    direction LR
    C1["`**mt_contrast()**
    differences and ratios
    between strata`"]
    C2["`**mt_inequality()**
    slope and relative index
    of inequality`"]
    C3["`**mt_decompose()**
    differences by age, cause
    and subgroup composition`"]
  end
  subgraph AGG["Aggregation"]
    direction LR
    G1["`**mt_aggregate_draws()**
    totals, wider age groups;
    mt_coherence() checks sums`"]
    G2["`**mt_disaggregate()**
    finer ages (mt_shape_rates(),
    mt_disaggregate_check())`"]
  end
  subgraph OUT["Output"]
    direction LR
    X1["`**mt_summarise()**
    quantiles of rates`"]
    X2["`**mt_plot_forecast(), mt_plot_age(),
    mt_plot_life_expectancy(), mt_plot_shares()**`"]
    X3["`**mt_export()**
    CSV of draws or any table`"]
  end
  P --> AGG
  P --> IND
  O --> IND
  P --> COMP
  O --> COMP
  IND --> OUT
  COMP --> OUT
  P --> OUT
```

## Functions by purpose

| Step | Function | Purpose | Input | Output |
|----|----|----|----|----|
| Read | [`mt_data()`](../reference/mt_data.md) | read deaths, exposures, populations (long or wide, suppression symbols, shared populations) | data frames | `mt_data` |
| Read | [`mt_add()`](../reference/mt_add.md) | add or replace one table | `mt_data`, data frame | `mt_data` |
| Read | [`mt_read_genesis()`](../reference/mt_read_genesis.md), [`mt_read_gbe()`](../reference/mt_read_gbe.md) | read German official exports | files | `mt_data`, table |
| Read | [`mt_from_hmd()`](../reference/mt_from_hmd.md), [`mt_read_hmd()`](../reference/mt_from_hmd.md), [`as_mt_data()`](../reference/as_mt_data.md) | read or convert HMD, StMoMo, demography data | objects, files | `mt_data` |
| Read | [`parse_age_labels()`](../reference/parse_age_labels.md), [`age_groups()`](../reference/mt_axis.md), [`single_ages()`](../reference/mt_axis.md), [`single_years()`](../reference/mt_axis.md) | age and period axes from labels or bounds | labels | `mt_axis` |
| Prepare | [`mt_as_exposure()`](../reference/mt_as_exposure.md) | exposures from populations (`type` required) | `mt_data` | `mt_data` |
| Prepare | [`mt_subset()`](../reference/mt_subset.md) | restrict ages, years, strata | `mt_data` | `mt_data` |
| Prepare | [`mt_strata()`](../reference/mt_strata.md), [`mt_obs()`](../reference/mt_strata.md), [`mt_grid()`](../reference/mt_strata.md), [`mt_set_grid()`](../reference/mt_set_grid.md) | look at strata (with their level in a hierarchy), tables and the modelling grid | `mt_data` | table, grid |
| Prepare | [`mt_group_levels()`](../reference/mt_group_levels.md), [`mt_relabel()`](../reference/mt_relabel.md) | aggregation layers, recode labels | `mt_data` | `mt_data` |
| Prepare | [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) | estimate suppressed counts | `mt_data`, hierarchy or bounds | `mt_data` |
| Check | [`mt_validate()`](../reference/mt_validate.md) | plausibility and hierarchy checks | `mt_data` | findings |
| Harmonise | [`mt_regroup()`](../reference/mt_regroup.md), [`mt_ungroup()`](../reference/mt_ungroup.md), [`mt_smooth()`](../reference/mt_smooth.md) | common resolution, single years, smoothing over age, time or both | `mt_data` | `mt_data` |
| Harmonise | [`mt_smooth_check()`](../reference/mt_smooth_check.md), [`mt_plot_smooth()`](../reference/mt_smooth_check.md) | how much smoothing removed (dispersion, trend ratio) | `mt_data` pair | table, ggplot |
| Harmonise | [`mt_prep()`](../reference/mt_prep.md) with [`step_ungroup()`](../reference/mt_prep.md), [`step_smooth()`](../reference/mt_prep.md), [`step_regroup()`](../reference/mt_prep.md), [`step_grid()`](../reference/mt_prep.md), [`step_shocks()`](../reference/mt_prep.md), [`step_window()`](../reference/mt_prep.md) | preprocessing inside fits and backtests | steps | pipeline |
| Model | [`lee_carter()`](../reference/lee_carter.md), [`rwd()`](../reference/rwd.md), [`li_lee()`](../reference/li_lee.md), [`product_ratio()`](../reference/li_lee.md), [`hyndman_ullah()`](../reference/hyndman_ullah.md), [`stmomo()`](../reference/stmomo.md), [`arima_series()`](../reference/arima_series.md), [`ets_series()`](../reference/arima_series.md), [`apc_bayes()`](../reference/apc_bayes.md), [`tucker3()`](../reference/tucker3.md), [`cause_shares()`](../reference/cause_shares.md) | model specifications | settings | `mt_spec` |
| Model | [`apc_select()`](../reference/apc_select.md), [`mcmc_control()`](../reference/mcmc_control.md), [`mt_tune()`](../reference/mt_tune.md) | choose random-walk orders, MCMC settings, settings by inner backtest | data, grid | table, `mt_spec` |
| Combine | [`ensemble()`](../reference/ensemble.md) (optionally `top_k`), [`reconciled()`](../reference/reconciled.md), [`bootstrapped()`](../reference/bootstrapped.md), [`tuned()`](../reference/mt_tune.md) | mixtures, coherence, parameter uncertainty (and of smoothing, `type = "data"`), tuning | `mt_spec` | `mt_spec` |
| Combine | [`mt_order_strata()`](../reference/mt_order_strata.md) | keep an ordering between strata (for example women below men) path by path | `mt_draws` | `mt_draws` |
| Calibrate | [`calibrated()`](../reference/calibrated.md) | rescale the spread (optionally the drift) of any model’s draws to its forecast errors in the last years | `mt_spec` | `mt_spec` |
| Jump-off | [`anchored()`](../reference/anchored.md), [`mt_jump_off()`](../reference/mt_jump_off.md) | start forecasts at the observed level of the last years; report the first-year step | `mt_spec`, `mt_draws` | `mt_spec`, table |
| Coherence | [`mt_hierarchy()`](../reference/mt_hierarchy.md), [`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md), [`mt_reconcile()`](../reference/mt_reconcile.md) (with `fixed` nodes; `"proportional"` scales the parts of fixed nodes), [`mt_aggregate_draws()`](../reference/mt_aggregate_draws.md) | define and apply totals of strata | hierarchy | `mt_data`, `mt_draws` |
| Coherence | [`mt_coherence()`](../reference/mt_coherence.md), [`mt_plot_coherence()`](../reference/mt_coherence.md) | sum of parts over aggregate, per year | `mt_draws`, hierarchy | table, ggplot |
| Fit | [`mt_fit()`](../reference/mt_fit.md) | estimate | `mt_spec`, `mt_data` | `mt_fit` |
| Inspect | [`summary()`](https://rdrr.io/r/base/summary.html), [`mt_parameters()`](../reference/mt_parameters.md), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`mt_plot_residuals()`](../reference/mt_plot_residuals.md) | model, assumptions, estimates, fit | `mt_fit` | text, table, plot |
| Forecast | [`mt_forecast()`](../reference/mt_forecast.md), [`mt_summarise()`](../reference/mt_draws.md) | simulated paths and their quantiles | `mt_fit` | `mt_draws`, table |
| Evaluate | [`mt_backtest()`](../reference/mt_backtest.md), [`mt_score()`](../reference/mt_score.md) | rolling-origin scores of rates, life expectancy or ASR | `mt_spec`, `mt_data` | table |
| Describe | [`mt_observed()`](../reference/mt_observed.md) | observed rates as draws | `mt_data` | `mt_draws` |
| Report | [`life_expectancy()`](../reference/life_expectancy.md), [`asr()`](../reference/derived.md), [`prob_dying()`](../reference/derived.md) | indicators path by path, observed or forecast | `mt_draws` or `mt_data` | table |
| Report | [`yll()`](../reference/derived.md), [`mt_deaths()`](../reference/mt_deaths.md) | years of life lost, deaths with projected populations | `mt_draws` | table, array |
| Report | [`life_table()`](../reference/life_table.md), [`standard_population()`](../reference/standard_population.md) | life tables from rates; standard populations esp2013, esp1976, who2000 | rates, name | table |
| Report | [`mt_contrast()`](../reference/mt_contrast.md), [`mt_inequality()`](../reference/mt_inequality.md) | gaps, ratios, SII, RII | `mt_draws` or `mt_data` | table |
| Report | [`mt_decompose()`](../reference/mt_decompose.md), [`mt_plot_decomposition()`](../reference/mt_plot_decomposition.md) | differences in life expectancy (or other summaries) by age, cause and subgroup composition: Arriaga, Horiuchi, stepwise replacement, Kitagawa, Arriaga meets Kitagawa; jump-off share and changes of gaps | `mt_draws` or `mt_data` (or both; also census years with gaps) | `mt_decomposition`, table, ggplot |
| Report | [`mt_plot_forecast()`](../reference/mt_plot.md), [`mt_plot_age()`](../reference/mt_plot.md), [`mt_plot_life_expectancy()`](../reference/mt_plot.md), [`mt_plot_shares()`](../reference/mt_plot_more.md), [`mt_plot_rates()`](../reference/mt_plot.md), [`mt_plot_ungroup()`](../reference/mt_plot.md), [`mt_plot_tucker()`](../reference/mt_plot_more.md), [`mt_export()`](../reference/mt_export.md) | figures and CSV files | draws, data, tables | ggplot, file |
| Data | `ons_data`, `ons_cause_groups`, `de_deaths`, `de_population`, `de_cause_groups`, [`mt_example()`](../reference/mt_example.md), [`mt_example_causes()`](../reference/mt_example_causes.md), [`mt_example_regions()`](../reference/mt_example_regions.md) | example data (published and simulated) |  | tables, `mt_data` |
| Report | [`theme_mortensor()`](../reference/mt_plot_style.md), [`scale_colour_mortensor()`](../reference/mt_plot_style.md), [`mt_palette()`](../reference/mt_plot_style.md) | the package’s plot style for your own figures |  | theme, scales |
| Disaggregate | [`mt_disaggregate()`](../reference/mt_disaggregate.md), [`mt_shape_rates()`](../reference/mt_shape_rates.md), [`mt_disaggregate_check()`](../reference/mt_disaggregate_check.md), [`mt_plot_disaggregation()`](../reference/mt_disaggregate_check.md) | forecasts on finer ages, with fixed or changing age shapes | `mt_draws` | `mt_draws`, table, ggplot |
| Exchange | [`as_demogdata()`](../reference/as_demogdata.md), [`as_stmomodata()`](../reference/as_demogdata.md) | cross-checks with demography and StMoMo | `mt_data` | `demogdata`, `StMoMoData` |

## What is not covered

These steps of a mortality analysis need other tools or are open
development points:

- **Population projections.** mortensor forecasts death rates; future
  exposures (for expected deaths) come from a population projection, for
  example of the statistical office, passed as `exposure`.
- **Spatial models for small areas.** There is no model with spatial
  smoothing between neighbouring regions (such as BYM models); pooling
  works through common trends ([`li_lee()`](../reference/li_lee.md));
  Bayesian APC models smooth over age, period and cohort within each
  stratum.
- **Uncertainty of score differences.**
  [`mt_score()`](../reference/mt_score.md) reports mean scores; standard
  errors over origins and strata or paired tests of two models are not
  computed.
- **Lifespan variation.** Lifespan variation (e-dagger, Gini) is not a
  built-in quantity; [`mt_decompose()`](../reference/mt_decompose.md)
  accepts it as a function of the rates.
- **Multiple imputation of suppressed cells.**
  [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) returns
  one estimate; the uncertainty of suppressed counts is not carried into
  forecasts.
- **Causes of death with changing coding.** ICD revisions and coding
  changes (for example of dementia) are not bridged; cause groups that
  span the affected codes reduce the problem.
