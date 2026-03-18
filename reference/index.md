# Package index

## All functions

- [`accuracy()`](https://kalleEF.github.io/mortensor/reference/accuracy.md)
  : Evaluate forecast accuracy
- [`add_pi()`](https://kalleEF.github.io/mortensor/reference/add_pi.md)
  : Add prediction intervals to a mort_forecast
- [`arima_forecast()`](https://kalleEF.github.io/mortensor/reference/arima_forecast.md)
  : Per-series ARIMA forecast
- [`as.data.table(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/as.data.table.mort_forecast.md)
  : Convert mort_forecast to a data.table
- [`blend_parametric_tail()`](https://kalleEF.github.io/mortensor/reference/blend_parametric_tail.md)
  : Blend GAM predictions with a parametric tail in log-space
- [`bootstrap_coherent_forecast()`](https://kalleEF.github.io/mortensor/reference/bootstrap_coherent_forecast.md)
  : Bootstrap PIs for Tucker4 coherent forecast
- [`bootstrap_forecast()`](https://kalleEF.github.io/mortensor/reference/bootstrap_forecast.md)
  : Bootstrap prediction intervals for Tucker3 forecasts
- [`bootstrap_gam_bayes()`](https://kalleEF.github.io/mortensor/reference/bootstrap_gam_bayes.md)
  : GAM Bayesian posterior propagation prediction intervals
- [`bootstrap_lc_fdm()`](https://kalleEF.github.io/mortensor/reference/bootstrap_lc_fdm.md)
  : Bootstrap PIs for Lee-Carter / FDM models
- [`bootstrap_stmomo()`](https://kalleEF.github.io/mortensor/reference/bootstrap_stmomo.md)
  : Bootstrap PIs for StMoMo-APC models
- [`build_array()`](https://kalleEF.github.io/mortensor/reference/build_array.md)
  : Build a 3D mortality array from a data.table
- [`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md)
  : Calibrate and compare mortality forecasting models
- [`calibrate_all_models()`](https://kalleEF.github.io/mortensor/reference/calibrate_all_models.md)
  : Calibrate all non-Tucker3 models
- [`calibrate_tucker3()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker3.md)
  : Calibrate Tucker3 via grid search
- [`calibrate_tucker4()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker4.md)
  : Calibrate Tucker4 via two-stage grid search
- [`check_disaggregation()`](https://kalleEF.github.io/mortensor/reference/check_disaggregation.md)
  : Check disaggregation quality
- [`check_smoothing()`](https://kalleEF.github.io/mortensor/reference/check_smoothing.md)
  : Check smoothing quality
- [`coef(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/coef.mort_forecast.md)
  : Extract model coefficients from a mort_forecast
- [`compare_accuracy()`](https://kalleEF.github.io/mortensor/reference/compare_accuracy.md)
  : Compare accuracy across multiple forecasts
- [`compare_coherent_vs_independent()`](https://kalleEF.github.io/mortensor/reference/compare_coherent_vs_independent.md)
  : Compare coherent Tucker4 vs independent per-disease Tucker3
- [`compare_models()`](https://kalleEF.github.io/mortensor/reference/compare_models.md)
  : Compare all candidate models on a train/test split
- [`compute_cause_shares()`](https://kalleEF.github.io/mortensor/reference/compute_cause_shares.md)
  : Compute cause shares from disease-specific and all-cause mortality
- [`compute_crps_empirical()`](https://kalleEF.github.io/mortensor/reference/compute_crps_empirical.md)
  : Empirical CRPS from bootstrap simulations (log scale)
- [`compute_crps_gaussian()`](https://kalleEF.github.io/mortensor/reference/compute_crps_gaussian.md)
  : CRPS for a Gaussian predictive distribution (closed-form)
- [`compute_dx()`](https://kalleEF.github.io/mortensor/reference/compute_dx.md)
  : Convert mortality rates to a death distribution (dx)
- [`compute_gam_ci_arrays()`](https://kalleEF.github.io/mortensor/reference/compute_gam_ci_arrays.md)
  : Compute pointwise Bayesian 95% credible intervals from GAM
  posteriors
- [`compute_interval_score()`](https://kalleEF.github.io/mortensor/reference/compute_interval_score.md)
  : Interval Score
- [`compute_interval_score_array()`](https://kalleEF.github.io/mortensor/reference/compute_interval_score_array.md)
  : Interval Score for 3D mortality arrays
- [`compute_log_score_gaussian()`](https://kalleEF.github.io/mortensor/reference/compute_log_score_gaussian.md)
  : Negative log-score for Gaussian predictive distribution
- [`compute_mae()`](https://kalleEF.github.io/mortensor/reference/compute_mae.md)
  : Mean Absolute Error
- [`compute_mase()`](https://kalleEF.github.io/mortensor/reference/compute_mase.md)
  : Mean Absolute Scaled Error (MASE) on log(mx)
- [`compute_mase_array()`](https://kalleEF.github.io/mortensor/reference/compute_mase_array.md)
  : Cell-level MASE across a 3D array
- [`compute_mase_combined()`](https://kalleEF.github.io/mortensor/reference/compute_mase_combined.md)
  : Combined MASE (cell + trend)
- [`compute_mase_group_trend()`](https://kalleEF.github.io/mortensor/reference/compute_mase_group_trend.md)
  : Group-trend MASE
- [`compute_rmse_log()`](https://kalleEF.github.io/mortensor/reference/compute_rmse_log.md)
  : RMSE on the log scale
- [`compute_wis()`](https://kalleEF.github.io/mortensor/reference/compute_wis.md)
  : Weighted Interval Score (WIS)
- [`correct_group_means()`](https://kalleEF.github.io/mortensor/reference/correct_group_means.md)
  : Smooth iterative group-mean correction
- [`cv_expanding_model()`](https://kalleEF.github.io/mortensor/reference/cv_expanding_model.md)
  : Expanding-window cross-validation for a generic model
- [`cv_rolling_origin()`](https://kalleEF.github.io/mortensor/reference/cv_rolling_origin.md)
  : Cross-validate Tucker3 with rolling origin
- [`cv_rolling_origin_tucker4()`](https://kalleEF.github.io/mortensor/reference/cv_rolling_origin_tucker4.md)
  : Rolling-origin CV for Tucker4 coherent forecast
- [`default_configs()`](https://kalleEF.github.io/mortensor/reference/default_configs.md)
  : Create default cause-specific configurations
- [`default_disease_cols()`](https://kalleEF.github.io/mortensor/reference/default_disease_cols.md)
  : Default disease columns (ONS UK data)
- [`default_disease_order()`](https://kalleEF.github.io/mortensor/reference/default_disease_order.md)
  : Default disease ordering for coherent ILR composition
- [`disaggregate()`](https://kalleEF.github.io/mortensor/reference/disaggregate.md)
  : Disaggregate a mort_forecast to single-year ages
- [`disaggregate_forecast_gam()`](https://kalleEF.github.io/mortensor/reference/disaggregate_forecast_gam.md)
  : Temporally-consistent GAM-based disaggregation with logistic blend
- [`disaggregate_hyman_kannisto()`](https://kalleEF.github.io/mortensor/reference/disaggregate_hyman_kannisto.md)
  : Disaggregate grouped mortality to single-year ages (Hyman +
  Kannisto)
- [`disaggregate_natural_gompertz()`](https://kalleEF.github.io/mortensor/reference/disaggregate_natural_gompertz.md)
  : Disaggregate grouped mortality to single-year ages (Natural +
  Gompertz)
- [`disaggregate_pi_ratios()`](https://kalleEF.github.io/mortensor/reference/disaggregate_pi_ratios.md)
  : Disaggregate prediction interval bounds using ratio interpolation
- [`dispatch`](https://kalleEF.github.io/mortensor/reference/dispatch.md)
  : Dispatch forecast or bootstrap for any supported model
- [`dispatch_forecast()`](https://kalleEF.github.io/mortensor/reference/dispatch_forecast.md)
  : Dispatch a full-horizon forecast for any model
- [`dispatch_model_sims()`](https://kalleEF.github.io/mortensor/reference/dispatch_model_sims.md)
  : Generate bootstrap simulations for any model
- [`dx_to_mx()`](https://kalleEF.github.io/mortensor/reference/dx_to_mx.md)
  : Convert a death distribution back to mortality rates
- [`enforce_coherence()`](https://kalleEF.github.io/mortensor/reference/enforce_coherence.md)
  : Enforce sex-differential coherence on forecasts
- [`ensemble_bootstrap()`](https://kalleEF.github.io/mortensor/reference/ensemble_bootstrap.md)
  : MASE-weighted ensemble bootstrap
- [`ets_forecast()`](https://kalleEF.github.io/mortensor/reference/ets_forecast.md)
  : Damped ETS forecast (per-series)
- [`example_mortality`](https://kalleEF.github.io/mortensor/reference/example_mortality.md)
  : Example stratified mortality data
- [`export()`](https://kalleEF.github.io/mortensor/reference/export.md)
  : Export mortality objects to files
- [`extrapolate_gompertz()`](https://kalleEF.github.io/mortensor/reference/extrapolate_gompertz.md)
  : Gompertz extrapolation: log-linear in age
- [`fdm_forecast()`](https://kalleEF.github.io/mortensor/reference/fdm_forecast.md)
  : Hyndman-Ullah Functional Demographic Model (FDM) forecast
- [`fitted(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/fitted.mort_forecast.md)
  : Extract fitted values from a mort_forecast
- [`forecast()`](https://kalleEF.github.io/mortensor/reference/forecast.md)
  : Forecast mortality from a calibrated model
- [`forecast(`*`<mort_calibration>`*`)`](https://kalleEF.github.io/mortensor/reference/forecast.mort_calibration.md)
  : Produce mortality forecasts from a calibrated model
- [`ilr_cause_transform()`](https://kalleEF.github.io/mortensor/reference/ilr_cause_transform.md)
  : ILR transform along the 4th mode of a shares array
- [`inv_clr_to_dx()`](https://kalleEF.github.io/mortensor/reference/inv_clr_to_dx.md)
  : Inverse CLR transform to death distribution
- [`lca_forecast()`](https://kalleEF.github.io/mortensor/reference/lca_forecast.md)
  : Lee-Carter forecast via demography package
- [`life_expectancy()`](https://kalleEF.github.io/mortensor/reference/life_expectancy.md)
  : Compute life expectancy
- [`life_table()`](https://kalleEF.github.io/mortensor/reference/life_table.md)
  : Compute a standard abridged life table
- [`methods`](https://kalleEF.github.io/mortensor/reference/methods.md)
  : S3 methods for mort_data and mort_forecast objects
- [`model_colors`](https://kalleEF.github.io/mortensor/reference/model_colors.md)
  : Standard colour palette for model comparison
- [`mort_calibration()`](https://kalleEF.github.io/mortensor/reference/mort_calibration.md)
  : Create a mortality calibration object
- [`mort_config()`](https://kalleEF.github.io/mortensor/reference/mort_config.md)
  : Create a validated configuration object
- [`mort_data()`](https://kalleEF.github.io/mortensor/reference/mort_data.md)
  : Create a mortality data object
- [`mort_forecast()`](https://kalleEF.github.io/mortensor/reference/mort_forecast.md)
  : Create a mortality forecast object
- [`mort_shared_config()`](https://kalleEF.github.io/mortensor/reference/mort_shared_config.md)
  : Create a shared configuration for multi-disease pipelines
- [`pipeline`](https://kalleEF.github.io/mortensor/reference/pipeline.md)
  : Master pipeline functions for cause-specific mortality forecasting
- [`plot(`*`<mort_calibration>`*`)`](https://kalleEF.github.io/mortensor/reference/plot.mort_calibration.md)
  : Plot a mort_calibration object
- [`plot(`*`<mort_data>`*`)`](https://kalleEF.github.io/mortensor/reference/plot.mort_data.md)
  : Plot a mort_data object
- [`plot(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/plot.mort_forecast.md)
  : Plot a mort_forecast object
- [`plot_B_components()`](https://kalleEF.github.io/mortensor/reference/plot_B_components.md)
  : Plot B-mode time components
- [`plot_consistency_check()`](https://kalleEF.github.io/mortensor/reference/plot_consistency_check.md)
  : Consistency check boxplot (sum of diseases / all-cause)
- [`plot_disagg_comparison()`](https://kalleEF.github.io/mortensor/reference/plot_disagg_comparison.md)
  : Disaggregation method comparison plot
- [`plot_disagg_correction()`](https://kalleEF.github.io/mortensor/reference/plot_disagg_correction.md)
  : Group-mean correction diagnostic
- [`plot_disagg_profile()`](https://kalleEF.github.io/mortensor/reference/plot_disagg_profile.md)
  : Disaggregation profile diagnostic
- [`plot_disagg_smoothness()`](https://kalleEF.github.io/mortensor/reference/plot_disagg_smoothness.md)
  : Disaggregation smoothness heatmap
- [`plot_disagg_temporal()`](https://kalleEF.github.io/mortensor/reference/plot_disagg_temporal.md)
  : Temporal consistency diagnostic
- [`plot_disaggregation()`](https://kalleEF.github.io/mortensor/reference/plot_disaggregation.md)
  : Plot disaggregated single-year mortality
- [`plot_forecast_age()`](https://kalleEF.github.io/mortensor/reference/plot_forecast_age.md)
  : Forecast age profile for a single year
- [`plot_forecast_trends()`](https://kalleEF.github.io/mortensor/reference/plot_forecast_trends.md)
  : Forecast trend plot (age-averaged, by group)
- [`plot_forecast_trends_by_age()`](https://kalleEF.github.io/mortensor/reference/plot_forecast_trends_by_age.md)
  : Forecast trends faceted by age group
- [`plot_model_comparison()`](https://kalleEF.github.io/mortensor/reference/plot_model_comparison.md)
  : Bar chart of model MASE scores
- [`plot_obs_vs_smooth()`](https://kalleEF.github.io/mortensor/reference/plot_obs_vs_smooth.md)
  : Observed vs smoothed comparison plot
- [`plot_share_evolution()`](https://kalleEF.github.io/mortensor/reference/plot_share_evolution.md)
  : Stacked area chart of cause shares (coherent mode)
- [`plot_smoothing_comparison()`](https://kalleEF.github.io/mortensor/reference/plot_smoothing_comparison.md)
  : Scatterplot of raw vs smoothed log(mx)
- [`plotting`](https://kalleEF.github.io/mortensor/reference/plotting.md)
  : Plotting functions for mortality forecasts
- [`read_mortality()`](https://kalleEF.github.io/mortensor/reference/read_mortality.md)
  : Read mortality data from a long-format table
- [`read_mortality_from_arrays()`](https://kalleEF.github.io/mortensor/reference/read_mortality_from_arrays.md)
  : Create a mort_data object from pre-built arrays
- [`read_ons_excel()`](https://kalleEF.github.io/mortensor/reference/read_ons_excel.md)
  : Read mortality data from ONS-format Excel files
- [`reconcile()`](https://kalleEF.github.io/mortensor/reference/reconcile.md)
  : Reconcile multiple cause-specific forecasts to sum to an all-cause
  forecast
- [`reconstruct_array()`](https://kalleEF.github.io/mortensor/reference/reconstruct_array.md)
  : Reconstruct a 3D array from Tucker3 factors
- [`reconstruct_array4()`](https://kalleEF.github.io/mortensor/reference/reconstruct_array4.md)
  : 4D Tucker reconstruction
- [`residuals(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/residuals.mort_forecast.md)
  : Extract residuals from a mort_forecast
- [`run_coherent_forecast()`](https://kalleEF.github.io/mortensor/reference/run_coherent_forecast.md)
  : Run the coherent multi-disease forecast pipeline (Tucker4-ILR)
- [`run_forecast()`](https://kalleEF.github.io/mortensor/reference/run_forecast.md)
  : Run the full mortality forecast pipeline for a single cause
- [`run_reconciled_forecast()`](https://kalleEF.github.io/mortensor/reference/run_reconciled_forecast.md)
  : Run the reconciled multi-disease forecast pipeline
- [`rwd_forecast()`](https://kalleEF.github.io/mortensor/reference/rwd_forecast.md)
  : Random Walk with Drift forecast (per-series)
- [`smooth_mortality()`](https://kalleEF.github.io/mortensor/reference/smooth_mortality.md)
  : Smooth mortality data
- [`stmomo_forecast()`](https://kalleEF.github.io/mortensor/reference/stmomo_forecast.md)
  : Stochastic Mortality Model (StMoMo) forecast
- [`` `[`( ``*`<mort_data>`*`)`](https://kalleEF.github.io/mortensor/reference/sub-.mort_data.md)
  : Subset a mort_data object by group indices
- [`summary(`*`<mort_forecast>`*`)`](https://kalleEF.github.io/mortensor/reference/summary.mort_forecast.md)
  : Summary of a mort_forecast object
- [`to_data_table()`](https://kalleEF.github.io/mortensor/reference/to_data_table.md)
  [`as.data.table(`*`<mort_data>`*`)`](https://kalleEF.github.io/mortensor/reference/to_data_table.md)
  : Convert a mort_data or mort_forecast object to a long-format
  data.table
- [`tucker3_coda_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_coda_forecast.md)
  : Tucker3-CoDa forecast
- [`tucker3_direct_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_direct_forecast.md)
  : Tucker3-Direct forecast
- [`tucker3_info_criterion()`](https://kalleEF.github.io/mortensor/reference/tucker3_info_criterion.md)
  : Information criterion for Tucker3 rank selection
- [`tucker4`](https://kalleEF.github.io/mortensor/reference/tucker4.md)
  : Tucker4-ILR coherent multi-disease forecasting
- [`tucker4_coherent_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker4_coherent_forecast.md)
  : Tucker4 forecast on ILR cause-shares
- [`wls_reconcile()`](https://kalleEF.github.io/mortensor/reference/wls_reconcile.md)
  : WLS reconciliation of independent disease forecasts
- [`workflow`](https://kalleEF.github.io/mortensor/reference/workflow.md)
  : Step-by-step mortality forecasting workflow
