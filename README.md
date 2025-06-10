# OlinkAnalysisHEDIMED
An R-package to speed up analysing Olink data - a custom companion package to supplement the official OlinkAnalyze package.

Developed by Jonathan Thorsen 2023-25

## List of functions

These functions assume that Olink data is supplied in long format with relevant covariates joined.
Some functions expect target96 immune protein names. Feel free to fork and adapt to other datasets.
Olink-related column names are assumed to be lowecase.

### Utilities

show_outcome_distribution(olink_data, outcome, annotation = NULL)
- Designed to count outcome prevalence in different datasets

make_pca(npx_data, getWide = FALSE)
- Quickly make a pca from an olink dataset

get_calibrated_residuals(formula, data)
- Do a batch calibration using supplied covariates

### Differential expression

olink_lm(olink_data, outcome, covariates = "1", annotation = "")
- Do a linear model with NPX values on left hand side and outcome (and any covariates) on right hand side. Optionally annotate with model description, useful if aggregating may such calls.

olink_lmer(olink_data, outcome, covariates = "1", annotation = "")
- The same as above, but with a linear mixed model using lmer. See lme4 documentation for how to do random effects.

olink_volcano_jt(olink_res)
- Visualize the results from olink_lm or olink_lmer as a ggplot2 based volcano-plot.

olink_lm_corrs(reslist, method = "spearman")
- Correlate npx estimates between two or more lm / lmer runs

olink_corrs_heatmap(corrdata, xtextangle = 30)
- Visualize such a correlation test of estimates between lm runs in a heatmap

### Supervised analysis

olink_fit_sPLS(olink_data, outcome, seed = 123, fixX = c(), annotation = "sPLS analysis")
- Run a sparse Partial Least Squares model from the mixOmics package using repeated cross-validation from the caret-package.

make_cross_pls_list(pls_objects = list())
- Supplying a list of sPLS objects from olink_fit_sPLS, do a cross-cohort prediction of all vs all and get results back as a list

make_cross_pls_list_plotdat(cross_pls_list = list())
- Convert list from make_cross_pls_list to a data.frame suitable for plotting

plot_results_of_pls_models(pls_model_list = list())
- Plot a summary of resuls from a list of sPLS objects from olink_fit_sPLS

plot_cross_auc_heatmap(cross_results = data.frame())
- Plot a heatmap with AUC values from cross-pls comparison results, using output from make_cross_pls_list_plotdat


