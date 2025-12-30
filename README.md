# fQGAM: Functional Quantile Regression via Generalized Additive Models

**fQGAM** implements functional quantile regression for longitudinal functional data using generalized additive models. The package provides:

- Tools for fitting quantile regression models with functional covariates via `qgam`
- **Block bootstrap** for variance estimation (resamples complete subject trajectories)
- **Wild bootstrap** for bias adjustment (perturbs residuals and resamples random effects)
- Helper functions for data preparation, prediction, and inference


The methodology is described in two papers:

> Battagliola, M.L., Sørensen, H., Tolver, A., & Staicu, A.-M. (2024). *Quantile regression for longitudinal functional data with application to feed intake of lactating sows*. Journal of Agricultural, Biological, and Environmental Statistics, 30(1), 211-230. [doi:10.1007/s13253-024-00601-5](https://doi.org/10.1007/s13253-024-00601-5)

> Battagliola, M.L., Sørensen, H., Tolver, A., & Staicu, A.-M. (2022). *A bias-adjusted estimator in quantile regression for clustered data*. Econometrics and Statistics, 23, 165-186. [doi:10.1016/j.ecosta.2021.07.003](https://doi.org/10.1016/j.ecosta.2021.07.003)

## Overview

Given longitudinal functional data $\{(Y_{ij}, X_{ij}(\cdot), t_{ij})\}$ for subjects $i = 1, \ldots, N$ with $j = 1, \ldots, n_i$ observations each, the functional quantile regression model is:

$$Q_{Y_{ij}| X_{ij}, u_i}^\tau(t_{ij}) = \alpha^\tau(t_{ij}) + \int_{\mathcal{S}} \beta^\tau(s,t_{ij}) X_{ij}(s)ds + u_i$$

where:
- $\tau \in (0,1)$ is the quantile level
- $\alpha^\tau(t)$ is the time-varying intercept
- $\beta^\tau(s,t)$ is the bivariate functional coefficient
- $u_i$ are subject-specific random effects with $E[u_i] = 0$


For estimation and inference on the coefficients of the model, the recommended approach (from Battagliola et al. 2024) combines two bootstrap schemes:

### Block Bootstrap (`boot_pred_block`)
- **Purpose**: Variance estimation and standard error computation
- **Method**: Resamples complete subject trajectories with replacement
- **Use case**: Constructing confidence intervals

### Wild Bootstrap (`boot_pred_wild`)
- **Purpose**: Bias adjustment
- **Method**: Combines resampling of estimated random effects with wild bootstrap of residuals
- **Key feature**: Generates bootstrap data where the true parameter equals the original estimate, allowing bias measurement



## Installation and quick start:

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("LauraBattagliola/fQGAM")
```

To run the package, you will need \texttt{R} packages \texttt{qgam} and \texttt{mgcv}.

Example:

```r
library(fQGAM)
library(refund)
library(qgam)

# Load and prepare DTI data
data(DTI)
myData <- DTI[complete.cases(DTI$cca), ]
myData <- subset(myData, Nscans > 1)

# Prepare data for fQGAM
myDataFit <- prep_fqgam_data(
  y = myData$pasat,
  id = myData$ID,
  time = myData$visit.time,
  Xobs = myData$cca
)

# Fit the model at tau = 0.1
tau <- 0.1
formula <- y ~ s(time, bs = "cr", k = 10) + 
              te(sMat, tMat, by = Xobs, bs = c("cr", "cr"), k = c(10, 10)) + 
              s(id, bs = "re")

fit <- qgam(formula, qu = tau, data = myDataFit)

# Get reference curves (20% and 80% pointwise quantiles)
cca20 <- apply(myData$cca, 2, quantile, probs = 0.2)
cca80 <- apply(myData$cca, 2, quantile, probs = 0.8)

# Run bootstrap methods (use B >= 100 in practice)
bb <- boot_pred_block(myDataFit, B = 100, seed = 1234, 
                      tau = tau, X20 = cca20, X80 = cca80, model = fit)

wb <- boot_pred_wild(myDataFit, B = 100, seed = 4321,
                     tau = tau, X20 = cca20, X80 = cca80, model = fit)

# Compute original predictions
model_results <- compute_model_se(fit, cca20, cca80)
predDiff <- model_results$predDiff

# Compute bias-adjusted confidence intervals
ci <- boot_ci(block_boot = bb, wild_boot = wb, 
              original_pred = predDiff, alpha = 0.05)

print(ci)
```


## Citation

If you use this package in your research, please cite both papers:

```bibtex
@article{battagliola2024quantile,
  title={Quantile regression for longitudinal functional data with application to feed intake of lactating sows},
  author={Battagliola, Maria Laura and S{\o}rensen, Helle and Tolver, Anders and Staicu, Ana-Maria},
  journal={Journal of Agricultural, Biological, and Environmental Statistics},
  volume={30},
  number={1},
  pages={211--230},
  year={2024},
  doi={10.1007/s13253-024-00601-5}
}

@article{battagliola2022bias,
  title={A bias-adjusted estimator in quantile regression for clustered data},
  author={Battagliola, Maria Laura and S{\o}rensen, Helle and Tolver, Anders and Staicu, Ana-Maria},
  journal={Econometrics and Statistics},
  volume={23},
  pages={165--186},
  year={2022},
  doi={10.1016/j.ecosta.2021.07.003}
}
```

## License

MIT © Maria Laura Battagliola, Helle Sørensen, Anders Tolver, Ana-Maria Staicu
