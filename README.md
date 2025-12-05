# VARtime: A Multivariate Framework for Tail Risk Forecasting

`VARtime` implements a multivariate dynamic framework for forecasting extreme financial risk.  
The model integrates information from both extreme and non-extreme events and allows the inclusion of high-frequency realized measures—such as realized volatility, bipower variation, and semivariances—to better capture tail-risk dynamics.

This repository contains the R code required to estimate VARtime models and to replicate the empirical results of the accompanying research article:

**_VARtime: A Multivariate Framework for Tail Risk Forecasting_ (2025)**

---

## Features

- Multivariate tail-risk modeling using Weibull, Burr, or Pareto tail specifications  
- Vector autoregressive tail-index dynamics  
- Optional inclusion of realized measures:  
  - Realized volatility (RV)  
  - Bipower variation (BPV) and jumps  
  - Positive/negative semivariances  
- One-step-ahead forecasting of Value at Risk (VaR) and Expected Shortfall (ES)  
- Replication scripts for stock-market empirical applications  

---

## Installation

```r
devtools::install_github("yourusername/VARtime")
```

---

## Basic Example

```r
library(VARtime)

# Example dataset (placeholder)
data("tech_returns")

# Fit a VARtime model
fit <- vartime_fit(
  data     = tech_returns,
  tail     = "burr",
  realized = NULL   # options: "rv", "bj", "sv", "sv2"
)

summary(fit)
```

---

## Tail Risk Forecasting

```r
# Compute VaR and ES forecasts
forecast <- vartime_forecast(fit, horizon = 1)

forecast$VaR
forecast$ES
```

---

## Using Realized Measures

```r
fit_rv <- vartime_fit(
  data     = tech_returns,
  tail     = "weibull",
  realized = "rv"
)
```

---

## Reproducibility

This repository includes:

- Model estimation functions  
- Forecasting utilities  
- Scripts to replicate empirical results from the VARtime article  

More documentation will be added as the package evolves.

---

## Citation

If you use this code, please cite:

**Herrera et al. (2025). _VARtime: A Multivariate Framework for Tail Risk Forecasting_.**

---

## License

MIT License © Your Name
