# VARtime: A Multivariate Framework for Tail Risk Forecasting 📉📈

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.0.0](https://img.shields.io/badge/R-%3E%3D%204.0.0-276DC3.svg)](https://www.r-project.org/)
[![GitHub issues](https://img.shields.io/github/issues/rodherrera/VARtime.svg)](https://github.com/rodherrera/VARtime/issues)
[![GitHub last commit](https://img.shields.io/github/last-commit/rodherrera/VARtime.svg)](https://github.com/rodherrera/VARtime)
[![Made with R](https://img.shields.io/badge/Made%20with-R-276DC3.svg)](https://www.r-project.org/)

`VARtime` implements a multivariate dynamic framework for forecasting extreme financial risk.  
The model integrates information from both extreme and non-extreme events and incorporates high-frequency realized measures—such as realized volatility, bipower variation, and semivariances—to capture shifts in tail-risk dynamics with greater accuracy.

This repository contains the full R implementation used to estimate VARtime models, evaluate forecasting performance, and replicate the empirical results of the research article:

**Candia, C., Herrera, R., & Clements, A. (2025). _VARtime: A Multivariate Framework for Tail Risk Forecasting_. Submitted.**

---

## 🚀 Features

- Multivariate tail-risk modeling using Weibull, Burr, or Pareto tail specifications  
- Vector autoregressive tail-index dynamics  
- Optional inclusion of realized measures:  
  - Realized volatility (RV)  
  - Bipower variation (BPV) and jump components  
  - Positive and negative semivariances  
- One-step-ahead forecasts of Value at Risk (VaR) and Expected Shortfall (ES)  
- Full replication of empirical results for major U.S. technology stocks  
- Construction of Model Confidence Sets (MCS) and scoring-function evaluations  

---

## 🗂️ Repository Structure

### **Main Estimation Folder** 📦  
The core implementation of the VARtime framework is organized into four files:

---

### **1. VARtime** (Main Script)
- Loads return data and realized measures  
- Splits the sample into in-sample and out-of-sample periods  
- Estimates all VARtime specifications  
- Runs VaR and ES accuracy tests  
- Computes scoring-function evaluations  
- Constructs the **Model Confidence Set (MCS)**  
- Generates the main figure presented in the paper  
- Includes the extended specification using bipower variation and jump components under a Pareto tail  

---

### **2. filter.arma**
- Implements a first-order autoregressive mechanism  
- Governs the evolution of multivariate tail losses  
- Captures persistence and cross-asset spillovers in extremes  

---

### **3. Likelihood**
- Contains the joint likelihood function used for parameter estimation  
- Supports all tail specifications (Weibull, Burr, Pareto)  
- Computes standard errors for all estimated parameters  

---

### **4. Accuracy Tests**
- Implements backtesting procedures for **VaR** and **ES**  
- Includes strictly consistent scoring functions  
- Implements regression-based ES tests (ERS, sCC, iERS, etc.)

---

### **Plots Folder** 🖼️  
Contains scripts that reproduce the figures used in the paper, including the main comparative forecasting figure also produced in the VARtime script.

---

## 🔁 Reproducibility

This repository includes:

- Model estimation functions  
- Forecasting utilities  
- Backtesting and scoring-function procedures  
- Scripts to fully replicate the empirical analysis  

Additional documentation will be added as the package evolves.

---

## 📚 Citation

**Candia, C., Herrera, R., & Clements, A. (2025). _VARtime: A Multivariate Framework for Tail Risk Forecasting_. Submitted.**

---

## 📝 License

MIT License © Rodrigo Herrera
