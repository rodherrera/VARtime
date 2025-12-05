
1. The VARtime model and the backtesting procedures are implemented using four files located in the Main estimation folder.

2. This folder contains four files: one main script, VARtime, and three auxiliary scripts: filter.arma, Likelihood, and Accuracy Tests.

3. The filter.arma file includes the code for a first-order autoregressive specification that describes how tail losses evolve over time in a multivariate setting.

4. The Likelihood file contains the code for the joint likelihood function, which is used to estimate the model parameters and their standard errors.

5. The Accuracy Tests file implements the backtesting procedures for both Value-at-Risk (VaR) and Expected Shortfall (ES), as well as the corresponding scoring functions.

6. The VARtime file loads the return data and realized measures, splits the sample into in-sample and out-of-sample periods, and then estimates the models. It also runs the accuracy tests and scoring-function estimation. Finally, it constructs the Model Confidence Set and generates the key figure presented in the paper. (This figure can also be produced using the files in the Plots folder.)
In addition, VARtime provides estimates for the extension that incorporates bipower variation and jump components using a Pareto regularly varying function.