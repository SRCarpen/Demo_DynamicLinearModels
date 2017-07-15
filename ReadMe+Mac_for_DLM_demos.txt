Demonstrations of Dynamic Linear Models
Stephen R. Carpenter July 2017

Two examples (one a simulated paleo data set, the other whole-lake experimental data)
are analyzed with autoregressions fit as Dynamic Linear Models using two methods: 
Online fitting (Pole et al. 1994) and Maximum Likelihood using a state space model
(Holmes et al. 2014).

The rationale for using time-varying regressions (such as DLMs) to estimate eigenvalues
of time series is based on Ives and Dakos (2012).

Whole-lake data come from Pace et al. 2017.

Data sets:

Core_simulated_CV45pct_2017-07-01.Rdata is a simulated paleo record as a lake becomes
eutrophic:

PeterL_2015_BGA_HYLB.Rdata is daily average phycocyanin data from Peter Lake during
the experiment described by Pace et al. 2017.

Functions for fitting AR(p) models as DLMs:

MLE-AR-DLM_2017-07-01.R is a maximum-likelihood state-space method using MARSS
(Holmes et al. 2014). TO RUN ON McINTOSH: Replace "windows" with "quartz" globally.

ODLMAR_2017-07-01.R is the online method of Pole et al. (1994).
TO RUN ON McINTOSH: Replace "windows" with "quartz" globally.

Demonstration programs:

Demo_LakeExpt_MLE-AR-DLM-function_2017-07-01.R analyzes the Peter Lake time series
using the maximum likelihood state-space method.

Demo_LakeExpt_ODLMAR_2017-07-01.R analyzes the Peter Lake time series using
the online method.

Demo_Paleo_MLE-AR-DLM-function_2017-07-01.R analyzes the simulated paleo core
using the maximum likelihood state space method.

Demo_Paleo_ODLMAR_2017-07-01.R analyzes the simulated paleo core using the
online method.

References:

Holmes EE, Ward EJ, Scheuerell MD, 2014, Analysis of Multivariate Time Series Using 
the MARSS package. V 3.9. Northwest Fisheries Science Center, NOAA, Seattle, WA, USA

Ives, A. R., and V. Dakos. 2012. Detecting dynamical changes in nonlinear time series 
using locally linear state-space models. Ecosphere 3:art58.

Pace, M. L., R. D. Batt, C. D. Buelo, S. R. Carpenter, J. J. Cole, J. T. Kurtzweil, 
and G. M. Wilkinson. 2017. Reversal of a cyanobacterial bloom in response to 
early warnings. Proceedings of the National Academy of Sciences 114:352-357.

Pole, A., M. West and J. Harrison. 1994. Applied Bayesian Forecasting and Time 
Series Analysis. Chapman and Hall, N.Y. 409 pp.
