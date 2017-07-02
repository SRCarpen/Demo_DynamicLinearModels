# Try time-varying AR via ML DLM with MARSS on simulated cores
# SRC 2017-05-14
# (c) Stephen R. Carpenter 2017

rm(list = ls())
graphics.off()

library('MARSS')

source('MLE-AR-DLM_2017-07-01.R')

# Load data simulated by Eutrophication_PNAS-1D_ShortCoreTest_V0_2016-08-27.R
# Tvec and Wsim are time and true biomass; 
# T.core and core are the mixed & compressed core sample
# 400 years, switch occurs year 300
# save(Tvec,Wsim,T.core,core,file='Core_simulated.Rdata')
load(file='Core_simulated_CV45pct_2017-05-12.Rdata') # file with CV = 0.45 & flickering
L.core = length(core)

# Set up for DLMs
# remember that the core vector starts at the surface sample and runs to the deepest sample

# REVERSING TIME
core=rev(core)
T.core = rev(T.core)

# Detrend and standardize; may be helpful for the MLE method
lm.revcore = lm(core ~ T.core)
print('linear regression for detrending and standardizing:',quote=F)
print(summary(lm.revcore))
# To prevent detrending, comment out the line below
core = lm.revcore$residuals / sd(lm.revcore$residuals)

title = c('Sim. Core, CV=45 pct, MLE DLM')
nobs = length(core)

# set up for function call
nl = 3 # number of lags = AR order

# Compute dlm
dlmfit = DLM.MLE(nl,nobs,T.core,core,title)
Yyhat = dlmfit[[1]]
EigenVals = dlmfit[[2]]
B.est = dlmfit[[3]]
B.se = dlmfit[[4]]
dlm.struc = dlmfit[[5]]

