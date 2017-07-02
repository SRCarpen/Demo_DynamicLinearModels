# Try time-varying AR via online DLM on simulated cores
# SRC 2017-05-14
# (c) Stephen R. Carpenter 2017

rm(list = ls())
graphics.off()

source('ODLMAR_2017-07-01.R')

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

# Detrend and standardize; may not be useful for the online method
#lm.revcore = lm(core ~ T.core)
#print('linear regression for detrending and standardizing:',quote=F)
#print(summary(lm.revcore))
# To prevent detrending, comment out the line below:
#core = lm.revcore$residuals / sd(lm.revcore$residuals)

title = c('Sim. Core, CV=45 pct, Online DLM')
nobs = length(core)

# START PROTOTYPE SHELL
# USER MUST INPUT: nl; delta; core; T.core; title

nl = 1 # number of lags
delta = 0.9 # 0<delta<1; see advice in functions

ODL.out = ODLMAR(nl,delta,core,T.core,title)

Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]
B.sd = ODL.out[[4]]
