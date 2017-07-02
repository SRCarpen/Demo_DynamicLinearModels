# Try time-varying AR on Cascade data using ML DLM method with MARSS
# SRC 2017-07-01
# (c) Stephen R. Carpenter 2017-07-01

rm(list = ls())
graphics.off()

library('MARSS')

source('MLE-AR-DLM_2017-07-01.R')

# Load Peter Lake demo data set
# save(x.full,T.full,file='PeterL_2015_manual_chl.Rdata')
load(file='PeterL_2015_BGA_HYLB.Rdata')

# log transform? Optional. Comment out to use natural unit
x.full = log10(x.full)
title = c('Peter 2015 log10 Phycocyanin (HYLB), ML DLM AR(p)')
nobs = length(x.full)

# Detrend and standardize (optional)
linmod = lm(x.full ~ T.full)
err.std = linmod$residuals/sd(linmod$residuals)
x.full = err.std  # comment out this line to NOT detrend

# set up for function call
nl = 2 # number of lags = AR order
X = x.full
timevec = T.full

# Compute dlm
dlmfit = DLM.MLE(nl,nobs,timevec,X,title)
Yyhat = dlmfit[[1]]
EigenVals = dlmfit[[2]]
B.est = dlmfit[[3]]
B.se = dlmfit[[4]]
dlm.struc = dlmfit[[5]]

