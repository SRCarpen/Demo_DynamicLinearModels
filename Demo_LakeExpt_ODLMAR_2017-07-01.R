# Try time-varying AR on Cascade data via online DLM method
# SRC 2017-07-01
# (c) Stephen R. Carpenter 2017-07-01

rm(list = ls())
graphics.off()

source('ODLMAR_2017-07-01.R')

# Load Peter Lake demo data set
# save(x.full,T.full,file='PeterL_2015_manual_chl.Rdata')
load(file='PeterL_2015_BGA_HYLB.Rdata')

# log transform? Optional. Comment out to use natural unit
x.full = log10(x.full)
title = c('Peter 2015 log10 Phycocyanin (HYLB), Online DLM AR(p)')
nobs = length(x.full)

# START PROTOTYPE SHELL
# USER MUST INPUT: nl; delta; x.full; T.full; title

nl = 1 # number of lags
delta = 0.9 # 0<delta<1; see advice in functions

ODL.out = ODLMAR(nl,delta,x.full,T.full,title)

Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]
B.sd = ODL.out[[4]]
