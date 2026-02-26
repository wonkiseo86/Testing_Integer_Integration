## Description: codes for Figures 1 and 2.

## Figure 1: Canadian Yield Curves ##

source("load_packages.r")
load("CV1.Rdata") 
data_raw   = read.csv("data/Canadian_daily_yields.csv")
data_raw=as.matrix(data_raw[,2:ncol(data_raw)])

inner = dget("auxiliary/inprod.R")
lrvar = dget("auxiliary/lr_var_v2_for_fractional.R")

uband=1/5

ql=0.04467116 
qu=2.12588475 

T_sample=nrow(data_raw)
TTT=T_sample

x_mat=t(data_raw)
nt=nrow(x_mat)

# Figure 1
plot(fts(grid_point, x_mat), xlab = "Maturity", ylab = "Zero-coupon bond")




## Figure 2: French Mortality Rates and functional ACF ##
