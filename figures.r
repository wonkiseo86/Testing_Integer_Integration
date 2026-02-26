## Description: codes for Figures 1 and 2.

#########################################
#### Figure 1: Canadian Yield Curves ####
#########################################
source("load_packages.r") 
# - Note: This script only loads standard CRAN packages (e.g., fda, etc.) required for the analysis. No custom functions are defined herein.

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






#########################################################
## Figure 2: French Mortality Rates and functional ACF ##
#########################################################
source("load_packages.r") 
# - Note: This script only loads standard CRAN packages (e.g., fda, etc.) required for the analysis. No custom functions are defined herein.

n_sub=22
fseries = array(0,dim = c(121,106,n_sub))
mseries = array(0,dim = c(121,106,n_sub))

index=3 ## read mx	
female=read.table( "data/NUT2/1f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/1m.txt",header=TRUE,sep=";",dec=","); hhj=1;    miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/2f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/2m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/3f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/3m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/4f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/4m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/5f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/5m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/6f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/6m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/7f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/7m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/8f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/8m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "data/NUT2/9f.txt",header=TRUE,sep=";",dec=","); male=read.table( "data/NUT2/9m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/10f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/10m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/11f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/11m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/12f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/12m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/13f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/13m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/14f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/14m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/15f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/15m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/16f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/16m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/17f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/17m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/18f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/18m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/19f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/19m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/20f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/20m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/21f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/21m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("data/NUT2/22f.txt",header=TRUE,sep=";",dec=","); male=read.table("data/NUT2/22m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2

inner = dget("auxiliary/inprod.R")
lrvar = dget("auxiliary/lr_var_v2_for_fractional.R")

# plot
require(LaplacesDemon)
require(fdaACF)

plot(fts(0:105,t(logit(fseries[,,1]))), xlab = "Age", ylab = "Logit transformation of mortality rate")  ## Figure 2 Left
obtain_FACF(logit(replace(fseries[,,1], which(fseries[,,1]==0), 10^-6)), v = 0:105, nlags = 20) ## Figure 2 Right
