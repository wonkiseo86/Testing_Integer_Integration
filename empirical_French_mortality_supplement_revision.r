## Loading packages
source("load_packages.r") 
# - Note: This script only loads standard CRAN packages (e.g., fda, etc.) required for the analysis. No custom functions are defined herein.


# - Note: Load data 
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
# - usage: inner(a,b,c)
# - Description: Approximates the L2 inner product \int a(u)b(u) du using the Trapezoidal Rule for numerical integration.
# - Inputs:
#    a, b: Function values evaluated on grid 'c'.
#    c: Regularly spaced grid points (Numeric vector).
# - Output: Numeric scalar (Approximated integral).
# - Assumptions: 'c' is a constant-step grid; 'a' and 'b' have the same length.

lrvar = dget("auxiliary/lr_var_v2_for_fractional.R")
# - usage: lr_var(u, kernel)
# - Desc: Computes the long-run covariance matrix (Omega) required for test statistics.
# - Inputs: 
#    u: Input data matrix
#    kernel: Kernel type index 
# - Output: A list containing 'omega' (the estimated long-run covariance matrix).
# - Assumptions: Bartlett kernel is used. 

uband=1/5
ql=0.04467116 
qu=2.12588475 

# Transformations#####################################################################
transformation=1  ## 1: logit, 2: probit, 3: no transformation, 4: log-transformation


######################################################################################
### Section 1: V0 and V1 tests for the NINTH subregion (female & male)####
######################################################################################
# - Note: Setup for female data
X_series=fseries[21:121,,9]
findex=NULL
mindex=NULL

transformation=4 ## 1 for the results for Alsace in Table 5 of the main manuscript ; 3 (resp. 4) for the results for Alsace in Table 1 (resp. Table 2) of the Supplementary Material; 

result=NULL 

# - Note: Data preprocessing
x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)



# - Note: Transformation of the female data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note:Implementatino of V0 test (female data)
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)

if (teststat < ql) {result=append(result,-1)}
if (teststat > ql & teststat <qu) {result=append(result,0)}

if (teststat > qu){
# - Note: V1 test for female data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result=append(result,1)}

if (teststat2 > qu) {result=append(result,1.5)}
if (teststat2 < ql) {result=append(result,0.5)}





# - Note: Setup for male data
X_series=mseries[21:121,,9]
result2=NULL 

# - Note: Data preprocessing
x_series = X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)


# - Note: Transformation of the female data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note: Implementation of the V0 test (male data)
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)


if (teststat < ql) {result2=append(result2,-1)}
if (teststat > ql & teststat <qu) {result2=append(result2,0)}

if (teststat > qu){
# - Note: V1 test for male data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result2=append(result2,1)}

if (teststat2 > qu) {result2=append(result2,1.5)}
if (teststat2 < ql) {result2=append(result2,0.5)}


# - Note: Report results collectively.
result
result2







######################################################################################
######################################################################################
### Section 2: V0,V1 and V2 tests together for the FIRST subregion (female)####
######################################################################################
######################################################################################
## This part is to just confirm that d = 2 or less for the data set of log mortality rates

# - Note: Basic setup and data preprocessing
X_series=fseries[1:121,,1]
findex=NULL
mindex=NULL

transformation = 1
result=NULL 
x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)



# - Note: Transformation of the female data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note: Implementation of the V0 test (female)
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)

if (teststat < ql) {result=append(result,-1)}
if (teststat > ql & teststat <qu) {result=append(result,0)}


if (teststat > qu){
# - Note: V1 test for female data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result=append(result,1)}
if (teststat2 > qu) {result=append(result,1.5)}
if (teststat2 < ql) {result=append(result,0.5)}


# - Note: V2 test for female data when V1 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat3=(eval_vx1/eval_vx0)
teststat3 
if (teststat3 > qu) {print("Reject at upper tail")}
if (teststat3 < ql) {print("Reject at lower tail")}
if (teststat3 > ql & teststat3 <qu) {print("Accept")}





######################################################################################
######################################################################################
### Section 3: V0,V1 and V2 tests together for the SEVENTH subregion (female)####
######################################################################################
######################################################################################
## This part is to just confirm that d = 2 or less for the data set of log mortality rates

# - Note: Basic setup and data preprocessing
X_series=fseries[1:121,,7]
findex=NULL
mindex=NULL

transformation = 1

result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)



# - Note: Transformation of the female data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note: Implementation of V0 test for female data 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)

if (teststat < ql) {result=append(result,-1)}
if (teststat > ql & teststat <qu) {result=append(result,0)}


if (teststat > qu){
# - Note: V1 test for female data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result=append(result,1)}
if (teststat2 > qu) {result=append(result,1.5)}
if (teststat2 < ql) {result=append(result,0.5)}


# - Note: V2 test for female data when V1 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat3=(eval_vx1/eval_vx0)
teststat3 
if (teststat3 > qu) {print("Reject at upper tail")}
if (teststat3 < ql) {print("Reject at lower tail")}
if (teststat3 > ql & teststat3 <qu) {print("Accept")}






######################################################################################
######################################################################################
### Section 4: V0,V1 and V2 tests together for the 21TH subregion (male)####
######################################################################################
######################################################################################
## This part is to just confirm that d = 2 or less for the data set of log mortality rates


# - Note: Basic parameter setup and data preprocessing
X_series=mseries[1:121,,21]
findex=NULL
mindex=NULL

transformation = 1
result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)


# - Note: Transformation of the male data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note: Implementation of the V0 test (male data)
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)

if (teststat < ql) {result=append(result,-1)}
if (teststat > ql & teststat <qu) {result=append(result,0)}


if (teststat > qu){
# - Note: V1 test for male data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result=append(result,1)}
if (teststat2 > qu) {result=append(result,1.5)}
if (teststat2 < ql) {result=append(result,0.5)}


# - Note: V2 test for female data when V1 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat3=(eval_vx1/eval_vx0)
teststat3 
if (teststat3 > qu) {print("Reject at upper tail")}
if (teststat3 < ql) {print("Reject at lower tail")}
if (teststat3 > ql & teststat3 <qu) {print("Accept")}






######################################################################################################
######################################################################################################
### Section 5: V0,V1 and V2 tests for log-transformation (revision, for the Supplementary Appendix) ##
######################################################################################################
######################################################################################################
## This part is to just confirm that d = 2 or less for the data set of log mortality rates
region_index=21   ## Change this to "1,5,7,9,13,21 for female" & "1,21 for male". 
gender=2 # 1 for female, 2 for male
transformation=4 

# - Note: Basic parameter setup and data preprocessing
if(gender==1 & region_index==9){X_series=fseries[21:121,,region_index]}
if(gender==1 & region_index!=9){X_series=fseries[1:121,,region_index]}
if(gender==2 & region_index==9){X_series=mseries[21:121,,region_index]}
if(gender==2 & region_index!=9){X_series=mseries[1:121,,region_index]}

findex=NULL
mindex=NULL

result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]

x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)


# - Note: Transformation of the female data following the "transformation" parameter
if (transformation==1){
x_mat=t(log(x_series/(1-x_series)))}

if (transformation==2){
x_mat=t(qnorm(x_series,0,1))}

if (transformation==3){
x_mat=t(x_series)}

if (transformation==4){
x_mat=t(log(x_series))}

T_sample=ncol(x_mat)
TTT=T_sample
nt=nrow(x_mat)


# - Note: Construct basis functions to be used
lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


# - Note: Implementatino of V0 test 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh; xcoef=xcoef[,2:ncol(xcoef)]-xcoef[,1]
xcoef=t(xcoef)

T_sample=nrow(xcoef); TTT=T_sample

bandw=floor(TTT^(uband))
zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eelements=eigen(vx1)
ev0=eelements$vectors[,1]  

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0
teststat=(eval_vx1/eval_vx0)

if (teststat < ql) {result=append(result,-1)}
if (teststat > ql & teststat <qu) {result=append(result,0)}

if (teststat > qu){
# - Note: V1 test for data when V0 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat2=(eval_vx1/eval_vx0)}

if (teststat2 > ql & teststat2 <qu) {result=append(result,1)}
if (teststat2 > qu) {result=append(result,1.5)}
if (teststat2 < ql) {result=append(result,0.5)}


# - Note: V2 test for data when V1 test is rejected at an upper tail
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 

zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat3=(eval_vx1/eval_vx0)
teststat3 
if (teststat3 > qu) {print("Reject at upper tail")}
if (teststat3 < ql) {print("Reject at lower tail")}
if (teststat3 > ql & teststat3 <qu) {print("Accept")}






