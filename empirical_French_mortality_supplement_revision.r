## Loading packages
library(fda);library(tseries);library(sandwich);library(sde);library(variables);library(basefun);library(polynom);library(fracdiff);library(LongMemoryTS);library(arfima)

## Relevant path needs to be set.
setwd("Path/data/NUT2")


n_sub=22
fseries = array(0,dim = c(121,106,n_sub))
mseries = array(0,dim = c(121,106,n_sub))

index=3 ## read mx	
female=read.table( "1f.txt",header=TRUE,sep=";",dec=","); male=read.table( "1m.txt",header=TRUE,sep=";",dec=","); hhj=1;    miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "2f.txt",header=TRUE,sep=";",dec=","); male=read.table( "2m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "3f.txt",header=TRUE,sep=";",dec=","); male=read.table( "3m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "4f.txt",header=TRUE,sep=";",dec=","); male=read.table( "4m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "5f.txt",header=TRUE,sep=";",dec=","); male=read.table( "5m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "6f.txt",header=TRUE,sep=";",dec=","); male=read.table( "6m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "7f.txt",header=TRUE,sep=";",dec=","); male=read.table( "7m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "8f.txt",header=TRUE,sep=";",dec=","); male=read.table( "8m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table( "9f.txt",header=TRUE,sep=";",dec=","); male=read.table( "9m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("10f.txt",header=TRUE,sep=";",dec=","); male=read.table("10m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("11f.txt",header=TRUE,sep=";",dec=","); male=read.table("11m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("12f.txt",header=TRUE,sep=";",dec=","); male=read.table("12m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("13f.txt",header=TRUE,sep=";",dec=","); male=read.table("13m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("14f.txt",header=TRUE,sep=";",dec=","); male=read.table("14m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("15f.txt",header=TRUE,sep=";",dec=","); male=read.table("15m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("16f.txt",header=TRUE,sep=";",dec=","); male=read.table("16m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("17f.txt",header=TRUE,sep=";",dec=","); male=read.table("17m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("18f.txt",header=TRUE,sep=";",dec=","); male=read.table("18m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("19f.txt",header=TRUE,sep=";",dec=","); male=read.table("19m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("20f.txt",header=TRUE,sep=";",dec=","); male=read.table("20m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("21f.txt",header=TRUE,sep=";",dec=","); male=read.table("21m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2
female=read.table("22f.txt",header=TRUE,sep=";",dec=","); male=read.table("22m.txt",header=TRUE,sep=";",dec=","); hhj=hhj+1;miny=min(female[,1]);maxy=max(female[,1]); data1=NULL;data2=NULL;yy=miny;for (jjj in miny:maxy){data1=rbind(data1,female[(female[,1]==jjj),index]);data2=rbind(data2,male[(male[,1]==jjj),index])};fseries[,,hhj]=data1; mseries[,,hhj]=data2

setwd("C:/Users/wseo2199/Dropbox/Shang_Seo/Integer_integration/code/")
inner = dget("inprod.R")
lrvar = dget("lr_var_v2_for_fractional.R")

uband=1/5
ql=0.04467116 
qu=2.12588475 

# Transformations#####################################################################
transformation=4  ## 1: logit, 2: probit, 3: no transformation, 4: log-transformation




######################################################################################
### V0 and V1 tests for the NINTH subregion (female & male)####
######################################################################################
X_series=fseries[21:121,,9]
findex=NULL
mindex=NULL


result=NULL 



x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]


x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)




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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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






X_series=mseries[21:121,,9]
result2=NULL 

x_series = X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]


x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)

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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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

result
result2





######################################################################################
######################################################################################
### V0,V1 and V2 tests together for the FIRST subregion (female)####
######################################################################################
######################################################################################
X_series=fseries[1:121,,1]
findex=NULL
mindex=NULL


result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]


x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)

#x_mat=t(x_series)
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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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

### V_2 test## 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

## second difference ##
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
######################################################

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
### V0,V1 and V2 tests together for the SEVENTH subregion (female)####
######################################################################################
######################################################################################
X_series=fseries[1:121,,7]
findex=NULL
mindex=NULL


result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]


x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)

#x_mat=t(x_series)
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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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

### V_2 test## 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

## second difference ##
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
######################################################

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
### V0,V1 and V2 tests together for the 21TH subregion (male)####
######################################################################################
######################################################################################
X_series=mseries[1:121,,21]
findex=NULL
mindex=NULL


result=NULL 

x_series =X_series
subindex=is.na(rowSums(x_series))
if (sum(subindex)>=1){rendpoint=max(which(subindex==1))} else{rendpoint=0}
x_series=x_series[(rendpoint+1):nrow(x_series),]


x_series = replace(x_series, which(x_series == 0), 10^-4)
x_series = replace(x_series, which(x_series >= 1), 1-10^-4)

#x_mat=t(x_series)
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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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

### V_2 test## 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

## second difference ##
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
######################################################

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
### V0,V1 and V2 tests together for log-transformation (revision, for the Supplementary Appendix) ####
######################################################################################################
######################################################################################################
region_index=9
gender=1 # 1 for female, 2 for male
transformation=4 


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

#x_mat=t(x_series)
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


lbnumber2=26;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}


### V_0 test## 
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
### V_1 test## 
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

### V_2 test## 
hh=t(LBF[1:nt,])%*%x_mat[1:nt,]*(t[2]-t[1])
xcoef=hh;
xcoef=t(xcoef)
T_sample=nrow(xcoef); TTT=T_sample
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
zz0=zz0 #-rowMeans(zz0)
zz0=t(zz0)
zz=zz0

## second difference ##
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
TTT=TTT-1
zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
######################################################

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

