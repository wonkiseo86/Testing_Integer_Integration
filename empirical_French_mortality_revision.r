## Loading packages
source("load_packages.r")

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

# plot

require(LaplacesDemon)
require(fdaACF)

#savefig("rainbow", width = 12, height = 10, toplines = 0.8, type = "png")
#plot(fts(0:105,t(logit(fseries[,,1]))), xlab = "Age", ylab = "Logit transformation of mortality rate")
#dev.off()

#savefig("acf", width = 12, height = 10, toplines = 0.8, type = "png")
#obtain_FACF(logit(replace(fseries[,,1], which(fseries[,,1]==0), 10^-6)), v = 0:105, nlags = 20)
#dev.off()

## Data preprocessing
X_series=fseries
findex=NULL
for (jj in 1:n_sub)
{
x_series = X_series[,,jj]
check=rowSums(x_series)
if (sum(is.na(check))>0){findex=append(findex,jj)}
}
X_series=mseries
mindex=NULL
for (jj in 1:n_sub)
{
x_series = X_series[,,jj]
check=rowSums(x_series)
if (sum(is.na(check))>0){mindex=append(mindex,jj)}
}


######################################################################################
### V0 and V1 tests for the subregions except 9th regtion (female & male)####
### Results for 9th region can be obtained using the supplementary code  ####
######################################################################################

#### Transformations#####################################################################
transformation=1  ## 1: logit, 2: probit, 3: no transformation, 4: log-transformation


uband=1/5
ql=0.04467116 
qu=2.12588475 

###########################################
X_series=fseries[,,setdiff(1:n_sub,findex)]

result=NULL 
for (jj in 1:length(setdiff(1:n_sub,findex)))
{

x_series = X_series[,,jj]
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
}






X_series=mseries[,,setdiff(1:n_sub,mindex)]
result2=NULL 
for (jj in 1:length(setdiff(1:n_sub,mindex)))
{

x_series = X_series[,,jj]
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
}

result
result2

## Reported numbers 0.5,1.0 and 1.5 represent d in (0,1), d=1, and d in (1,2), respectively.  
