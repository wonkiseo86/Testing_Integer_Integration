## Description: Simulation codes for eigenvalue ratio estimators and variance ratio tests.
## The code can produce the results given in Tables 1 and 2 of the paper.

## Loading packages

source("load_packages.r")
setwd("Path")

inner = dget("auxiliary/inprod.R")
lrvar = dget("auxiliary/lr_var_v2_for_fractional.R")



nonstat=1
addmargin=0.01
 bdd=0.15;  ###########################################################################################
 bdd2=0.15;
 
sim_DGP <- function(seed_number, sample_size, d, grid_number)
 {
 bunin=0
 T=sample_size+bunin; nt=grid_number; t = (0:(nt-1))/(nt-1); d2=d*(d>-1/2 & d<1/2) + (0.25)*(d>=1/2);  

 
 YY=NULL
 for (da in c(d2,runif(1,max(d2-0.2,-0.5),max(d2-0.1,-0.5+addmargin))))
 #for (da in c(d2,max(d2-0.2,-0.4)))
 {
 if (da==d2){randim=trunc(runif(1,1,4))}else{randim=trunc(runif(1,1,4))}
 for (addind in 1:randim)
 {
 #YY=rbind(YY,arfima.sim(n=T, model = list(phi = runif(1,-bdd,bdd), theta = runif(1,-bdd,bdd) , dfrac = da,  dint = 0, n.burn=200, mu=rnorm(1,0,1))))
 aatem=fracdiff.sim(T, ar = runif(1,-bdd,bdd), ma = runif(1,-bdd,bdd) , d = da, n.start =100)$series
 YY=rbind(YY,aatem + rnorm(1))
 }
 }
 for (addind in (nrow(YY)+1):26)
 {
 
 #YY=rbind(YY,arfima.sim(n=T, model = list(phi = runif(1,-bdd,bdd), theta = runif(1,-bdd,bdd) , dfrac = max(d2-0.15,-0.4),  dint = 0, n.burn=200, mu=rnorm(1,0,1))))
 aatem=fracdiff.sim(T, ar = runif(1,-bdd,bdd), ma = runif(1,-bdd,bdd) ,  d = runif(1,max(d2-0.5,-0.5),max(d2-0.2,-0.5+addmargin)), n.start =100)$series
 YY=rbind(YY,aatem+rnorm(1))
 }
 YY=YY[,(bunin+1):ncol(YY)]
 
 if(d>=1/2){ d3=d;
 for (dd in c(runif(1,max(d3-0.2,0.5),max(d3-0.1,0.5+addmargin)),d3))
 #for (dd in c(max(d3-0.2,0.5),d3))
 {
 if (dd==d3){randim=trunc(runif(1,1,4))}else{randim=trunc(runif(1,1,4))}
 #randim=trunc(runif(1,1,3))
 for (addind in 1:randim)
 {
 aatem=fracdiff.sim(T, ar = runif(1,-bdd,bdd), ma = runif(1,-bdd,bdd) ,  d = dd-1, n.start =100)$series
 YY=rbind(cumsum(aatem)+rnorm(1),YY)
 }
 YY=YY[1:26,(bunin+1):ncol(YY)]
 }
 }
 #varseq=c(runif(cutsd,1-margin,1+margin),runif(length((1:(nrow(YY)-cutsd))),1-margin,1+margin)*decrea^((1:(nrow(YY)-cutsd))))
 varseq=1/((1:26)^(2))
 #YY=t(t(YY) %*%diag((0.9)^{0:(nrow(YY)-1)}))
 YY=t(t(YY) %*%diag(varseq))
 set1=sample(1:5,5)
 set2=setdiff(1:26,set1)
 rpert=append(set1,set2)
 LBFP=LBF[,rpert]
 data1=LBFP%*%YY
 
  #return(data1)
 }
 



decrea=0.5
margin=0
cutsd=1
uband=1/5 #1/3 and 1/4


#ql=0.04467116 
#qu=2.12588475 

ql=0.03038333
qu=0.57504599
ql2=qnorm(0.025)
qu2=qnorm(0.975)
bw1=0.65




#d_sim = 0.6,0.8,1,1.2,1.4
#d_sim = -0.4,-0.2,0,0.2,0.4


Dresults=NULL

Dset1=c(-0.45,-0.3,-0.15,0,0.15,0.3,0.45)
Dset2=c(1-0.45,1-0.3,1-0.15,1,1.15,1.3,1.45)
 
if(nonstat==1){DDset=Dset2}else{DDset=Dset1} 

set.seed(99999999)
for(d_sim in DDset)  ##########################################
{
REPORT=NULL; REPORT2=REPORT ; REPORT3=REPORT


maxiter=2000

T_sample=c(125,250,500,750,1000)

T_max=max(T_sample)
TEST_RESULT=matrix(ncol=length(T_sample),nrow=maxiter); TEST_RESULT2=TEST_RESULT ; TEST_RESULT3=TEST_RESULT
TESTSTAT=matrix(ncol=length(T_sample),nrow=maxiter); TESTSTAT2=TESTSTAT; TESTSTAT3=TESTSTAT
lbnumber2=26; nt = 150 ; t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}

seed_number=1
for (seed_number in 1:maxiter)
{

x_mat=sim_DGP(1000*seed_number, T_max, d_sim ,nt)
x_mat
hh=t(LBF[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
#hh2=t(LBF[2:(nt),])%*%y_mat[2:(nt),]*(t[2]-t[1])
xcoef=hh
xcoef=t(xcoef)


for (TTT in (T_sample))
{
bandw=floor(TTT^(uband))

zz0=t(xcoef[1:TTT,])
#zz0_undemeaned=t(zz0)
zz0=zz0 -rowMeans(zz0)
zz0=t(zz0)
zz=zz0

zz2=apply(zz0,2, cumsum)
vx_eigen = crossprod(zz2) / (TTT^2)
eelements=eigen(vx_eigen)
ev0=eelements$vectors[,1]  


if(d_sim>1/2){zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] }
zz2=apply(zz0,2, cumsum)
v0=lrvar(zz,kernel=2)$omega / TTT
vx1 = crossprod(zz2) / (TTT^2)
vx0 = lrvar(zz,kernel=2)$omega / TTT

eval_vx0= t(ev0)%*%vx0%*%ev0
eval_vx1= t(ev0)%*%vx1%*%ev0

teststat=(eval_vx1/eval_vx0)


xxcoef=xcoef[1:TTT,]-rowMeans(xcoef[1:TTT,])


if(d_sim<1/2 & d_sim>-1/2){
if(d_sim>0){corrind=1}
if(d_sim<0){corrind=-1}
if(d_sim==0){corrind=0} }

if(d_sim>=1/2){
if(d_sim>1){corrind=1}
if(d_sim<1){corrind=-1}
if(d_sim==1){corrind=0} }


if (d_sim>1/2){ql=0.04467116;qu=2.12588475 }else{ql=0.03038333;qu=0.57504599} 



TESTSTAT[seed_number,which(TTT==T_sample)]=teststat

if (teststat > ql & teststat<qu){TEST_RESULT[seed_number,which(TTT==T_sample)]=0}
if (teststat < ql){TEST_RESULT[seed_number,which((TTT)==T_sample)]=-1}
if (teststat > qu){TEST_RESULT[seed_number,which((TTT)==T_sample)]=1}


}

if(seed_number%%200==0){print(c(d_sim,sum(TEST_RESULT[1:seed_number,1]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,2]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,3]==corrind)/seed_number, TESTSTAT[seed_number,]))}

}

REPORT=rbind(REPORT,c(sum(TEST_RESULT[1:seed_number,1]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,2]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,3]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,4]==corrind)/seed_number,sum(TEST_RESULT[1:seed_number,4]==corrind)/seed_number))
Dresults=rbind(Dresults,REPORT)
}

AA=t(Dresults)

AA[,4]=1-AA[,4]
round(AA,digits=3)
