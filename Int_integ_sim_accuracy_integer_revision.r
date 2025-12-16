## Description: Simulation codes for eigenvalue ratio estimators and variance ratio tests.
## The code can produce the results given in Tables 1 and 2 of the paper.

## Loading packages

source("load_packages.r")
setwd("Path")

inner = dget("auxiliary/inprod.R")
lrvar = dget("auxiliary/lr_var_v2_for_fractional.R")


 bdd=0.15;  ###########################################################################################
 bdd2=0.15;
addmargin=0.01
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
 YY=rbind(YY,aatem)
 }
 }
 for (addind in (nrow(YY)+1):26)
 {
 
 #YY=rbind(YY,arfima.sim(n=T, model = list(phi = runif(1,-bdd,bdd), theta = runif(1,-bdd,bdd) , dfrac = max(d2-0.15,-0.4),  dint = 0, n.burn=200, mu=rnorm(1,0,1))))
 aatem=fracdiff.sim(T, ar = runif(1,-bdd,bdd), ma = runif(1,-bdd,bdd) ,  d = runif(1,max(d2-0.5,-0.5),max(d2-0.2,-0.5+addmargin)), n.start =100)$series
 YY=rbind(YY,aatem)
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
 YY=rbind(cumsum(aatem),YY)
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
 
     return(data1)
 }
 

decrea=0.5
margin=0
cutsd=1
uband=1/5 #1/3 and 1/4


ql=0.04467116 
qu=2.12588475 

#ql=0.03038333
#qu=0.57504599
ql2=qnorm(0.025)
qu2=qnorm(0.975)
bw1=0.65




#d_sim = 0.6,0.8,1,1.2,1.4
#d_sim = -0.4,-0.2,0,0.2,0.4

set.seed(99999)
Dresults=NULL
maxiter=2000
d_rand=sample(c(0,1),maxiter,replace=TRUE)
T_sample=c(125,250,500,750,1000)
REPORT=NULL; REPORT2=REPORT ; REPORT3=REPORT


TEST_RESULT=matrix(ncol=length(T_sample),nrow=maxiter); TEST_RESULT2=TEST_RESULT ; TEST_RESULT3=TEST_RESULT
TESTSTAT=matrix(ncol=length(T_sample),nrow=maxiter); TESTSTAT2=TESTSTAT; TESTSTAT3=TESTSTAT
for(iijj in 1:maxiter)
{
d_sim=d_rand[iijj]
T_max=max(T_sample)
lbnumber2=26; nt = 150 ; t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}

seed_number=iijj
x_mat=sim_DGP(1000*seed_number, T_max, d_sim ,nt)
hh=t(LBF[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
#hh2=t(LBF[2:(nt),])%*%y_mat[2:(nt),]*(t[2]-t[1])
xcoef=hh
ycoef=t(apply(xcoef,1, cumsum))

xcoef=t(xcoef)
ycoef=t(ycoef)


for (TTT in (T_sample))
{
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

if (teststat > ql & teststat<qu){TEST_RESULT[iijj,which(TTT==T_sample)]=0}
if (teststat < ql){TEST_RESULT[iijj,which((TTT)==T_sample)]=-0.25}
if (teststat > qu){TEST_RESULT[iijj,which((TTT)==T_sample)]=0.25}



	if(teststat > qu){
	zz=zz[2:TTT,]-zz[1:(TTT-1),];zz0=zz0[2:TTT,]-zz0[1:(TTT-1),] 
	zz2=apply(zz0,2, cumsum)
	v0=lrvar(zz,kernel=2)$omega / TTT
	vx1 = crossprod(zz2) / (TTT^2)
	vx0 = lrvar(zz,kernel=2)$omega / TTT

	eval_vx0= t(ev0)%*%vx0%*%ev0
	eval_vx1= t(ev0)%*%vx1%*%ev0
	teststat=(eval_vx1/eval_vx0)

	if (teststat > ql & teststat<qu){TEST_RESULT2[iijj,which(TTT==T_sample)]=1}
	if (teststat < ql){TEST_RESULT2[iijj,which((TTT)==T_sample)]=0.75}
	if (teststat > qu){TEST_RESULT2[iijj,which((TTT)==T_sample)]=1.25}
	}else{TEST_RESULT2[iijj,which((TTT)==T_sample)]=-999}
	

}
if(seed_number%%100==0){print(seed_number)}
}

true_0=which(d_rand==0)
true_1=which(d_rand==1) 

correct_set1 = NULL;
correct_set2 = NULL;
correct_set3 = NULL;cs3=NULL
for (i in 1:length(T_sample))
{ 
cs1=NULL; cs2=NULL; cs3=NULL
for (j in 1:maxiter)
{
dett = (TEST_RESULT[j,i]==0) | (TEST_RESULT[j,i]!=0 & TEST_RESULT2[j,i]==1& TEST_RESULT2[j,i]!=-999)
dett2 = (TEST_RESULT[j,i]==0 & d_rand[j]==0)
dett3 = (TEST_RESULT[j,i]!=0 & TEST_RESULT2[j,i]==1 & TEST_RESULT2[j,i]!=-999 & d_rand[j]==1)
cs1=append(cs1,dett)
cs2=append(cs2,dett2)
cs3=append(cs3,dett3)
}
correct_set1=cbind(correct_set1,cs1)
correct_set2=cbind(correct_set2,cs2)
correct_set3=cbind(correct_set3,cs3)
}
round(colSums(correct_set1)/length(d_rand),digits=3)
round(colSums(correct_set2)/sum(d_rand==0),digits=3)
round(colSums(correct_set3)/sum(d_rand==1),digits=3)


