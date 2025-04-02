## Description: Critical values of the VR test can be computed from this code
## Specified quantiles for the max-test are stores in MAXCV
## Specified quantiles for the trace-test are stores in TRACECV

####################################################
## Key simulation paramters and related settings ###
####################################################
outeriter=200000      ## number of repetitions
etaquantile=0.95     ## quantile
leng=1000             ## sample points of the fractional Brownian motions
dim = 1:10           ## the set of value of q
library(sde)
####################################################
####################################################

MEAN=1

CV1=NULL
CV2=NULL

#SET=c(0.95)###################################################################


        aaa=NULL
		bbb=NULL
		
 for(nnn in 1:outeriter)
{
            atem=as.vector(BM(x=0, t0=0, T=1, N=leng-1))
            btem=as.vector(BBridge(x=0, t0=0, T=1, N=leng-1))
            CV1=append(CV1,mean(atem^2))
			CV2=append(CV2,mean(btem^2))
print(nnn)			
}

  
quantile(CV1,c(0.025,0.975))
quantile(CV2,c(0.025,0.975))




quantile_grid=quantile(CV1,seq(0,1,0.001))