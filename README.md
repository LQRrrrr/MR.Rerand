# MR.Rerand: Re-randomized Inverse-Variance Weighted Estimator in two-sample Mendelian Randomization with Summary-data and Mediation analysis in Mendelian Randomization with Summary-data


To install this package in R, run the following commands:


```R
library(devtools) 
install_github("LQRrrrr/MR.Rerand")
```

This package can implement (i) RIVW for estimating the causal effect $\beta$ from the exposure variable $X$ to the outcome variable $Y$, (ii)MAGIC for estimating the direct effects $\theta$ and $\tau_Y$, the mediation effect $\tau=\tau_X\tau_Y$, the causal effect $\tau_X$ from $X$ to $M$, and the total effect $\theta+\tau$ as parameters shown in the below figure.

Another important usage of this package is to generate the re-randomized effect size and the standard errors after Rao-Blackwellization for two-sample MR analysis to eliminate the winner's curse bias. The generated effect sizes and the standard errors can then be applied to two-sample MR methods with summary data.  See more details in reference 1.

- MAGIC implements MAGIC for mediation analysis using Mendelian Randomization with summary data
- RIVW implements RIVW for two-sample Mendelian Randomization with summary data

The effect sizes and standard errors after Rao-Blackwellization can be obtained by

- pre_screening rerandomizing in the exposure dataset and conduct Rao-Blackwellization to eliminate the winner's curse and the loser's curse

- pre_selection rerandomizing in both the exposure dataset and the mediator dataset and conduct Rao-Blackwellization to eliminate both the winner's curse and the loser's curse



<figure>
  <img src="man/figures/causal_diagram.jpeg" alt="Example Image">
  <figcaption>Causal Diagram for Mediation Analysis.</figcaption>
</figure>

## Example usage:

```R
library(MR.Rerand)
##Generating data
   M=100000
   pi1=pi2=0.02
   pi3=0.02
   nx=ny=nm=100000
   se_x=rep(sqrt(1/nx),M)
   se_y=rep(sqrt(1/ny),M)
   se_m=rep(sqrt(1/nm),M)
   sigma2x=sigma2y=0.5e-4
   theta=-0.2
   tauY=0.2
   tauX=0.6
  gamma=rep(0,M)
  ind1=sample(M,round(M*pi1))#valid for betaXj
  causalsnps=ind1
  ## simulation 1 (see more simulation settings in reference 2)
  ind3=sample(setdiff(1:M,causalsnps),round(M*pi3))# valid for sj
  causalsnps=c(causalsnps,ind3)
  alpha=rep(0,M)
  gamma[ind1]=rnorm(length(ind1),0,sd=sqrt(sigma2x))
  alpha[ind3]=rnorm(length(ind3),0,sd=sqrt(sigma2y))
  gammam=tauX*gamma+alpha# pleiotropy sj
  gammay=theta*gamma+tauY*gammam#no pleiotropy
      betax=gamma
      betam=gammam
      betay=gammay
      betahat_x=gamma+rnorm(M,mean=0,sd=sqrt(1/nx))
      betahat_m=gammam+rnorm(M,mean=0,sd=sqrt(1/nm))
      betahat_y=gammay+rnorm(M,mean=0,sd=sqrt(1/ny))
  ## Conducting mediation analysis  
  mediaiton.result=MAGIC(betahat_x,betahat_m,betahat_y,se_x,se_m,se_y)
  ## Conducting two-sample MR using RIVW with smoothing to decrease the variance
  RIVW.result,smoothing=RIVW(betahat_x,betahat_m,se_x,se_m,smoothing = TRUE)   
  ## Conducting two-sample MR using RIVW with smoothing to decrease the variance
  RIVW.result=RIVW(betahat_x,betahat_m,se_x,se_m)
  ## Just generate rerandomized effect size and standard errors
  rerandomized.GWAS=pre_screening(betahat_x,se_x)
```

### References

- Ma, X., Wang, J., and Wu, C. (2023a). “Breaking the winner’s curse in Mendelian randomization: rerandomized inverse variance weighted estimator,” Annals of Statistics, 51 (1), 211–232.

- Lyu, R. Q., Wu, C., Ma, X., and Wang, J. (2023). Mediation Analysis with Mendelian Randomization and Efficient Multiple GWAS Integration. arXiv preprint arXiv:2312.10563.
