#' Supplementary function for RIVW
#'
#' @param gamma1.exp  SNP effect size's vector of the exposure vairable
#' @param se1.exp    SNP effect size's standard errors of \code{beta.exposure}
#' @param etamean  rerandomized scale of exposure variable. Default is 0.5.
#' @param pthr The specified pre-screening threshold. Default is 5e-5. (corresponding lambda is 4.06)
#' @param seed  A random seed. Default is 0.
#' @param smoothing Whether to use smoothing to decrease variance . Default is FALSE.
#' @return A list
#' \describe{
#' \item{filter1}{Indexs of selected relevant IVS.}
#' \item{gamma_exp1}{Effect size in GWAS (I) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{se1}{Standard errors in GWAS (I) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{weights}{The weights for each SNP. If smoothing is False, weights are the same for each SNP.}
#' }
#'
#' @references   Xinwei Ma, Jingshen Wang, Chong Wu. (2023). Breaking the Winner’s Curse in Mendelian Randomization:Rerandomized Inverse Variance Weighted Estimator \url{https://projecteuclid.org/journals/annals-of-statistics/volume-51/issue-1/Breaking-the-winners-curse-in-Mendelian-randomization--Rerandomized-inverse/10.1214/22-AOS2247.full}.
#'
#' @import stats
#'
#' @export


pre_screening<-function(gamma1.exp,se1.exp, etamean = 0.5,pthr = 5e-5,seed = 0, smoothing=FALSE) {
  set.seed(seed)
  if(smoothing==FALSE)
  {p_x=pthr

  #rerandomization
  W = rnorm(length(gamma1.exp), 0, etamean)
  # Step 2: Select significant SNPs:
  C_sel_x = qnorm(p_x/2,lower.tail = FALSE)
  ind_filter1 = which(abs(gamma1.exp/se1.exp + W) >= C_sel_x)
  gamma1.exp_sel=gamma1.exp[ind_filter1]
  se1.exp_sel=se1.exp[ind_filter1]
  eta_sel=etamean
  alpha1=(-C_sel_x-gamma1.exp_sel/se1.exp_sel)/eta_sel
  alpha2=(C_sel_x-gamma1.exp_sel/se1.exp_sel)/eta_sel
  gamma1.carve=gamma1.exp_sel-(se1.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha1)+1-pnorm(alpha2)))
  sigma21.carve=(1-((alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1))-((dnorm(alpha2)-dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1)))^2)/eta_sel^2)*se1.exp_sel^2
  weights=rep(1,length(ind_filter1))


  }else
    {
      W = rnorm(length(gamma1.exp), 0, etamean)

      # Step 2: Select significant SNPs:
      C_sel = qnorm(pthr/2,lower.tail = FALSE)

      ind_filter = which(abs(gamma1.exp/se1.exp + W) >= 0) # there is no selections, we used all IVs in sRIVW

      gamma.exp_sel = gamma1.exp
      se.exp_sel = se1.exp

      # Step 3. Construct the unbiased carved estimator (also the UMVUE)
      alpha1 = (-C_sel - gamma.exp_sel/se.exp_sel) / etamean
      alpha2 = (C_sel - gamma.exp_sel/se.exp_sel) / etamean

      weights =  pnorm(alpha1) + 1 - pnorm(alpha2)

      gamma.carve = gamma.exp_sel - (se.exp_sel/etamean) * ( (dnorm(alpha2) - dnorm(alpha1)) / weights )
      sigma2.carve = (1 - ((alpha2*dnorm(alpha2) - alpha1*dnorm(alpha1)) / (1 - pnorm(alpha2) + pnorm(alpha1) ) - ((dnorm(alpha2) - dnorm(alpha1))/(1 - pnorm(alpha2) + pnorm(alpha1)))^2) / etamean^2 ) * se.exp_sel^2

      # For the nu
      warning("Some corrected IVs with variance <= 0\n")
      tmp = (1 + etamean^2 ) * se.exp_sel^2
      sigma2.carve[sigma2.carve <=0] = tmp[sigma2.carve <=0]
  }

  return(list(filter1=ind_filter1,gamma_exp1=gamma1.carve,se1=sqrt(sigma21.carve),weights=weights))
}




#' Main function for RIVW
#' @param beta.exposure  SNP effect size's vector of the exposure vairable (GWASI)
#' @param beta.outcome   SNP effect size's vector of the outcome vairable (GWASII)
#' @param se.exposure    SNP effect size's standard errors of \code{beta.exposure}
#' @param se.outcome     SNP effect size's standard errors of \code{beta.outcome}
#' @param Conf.level Confidence level. Default is 0.95.
#' @param smoothing  Whether to use smoothing to decrease variance. Default is FALSE.
#' @param pval.select The specified pre-screening threshold. Default is 5e-5. (corresponding lambda is 4.06)
#' @param eta A vector of rerandomized scale. Default is 0.5.
#' @param seed The value of random seed. Default is 0.
#' @return A list
#' \describe{
#' \item{beta.rerand}{Exposure dataset effect size after rerandomization.}
#' \item{se.rerand}{Exposure dataset standard errors after rerandomization.}
#' \item{beta.hat}{Estimated direct effect from exposure to outcome variable}
#' \item{beta.se}{Standard error of \code{beta.hat}}
#' \item{IV}{The index of IVs selected in Sx}
#' \item{n.IV}{Number of IVs used in exposure dataset}
#' \item{F}{The value of F-statistic}
#' \item{p.val}{The p-value of estimated causal effect}
#' \item{Conf.Interval}{Confidence interval of the causal effect given \code{Conf.level}}
#' }
#'
#' @references   Xinwei Ma, Jingshen Wang, Chong Wu. (2023). Breaking the Winner’s Curse in Mendelian Randomization:Rerandomized Inverse Variance Weighted Estimator \url{https://projecteuclid.org/journals/annals-of-statistics/volume-51/issue-1/Breaking-the-winners-curse-in-Mendelian-randomization--Rerandomized-inverse/10.1214/22-AOS2247.full}.
#'
#' @import stats
#'
#' @export
RIVW<-function(beta.exposure, beta.outcome, se.exposure, se.outcome, Conf.level=0.95, smoothing=FALSE, pval.select=5e-5,eta=0.5,seed=0)
{
  over_summary=pre_screening(gamma1.exp = beta.exposure, se1.exp=se.exposure,etamean = eta, pthr=pval.select,seed=seed,smoothing=smoothing)
  # *****Rerandomization****
  ind_filter=over_summary$filter1
  betanew_x=over_summary$gamma_exp1
  senew_x=over_summary$se1

  if(smoothing==FALSE)
  {
  beta.out_sel=beta.outcome[ind_filter]
  se.out_sel=se.outcome[ind_filter]
  beta = sum(betanew_x*beta.out_sel*(1/se.out_sel^2)) / sum((betanew_x^2 - senew_x^2)/se.out_sel^2)
  numIV_sel = length(ind_filter)
  RIVW.var = sum( (beta.out_sel * betanew_x - beta * (betanew_x^2 - senew_x^2) )^2 / se.out_sel^4) / (sum((betanew_x^2 - senew_x^2) / se.out_sel^2) )^2
  beta.se = sqrt(RIVW.var)
  F_stats = sum(beta.exposure[ind_filter]^2/se.exposure[ind_filter]^2)/length(beta.exposure[ind_filter]) - 1
  p = pnorm(abs(beta/beta.se), lower.tail = F) * 2
  }else
  {
    weights=over_summary$weights
    beta.out_sel=beta.outcome[ind_filter]
    se.out_sel=se.outcome[ind_filter]
    beta = sum(betanew_x * beta.out_sel * weights * (1/se.out_sel^2)) / sum((betanew_x^2 - senew_x^2) * weights / se.out_sel^2)
    # estimation based on regression residuals
    RIVW.var = sum( (beta.out_sel * betanew_x - beta * (betanew_x^2 - senew_x^2) )^2 * weights^2 / se.out_sel^4) / (sum((betanew_x^2 -senew_x^2) * weights / se.out_sel^2) )^2
    beta.se = sqrt(RIVW.var)
    F_stats = sum(beta.exposure[ind_filter]^2/se.exposure[ind_filter]^2)/length(beta.exposure[ind_filter]) - 1
    p = pnorm(abs(beta/beta.se), lower.tail = F) * 2
  }
  beta.rerand=betanew_x
  se.rerand=senew_x
  beta.hat=beta
  beta.se=beta.se
  n.IV=length(ind_filter)
  p.val=p
  z_score=qnorm((1-Conf.level)/2,lower.tail = FALSE)
  Conf.Interval=c(beta.hat-z_score*beta.se,beta.hat+z_score*beta.se)
  return(list(beta.rerand=beta.rerand,se.rerand=se.rerand,beta.hat=beta.hat,beta.se=beta.se,IV=ind_filter,n.IV=n.IV,F=F_stats,p.val=p.val,Conf.Interval=Conf.Interval))


  }
