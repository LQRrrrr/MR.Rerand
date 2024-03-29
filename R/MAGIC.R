#' Supplementary function for MAGIC to conduct rerandomization and Rao-Blackwellization
#' in GWAS (I) and GWAS (III).
#'
#'
#' Take GWAS (I) \eqn{(\hat{\beta}_{X_j}, {\sigma}_{X_j})} as an example, here
#' \eqn{\hat{\beta}_{X_j,\mathtt{RB}} =  \hat{\beta}_{X_j} - \frac{\sigma_{X_j}}{\eta}\Big(\phi\big(A_{j,+}\big) - \phi\big(A_{j,-}\big)\Big)\Big(\frac{\mathbf{1}_{(j\in\mathcal{S}_{\mathtt{x}})}}{1 - \Phi\big(A_{j,+}\big) + \Phi\big(A_{j,-}\big)} - \frac{\mathbf{1}_{(j\not\in\mathcal{S}_{\mathtt{x}})}}{ \Phi\big(A_{j,+}\big) - \Phi\big(A_{j,-}\big)} \Big)},
#' \eqn{\hat{\sigma}_{X_j,\mathtt{RB}}^{ \mathrm{2} }  = \sigma_{X_j}^2 \Bigg( 1 - \frac{1}{\eta^2}  \frac{A_{j,+}\phi(A_{j,+}) - A_{j,-}\phi(A_{j,-})}{1 - \Phi(A_{j,+}) + \Phi( A_{j,-} )} +\frac{1}{\eta^2} \Big(\frac{\phi(A_{j,+}) - \phi(A_{j,-}) }{1 - \Phi(A_{j,+}) + \Phi( A_{j, -} )}\Big)^2 \Bigg)},
#' and \eqn{\phi(\cdot)} is the density function of standard normal distrbution, \eqn{\Phi(\cdot)} is the cummulative function of standard normal distrbution.
#' @param gamma1.exp  SNP effect size's vector of the exposure vairable (GWAS(I))
#' @param gamma2.exp  SNP effect size's vector of the mediator vairable (GWAS(III))
#' @param se1.exp    SNP effect size's standard errors of \code{beta.exposure}
#' @param se2.exp    SNP effect size's standard errors of \code{beta.mediator}
#' @param etamean1  rerandomized scale of exposure variable. Default is 0.5.
#' @param etamean2  rerandomized scale of mediator variable. Default is 0.5.
#' @param pthr A vector of specified pre-screening threshold in the ordering of (exposure, mediator). Default is (5e-5, 5e-5). (corresponding \eqn{\lambda} is 4.06)
#' @param seed  The value of the random seed. Default is 0.
#' @return A list
#' \describe{
#' \item{filter1}{Indexs of selected relevant IVs in \eqn{\mathcal{S}_{\mathtt{x}}}}
#' \item{filter2}{Indexs of selected relevant IVs in \eqn{\mathcal{S}_{\mathtt{m}}}}
#' \item{gamma_exp1}{Effect size in GWAS (I) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{se1}{Standard errors in GWAS (I) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{gamma_exp2}{Effect size in GWAS (III) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{se2}{Standard errors in GWAS (III) after Rao-Blackwellization to eliminate the winner's curse}
#' \item{gamma_exp1.carve}{Effect size in GWAS (I) after Rao-Blackwellization to eliminate the loser's curse}
#' \item{gamma_exp2.carve}{Effect size in GWAS (III) after Rao-Blackwellization to eliminate the loser's curse}
#' }
#'
#' @references Rita Qiuran Lyu, Chong Wu, Xinwei Ma, Jingshen Wang (2023). Mediation Analysis with Mendelian Randomization and Efficient Multiple GWAS Integration. \url{https://arxiv.org/abs/2312.10563}.
#'
#' @import stats
#' @importFrom MASS mvrnorm
#' @export

pre_selection<-function(gamma1.exp,se1.exp,gamma2.exp,se2.exp, etamean1 = 0.5, etamean2=0.5,pthr = c(5e-5,5e-5),seed = 0) {

  p_x=pthr[1]
  p_m=pthr[2]
    set.seed(seed)
    #rerandomization
    W = mvrnorm(length(gamma1.exp), c(0,0), matrix(c(etamean1^2,0,0,etamean2^2),nrow=2,ncol=2,byrow=TRUE))
    # Step 2: Select significant SNPs:
    C_sel_x = qnorm(p_x/2,lower.tail = FALSE)
    C_sel_m = qnorm(p_m/2,lower.tail = FALSE)
    ind_filter1 = which(abs(gamma1.exp/se1.exp + W[,1]) >= C_sel_x)
    ind_filter2 = which(abs(gamma2.exp/se2.exp + W[,2]) >= C_sel_m)
    gamma1.exp_sel=gamma1.exp
    se1.exp_sel=se1.exp
    eta_sel=etamean1
    alpha1=(-C_sel_x-gamma1.exp_sel/se1.exp_sel)/eta_sel
    alpha2=(C_sel_x-gamma1.exp_sel/se1.exp_sel)/eta_sel
    gamma1.carve=gamma1.exp_sel-(se1.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha1)+1-pnorm(alpha2)))
    sigma21.carve=(1-((alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1))-((dnorm(alpha2)-dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1)))^2)/eta_sel^2)*se1.exp_sel^2

    warning("Some corrected IVs with variance <= 0\n")
       gamma1.carve2=gamma1.exp_sel+(se1.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))
    gamma2.exp_sel=gamma2.exp
    se2.exp_sel=se2.exp
    eta_sel=etamean2
    alpha1=(-C_sel_m-gamma2.exp_sel/se2.exp_sel)/eta_sel
    alpha2=(C_sel_m-gamma2.exp_sel/se2.exp_sel)/eta_sel
    gamma2.carve=gamma2.exp_sel-(se2.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha1)+1-pnorm(alpha2)))
    gamma2.carve2=gamma2.exp_sel+(se2.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))
    warning("Some corrected IVs with variance <= 0\n")
     sigma22.carve=(1-((alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1))-((dnorm(alpha2)-dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1)))^2)/eta_sel^2)*se2.exp_sel^2
  return(list(filter1=ind_filter1,filter2=ind_filter2,gamma_exp1=gamma1.carve,se1=sqrt(sigma21.carve),gamma_exp2=gamma2.carve,se2=sqrt(sigma22.carve),gamma_exp1.carve=gamma1.carve2,gamma_exp2.carve=gamma2.carve2))
}



#' Main function for using MAGIC framework to conduct mediation analysis.
#'
#' MAGIC stands for the Mediation Analysis framework through GWAS summary data Integration with the winner's (and the loser's) curse and measurement bias Correction.
#' Our framework efficiently integrates information stored in three independent GWAS summary data and mitigates the commonly encountered winner's curse and measurement error bias (a.k.a. instrument selection and weak instrument bias) in MR.
#' @param beta.exposure  SNP effect size's vector of the exposure vairable (GWAS(I))
#' @param beta.mediator  SNP effect size's vector of the mediator vairable (GWAS(III))
#' @param beta.outcome   SNP effect size's vector of the outcome vairable (GWAS(II))
#' @param se.exposure    SNP effect size's standard errors of \code{beta.exposure}
#' @param se.mediator    SNP effect size's standard errors of \code{beta.mediator}
#' @param se.outcome     SNP effect size's standard errors of \code{beta.outcome}
#' @param Conf.level Confidence level. Default is 0.95.
#' @param pval.select A vector of specified pre-screening threshold in the ordering of (exposure, mediator). Default is (5e-5, 5e-5). (corresponding \eqn{\lambda} is 4.06)
#' @param eta A vector of rerandomized scale in the ordering of  (exposure, mediator). Default is (0.5,0.5).
#' @param seed The value of the random seed. Default is 0.
#' @return A list
#' \describe{
#' \item{theta.hat}{Estimated direct effect from exposure to outcome variable}
#' \item{tauy.hat}{Estimated direct effect from mediator to outcome variable}
#' \item{taux.hat}{Estimated indirect effect from exposure to mediator variable}
#' \item{tau.hat}{Estimated mediation effect}
#' \item{tau_total.hat}{Estimated total effect}
#' \item{theta.se}{Standard error of \code{theta.hat}}
#' \item{tauy.se}{Standard error of \code{tauy.hat}}
#' \item{taux.se}{Standard error of \code{taux.hat}}
#' \item{tau.se}{Standard error of \code{tau.hat}}
#' \item{tau_total.se}{Standard error of \code{tau_total.hat}}
#' \item{n.IV.exp}{Number of IVs used in exposure dataset}
#' \item{n.IV.med}{Number of IVs used in mediator dataset}
#' \item{Conf.Interval}{Confidence interval given \code{Conf.level}}
#' \item{IV.exp}{The index of IVs selected in \eqn{\mathcal{S}_{\mathtt{x}}}}
#' \item{IV.med}{The index of IVs selected in \eqn{\mathcal{S}_{\mathtt{m}}}}
#' }
#'
#' @references Rita Qiuran Lyu, Chong Wu, Xinwei Ma, Jingshen Wang (2023). Mediation Analysis with Mendelian Randomization and Efficient Multiple GWAS Integration. \url{https://arxiv.org/abs/2312.10563}.
#'
#' @import stats
#' @importFrom MASS mvrnorm
#' @importFrom msm deltamethod
#' @export
#'
MAGIC<-function(beta.exposure, beta.mediator, beta.outcome, se.exposure, se.mediator, se.outcome, Conf.level=0.95, pval.select=c(5e-5,5e-5), eta=c(0.5,0.5),seed=0)
{
  eta_x=eta[1]
  eta_m=eta[2]
  pval_x=pval.select[1]
  pval_m=pval.select[2]
  # selection with rerandomization - no winner's curse
  over_summary=pre_selection(gamma1.exp = beta.exposure, se1.exp = se.exposure,gamma2.exp =beta.mediator,se2.exp = se.mediator,etamean1 = eta_x,etamean2 = eta_m,pthr=pval.select,seed=seed)
  #rerand
  # ********Proposed MR***********
  # *****Rerandomization****
  indx=over_summary$filter1
  indm=over_summary$filter2
  overlap_index=intersect(indx,indm)#overlapping snps (for covariance estimation)
  overlap_m=which(indm%in%overlap_index) # the order in selected snps Sm
  overlap_x=which(indx%in%overlap_index) #new beta and se

  betanew_x=over_summary$gamma_exp1[indx]
  betanew_m=over_summary$gamma_exp2[indm]
  senew_x=over_summary$se1[indx]
  senew_m=over_summary$se2[indm]
  ele1=sum((betanew_x^2-senew_x^2)/(se.outcome[indx]^2))
  ele4=sum((betanew_m^2-senew_m^2)/(se.outcome[indm]^2))
  ele5=sum((betanew_x^2-senew_x^2)/(se.mediator[indx]^2))
  ele2=sum(over_summary$gamma_exp1[overlap_index]*over_summary$gamma_exp2[overlap_index]/(se.outcome[overlap_index]^2))+sum(over_summary$gamma_exp1[setdiff(indx,overlap_index)]*over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)]/(se.outcome[setdiff(indx,overlap_index)]^2))
  ele3=sum(over_summary$gamma_exp1[overlap_index]*over_summary$gamma_exp2[overlap_index]/(se.outcome[overlap_index]^2))+sum(over_summary$gamma_exp1.carve[setdiff(indm,overlap_index)]*over_summary$gamma_exp2[setdiff(indm,overlap_index)]/(se.outcome[setdiff(indm,overlap_index)]^2))
  eleY1=sum(betanew_x*beta.outcome[indx]/(se.outcome[indx]^2))
  eleY2=sum(betanew_m*beta.outcome[indm]/(se.outcome[indm]^2))
  eleY3=sum(betanew_x[overlap_x]*over_summary$gamma_exp2[overlap_index]/(se.mediator[overlap_index]^2))+sum(over_summary$gamma_exp1[setdiff(indx,overlap_index)]*over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)]/(se.mediator[setdiff(indx,overlap_index)]^2))

  sum_mat=matrix(c(ele1,ele2,0,ele3,ele4,0,0,0,ele5),nrow=3,ncol=3,byrow=TRUE)
  sum_vector=c(eleY1,eleY2,eleY3)
  estimator_rerand=solve(sum_mat)%*%sum_vector
  thetaRIVW=estimator_rerand[1]
  tauyRIVW=estimator_rerand[2]
  tauxRIVW=estimator_rerand[3]
  tauRIVW=tauyRIVW*tauxRIVW
  tau_totalRIVW=tauRIVW+thetaRIVW
  regression_residuals1= sum( (beta.outcome[overlap_index] *over_summary$gamma_exp1[overlap_index]-thetaRIVW * (over_summary$gamma_exp1[overlap_index]^2 - over_summary$se1[overlap_index]^2) -tauyRIVW*over_summary$gamma_exp1[overlap_index]*over_summary$gamma_exp2[overlap_index])^2/ se.outcome[overlap_index]^4)+sum( (beta.outcome[setdiff(indx,overlap_index)] * over_summary$gamma_exp1[setdiff(indx,overlap_index)]-thetaRIVW * ( over_summary$gamma_exp1[setdiff(indx,overlap_index)]^2 -  over_summary$se1[setdiff(indx,overlap_index)]^2) -tauyRIVW*over_summary$gamma_exp1[setdiff(indx,overlap_index)]*over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)])^2/ se.outcome[setdiff(indx,overlap_index)]^4)
  regression_residuals2= sum( (beta.outcome[overlap_index] *over_summary$gamma_exp2[overlap_index]-tauyRIVW * (over_summary$gamma_exp2[overlap_index]^2 - over_summary$se2[overlap_index]^2) -thetaRIVW*over_summary$gamma_exp2[overlap_index]*over_summary$gamma_exp1[overlap_index])^2/ se.outcome[overlap_index]^4)+sum( (beta.outcome[setdiff(indm,overlap_index)] * over_summary$gamma_exp2[setdiff(indm,overlap_index)]-tauyRIVW * ( over_summary$gamma_exp2[setdiff(indm,overlap_index)]^2 -  over_summary$se2[setdiff(indm,overlap_index)]^2) -thetaRIVW*over_summary$gamma_exp2[setdiff(indm,overlap_index)]*over_summary$gamma_exp1.carve[setdiff(indm,overlap_index)])^2/ se.outcome[setdiff(indm,overlap_index)]^4)
  regression_residuals3= sum( (betanew_m[overlap_m] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2))^2 / se.mediator[overlap_index]^4)+sum( (over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)] * over_summary$gamma_exp1[setdiff(indx,overlap_index)] - tauxRIVW * (over_summary$gamma_exp1[setdiff(indx,overlap_index)]^2 - over_summary$se1[setdiff(indx,overlap_index)]^2))^2 / se.mediator[setdiff(indx,overlap_index)]^4)

  cov_12= sum( (beta.outcome[overlap_index] * betanew_x[overlap_x]-thetaRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2) -tauyRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(beta.outcome[overlap_index] * betanew_m[overlap_m] -tauyRIVW * (betanew_m[overlap_m]^2 - senew_m[overlap_m]^2) -thetaRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])/ se.outcome[overlap_index]^4)
  cov_13=sum( (beta.outcome[overlap_index] * betanew_x[overlap_x]-thetaRIVW * (betanew_x[overlap_x]^2 - senew_x[overlap_x]^2) -tauyRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(over_summary$gamma_exp2[overlap_index] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2)) / (se.outcome[overlap_index]^2*se.mediator[overlap_index]^2))+sum( (beta.outcome[setdiff(indx,overlap_index)] * over_summary$gamma_exp1[setdiff(indx,overlap_index)]-thetaRIVW * (over_summary$gamma_exp1[setdiff(indx,overlap_index)]^2 - over_summary$se1[setdiff(indx,overlap_index)]^2) -tauyRIVW*over_summary$gamma_exp1[setdiff(indx,overlap_index)]*over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)])*(over_summary$gamma_exp2.carve[setdiff(indx,overlap_index)] * over_summary$gamma_exp1[setdiff(indx,overlap_index)] - tauxRIVW * (over_summary$gamma_exp1[setdiff(indx,overlap_index)]^2 - over_summary$gamma_exp1[setdiff(indx,overlap_index)]^2))/ (se.outcome[setdiff(indx,overlap_index)]^2*se.mediator[setdiff(indx,overlap_index)]^2))
  cov_23=sum((beta.outcome[overlap_index] * betanew_m[overlap_m] -tauyRIVW * (betanew_m[overlap_m]^2 - senew_m[overlap_m]^2) -thetaRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(over_summary$gamma_exp2[overlap_index] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2))/(se.outcome[overlap_index]^2*se.mediator[overlap_index]^2))

  var_matrix=matrix(c(regression_residuals1,cov_12,cov_13,
                      cov_12,regression_residuals2,cov_23,
                      cov_13,cov_23,regression_residuals3),nrow=3,ncol=3,byrow=TRUE)
  var_matrix_final=solve(sum_mat)%*%var_matrix%*%t(solve(sum_mat))
  v1_theta=var_matrix_final[1,1]
  v2_tauy=var_matrix_final[2,2]
  v3_taux=var_matrix_final[3,3]
  #mediation effect
  vector_C=c(tauyRIVW,tauxRIVW)
  sub_variance_matrix=var_matrix_final[2:3,2:3]
  sd_tau=msm::deltamethod(~ x1*x2,vector_C,sub_variance_matrix)
  #total effect
  vector_C=c(thetaRIVW,tauyRIVW,tauxRIVW)
  sub_variance_matrix=var_matrix_final[1:3,1:3]
  sd_tau_total=msm::deltamethod(~ x1+x2*x3,vector_C,sub_variance_matrix)
  theta.hat=thetaRIVW
  tauy.hat=tauyRIVW
  taux.hat=tauxRIVW
  tau.hat=tauRIVW
  tau_total.hat=tau_totalRIVW
  theta.se=sqrt(v1_theta)
  tauy.se=sqrt(v2_tauy)
  taux.se=sqrt(v3_taux)
  tau.se=sd_tau
  tau_total.se=sd_tau_total
  n.IV.exp=length(betanew_x)
  n.IV.med=length(betanew_m)
  IV.exp=indx
  IV.med=indm
  z_score=qnorm((1-Conf.level)/2,lower.tail = FALSE)
  CI_theta=c(theta.hat-z_score*theta.se,theta.hat+z_score*theta.se)
  CI_tauy=c(tauy.hat-z_score*tauy.se,tauy.hat+z_score*tauy.se)
  CI_taux=c(taux.hat-z_score*taux.se,taux.hat+z_score*taux.se)
  CI_tau=c(tau.hat-z_score*tau.se,tau.hat+z_score*tau.se)
  CI_tau_total=c(tau_total.hat-z_score*tau_total.se,tau_total.hat+z_score*tau_total.se)
  tempdat=as.data.frame(rbind(CI_theta,CI_tauy, CI_taux, CI_tau, CI_tau_total))
  rownames(tempdat)=c("theta","tauy","taux","tau","tau.total")
  colnames(tempdat)=c("Lower Bound", "Upper Bound")
  Conf.Interval=tempdat
  return(list(theta.hat=theta.hat,theta.se=theta.se,tauy.hat=tauy.hat,tauy.se=tauy.se,taux.hat=taux.hat,taux.se=taux.se,tau.hat=tau.hat,tau.se=tau.se,tau_total.hat=tau_total.hat,
              tau_total.se=tau_total.se,n.IV.exp=n.IV.exp,n.IV.med=n.IV.med,IV.exp=IV.exp,IV.med=IV.med,Conf.Interval=Conf.Interval))

}
