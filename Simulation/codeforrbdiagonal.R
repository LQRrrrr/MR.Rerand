#no confounding
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(car)
rm(list = ls())
cover<-function(theta_rerandlist,theta_rerandse,theta)
{upper_theta=theta_rerandlist+1.96*theta_rerandse
lower_theta=theta_rerandlist-1.96*theta_rerandse
cover_theta=sum((upper_theta>=theta)&(lower_theta<=theta))/length(theta_rerandlist)
return(cover_theta)
}
abias<-function(theta_rerandlist,theta)
{

  return(abs(mean(theta_rerandlist)-theta))
}
pre_selection<-function(gamma1.exp,se1.exp,gamma2.exp,se2.exp, etamean1 = 0.5, etamean2=0.5,pthr = 5e-8,rerand=FALSE,seed = sample(1:100000,1), lambda = 0) {

  if (rerand==TRUE)
  {set.seed(seed)
    #rerandomization
    W = mvrnorm(length(gamma1.exp), c(0,0), matrix(c(etamean1^2,0,0,etamean2^2),nrow=2,ncol=2,byrow=TRUE))
    # Step 2: Select significant SNPs:
    C_sel = qnorm(pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    ind_filter1 = which(abs(gamma1.exp/se1.exp + W[,1]) >= C_sel)
    ind_filter2 = which(abs(gamma2.exp/se2.exp + W[,2]) >= C_sel)
    gamma1.exp_sel=gamma1.exp
    se1.exp_sel=se1.exp
    eta_sel=etamean1
    alpha1=(-C_sel-gamma1.exp_sel/se1.exp_sel)/eta_sel
    alpha2=(C_sel-gamma1.exp_sel/se1.exp_sel)/eta_sel
    gamma1.carve=gamma1.exp_sel-(se1.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha1)+1-pnorm(alpha2)))
    sigma21.carve=(1-((alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1))-((dnorm(alpha2)-dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1)))^2)/eta_sel^2)*se1.exp_sel^2
    gamma1.carve2=gamma1.exp_sel+(se1.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))


    gamma2.exp_sel=gamma2.exp
    se2.exp_sel=se2.exp
    eta_sel=etamean2
    alpha1=(-C_sel-gamma2.exp_sel/se2.exp_sel)/eta_sel
    alpha2=(C_sel-gamma2.exp_sel/se2.exp_sel)/eta_sel
    gamma2.carve=gamma2.exp_sel-(se2.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha1)+1-pnorm(alpha2)))
    gamma2.carve2=gamma2.exp_sel+(se2.exp_sel/eta_sel)*((dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))
    sigma22.carve=(1-((alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1))-((dnorm(alpha2)-dnorm(alpha1))/(1-pnorm(alpha2)+pnorm(alpha1)))^2)/eta_sel^2)*se2.exp_sel^2





  }
  else{
    C_sel = qnorm(pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    ind_filter1 = which(abs(gamma1.exp/se1.exp ) >= C_sel)
    gamma1.exp_sel=gamma1.exp[ind_filter1]
    se1.exp_sel=se1.exp[ind_filter1]
    gamma1.carve=gamma1.exp_sel
    sigma21.carve=se1.exp_sel^2
    gamma1.carve2=gamma1.exp_sel
    ind_filter2 = which(abs(gamma2.exp/se2.exp ) >= C_sel)
    gamma2.exp_sel=gamma2.exp[ind_filter2]
    se2.exp_sel=se2.exp[ind_filter2]
    gamma2.carve=gamma2.exp_sel
    sigma22.carve=se2.exp_sel^2
    gamma2.carve2=gamma2.exp_sel
  }
  return(list(filter1=ind_filter1,filter2=ind_filter2,gamma_exp1=gamma1.carve,se1=sqrt(sigma21.carve),gamma_exp2=gamma2.carve,se2=sqrt(sigma22.carve),gamma_exp1.carve=gamma1.carve2,gamma_exp2.carve=gamma2.carve2))
}






## oracle selection
##simulation setting 1 not overlapping 2. partially overlapping 3.the same region

theta_list=c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2)
tauY_list=c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2)

tauX_list=c(0,0.6)
M=100000
pi1_list=c(0.005,0.01,0.015)
pi3_list=c(0.01,0.015)
nx=ny=nm=100000
se_x=rep(sqrt(1/nx),M)
se_y=rep(sqrt(1/ny),M)
se_m=rep(sqrt(1/nm),M)
sigma2x=1e-4
sigma2y=0.5e-4
sigma2u=sigma2x_td=sigma2y_td=1e-4
simulation=c(1,2,3)
set.seed(10000)
for(tauX in tauX_list)
{for (theta in theta_list)
{
  for(tauY in tauY_list)
  {
    for(pi1 in pi1_list)
     {for(pi3 in pi3_list)
     {for (design in simulation)

      {tau=tauX*tauY


      twostep_tauY=c()
      twostep_tauX=c()
      twostep_tau=c()

      twostep_tauY_sd=c()
      twostep_tauX_sd=c()
      twostep_tau_sd=c()
      total_effects_list=c()
      difference_tau_list=c()
      total_effects_sd=c()
      mvmr_theta=c()
      mvmr_tauY=c()
      mvmr_tauX=c()
      mvmr_tau=c()
      mvmr_theta_ivw=c()
      mvmr_tauY_ivw=c()
      mvmr_tauX_ivw=c()
      mvmr_tau_ivw=c()
      #thetalist=c()
      ivw_theta_se=c()
      #tauYlist=c()
      ivw_tauY_se=c()
      #tauXlist=c()
      ivw_tauX_se=c()
      ivw_tau_se=c()
      theta_rerandlist=c()
      theta_rerandse=c()
      tauY_rerandlist=c()
      tauY_rerandse=c()
      tauX_rerandlist=c()
      tauX_rerandse=c()
      tau_rerandlist=c()
      tau_rerandse=c()

      # we use the same direct effect across all settings
      ele2_original=c()
      ele2_rb=c()
      print(paste0("theta=",theta,"tauX=",tauX,"tauY=",tauY))
      index_ratio=c()
      index_sx=c()
      index_sm=c()
      for (i in 1:1000)
      {
      print(i)
        ind1=sample(M,round(M*pi1))#valid for betaXj
        causalsnps=ind1
        gamma=rep(0,M)
        alpha=rep(0,M)
        alpha2=rep(0,M)

        # ****Simulation1****
        if(design==1)
       { ind3=sample(setdiff(1:M,causalsnps),round(M*pi3))# valid for sj
        causalsnps=c(causalsnps,ind3)

        alpha[ind3]=rnorm(length(ind3),0,sd=sqrt(sigma2y))
        }
        #ind2=sample(setdiff(1:M,causalsnps),round(M*pi3))#valid
        # ****Simulation2****
        if(design==2)
        {temp=sample(length(ind1),0.5*length(ind1))# overlappin=ng part
        ind3=c(ind1[temp],sample(setdiff(1:M,causalsnps),round(M*pi3-0.5*length(ind1))))
        causalsnps=c(causalsnps,ind3)
        alpha[ind3]=rnorm(length(ind3),0,sd=sqrt(sigma2y))
        }
        if(design==3)
        {
          if(pi1!=pi3) {next}else{
          if(M*pi3<M*pi1)
          {temp=sample(length(ind1),M*pi3)
          ind3=ind1[temp]
          causalsnps=c(causalsnps,ind3)
          alpha[ind3]=rnorm(length(ind3),0,sd=sqrt(sigma2y))}else{
            temp=sample(length(ind1),length(ind1))
            # overlappin=ng part
            ind3=c(ind1[temp],sample(setdiff(1:M,causalsnps),round(M*pi3-1*length(ind1))))
            causalsnps=c(causalsnps,ind3)
            alpha[ind3]=rnorm(length(ind3),0,sd=sqrt(sigma2y))
          }}
        }
        gamma[ind1]=rnorm(length(ind1),0,sd=sqrt(sigma2x))

        #alpha2[ind2]=rnorm(length(ind2),0,sqrt(sigma2y))

        gammam=tauX*gamma+alpha# pleiotropy sj
        gammay=theta*gamma+tauY*gammam#no pleiotropy
        ratiom=which(gammam!=0)
        ratiox=which(gamma!=0)
        index_sx=c(index_sx,length(setdiff(ratiox,ratiom)))
        index_sm=c(index_sm,length(setdiff(ratiom,ratiox)))
        ##perfect selection
        ##
        betax=gamma
        betam=gammam
        betay=gammay
        betahat_x=betax+rnorm(M,mean=0,sd=sqrt(1/nx))
        betahat_m=betam+rnorm(M,mean=0,sd=sqrt(1/nm))
        betahat_y=betay+rnorm(M,mean=0,sd=sqrt(1/ny))
        se_x=rep(sqrt(1/nx),M)
        se_y=rep(sqrt(1/ny),M)
        se_m=rep(sqrt(1/nm),M)
        pthr=5e-5
        etamean=0.5
        C_sel = qnorm(pthr/2,lower.tail = FALSE)
        # selection without rerandomization -winner's curse
        over1=pre_selection(gamma1.exp = betahat_x,se1.exp = se_x,gamma2.exp =betahat_m,se2.exp = se_m,etamean1 = 0.5,etamean2 = 0.5,pthr=5e-8,rerand = FALSE )
        # selection with rerandomization - no winner's curse
        over2=pre_selection(gamma1.exp = betahat_x,se1.exp = se_x,gamma2.exp =betahat_m,se2.exp = se_m,etamean1 = 0.5,etamean2 = 0.5,pthr=pthr,rerand =TRUE )

        # ********MVMR IVW corrected measurement error***********
        indx_mvmr=over1$filter1
        indm_mvmr=over1$filter2
        ind_union=union(indx_mvmr,indm_mvmr)
        overlap_index=intersect(indx_mvmr,indm_mvmr)
        indx_single=setdiff(indx_mvmr,indm_mvmr)
        indm_single=setdiff(indm_mvmr,indx_mvmr)

        #no rerand but correct measurement error selection set is the union set with pre-screening
        ele1_mvmr=sum((betahat_x[ind_union]^2-se_x[ind_union]^2)/(se_y[ind_union]^2))
        ele2_mvmr=sum(betahat_x[ind_union]*betahat_m[ind_union]/(se_y[ind_union]^2))
        ele3_mvmr=sum(betahat_x[ind_union]*betahat_m[ind_union]/(se_y[ind_union]^2))
        ele4_mvmr=sum((betahat_m[ind_union]^2-se_m[ind_union]^2)/(se_y[ind_union]^2))
        ele5_mvmr=sum((betahat_x[indx_mvmr]^2-se_x[indx_mvmr]^2)/(se_m[indx_mvmr]^2))
        eleY1_mvmr=sum(betahat_x[ind_union]*betahat_y[ind_union]/(se_y[ind_union]^2))
        eleY2_mvmr=sum(betahat_m[ind_union]*betahat_y[ind_union]/(se_y[ind_union]^2))
        eleY3_mvmr=sum(betahat_m[indx_mvmr]*betahat_x[indx_mvmr]/(se_m[indx_mvmr]^2))
        sum_mat_mvmr=matrix(c(ele1_mvmr,ele2_mvmr,0,ele3_mvmr,ele4_mvmr,0,0,0,ele5_mvmr),nrow=3,ncol=3,byrow=TRUE)
        sum_vector_mvmr=c(eleY1_mvmr,eleY2_mvmr,eleY3_mvmr)
        estimator_mvmr=solve(sum_mat_mvmr)%*%sum_vector_mvmr
        theta_mvmr=estimator_mvmr[1]
        tauy_mvmr=estimator_mvmr[2]
        taux_mvmr=estimator_mvmr[3]
        tau_mvmr=tauy_mvmr*taux_mvmr

        mvmr_theta=c(mvmr_theta,theta_mvmr)
        mvmr_tauY=c(mvmr_tauY,tauy_mvmr)
        mvmr_tauX=c(mvmr_tauX,taux_mvmr)
        mvmr_tau=c(mvmr_tau,tau_mvmr)
        ## two step MR
        ## first step to estimate tau X can include X
        onestep_taux=mr_ivw(betahat_x[indx_mvmr],betahat_m[indx_mvmr],se_x[indx_mvmr],se_m[indx_mvmr])
        onestep_tauy=mr_ivw(betahat_m[indm_single],betahat_y[indm_single],se_m[indm_single],se_y[indm_single])
        onestep_tau=onestep_taux$b*onestep_tauy$b
        twostep_tauY=c(twostep_tauY,onestep_tauy$b)
        twostep_tauX=c(twostep_tauX,onestep_taux$b)
        twostep_tau=c(twostep_tau,onestep_tau)
        twostep_tauY_sd=c(twostep_tauY_sd,onestep_tauy$se)
        twostep_tauX_sd=c(twostep_tauX_sd,onestep_taux$se)
        vector_C=c(onestep_tauy$b,onestep_taux$b)
        covariance_matrix=matrix(c(onestep_tauy$se^2,0,
                                   0,onestep_taux$se^2),nrow=2,ncol=2,byrow=TRUE)
        tau_sd=msm::deltamethod(~ x1*x2,vector_C,covariance_matrix)
        twostep_tau_sd=c(twostep_tau_sd,tau_sd)

        ##MVMR-IVW no measurement error correction
        BXGs=data.frame(cbind(betahat_x[ind_union],betahat_m[ind_union]))
        seBXGs=data.frame(cbind(se_x[ind_union],se_m[ind_union]))
        colnames(BXGs)=c("X","M")
        colnames(seBXGs)=c("X","M")
        input=format_mvmr(BXGs,betahat_y[ind_union],seBXGs,se_y[ind_union])
        ivw_mvmr_result=ivw_mvmr(input)
        ivw_tauX=mr_ivw(betahat_x[indx_mvmr],betahat_m[indx_mvmr],se_x[indx_mvmr],se_m[indx_mvmr])
        ivw_theta=ivw_mvmr_result[1,1]
        ivw_tauY=ivw_mvmr_result[2,1]
        MVMRivw_theta_se=ivw_mvmr_result[1,2]
        MVMRivw_tauY_se=ivw_mvmr_result[2,2]
        ivw_tau=ivw_tauX$b*ivw_tauY
        mvmr_theta_ivw=c(mvmr_theta_ivw,ivw_theta)
        mvmr_tauY_ivw=c(mvmr_tauY_ivw,ivw_tauY)
        mvmr_tauX_ivw=c(mvmr_tauX_ivw,ivw_tauX$b)
        mvmr_tau_ivw=c(mvmr_tau_ivw,ivw_tau)

        #thetalist=c()
        ivw_theta_se=c(ivw_theta_se,MVMRivw_theta_se)
        #tauYlist=c()
        ivw_tauY_se=c(ivw_tauY_se,MVMRivw_tauY_se)
        #tauXlist=c()
        ivw_tauX_se=c(ivw_tauX_se,ivw_tauX$se)
        ##difference of effects total effects of X on Y should select snps only associated with X theta can be estimated by mvmr
        total_effects=mr_ivw(betahat_x[indx_mvmr],betahat_y[indx_mvmr],se_x[indx_mvmr],se_y[indx_mvmr])
        differenc_tau=total_effects$b-ivw_theta
        total_effects_list=c(total_effects_list,total_effects$b)
        difference_tau_list=c(difference_tau_list,differenc_tau)
        total_effects_sd=c(total_effects_sd,total_effects$se)

        #rerand
        # ********Proposed MR***********
        # *****Rerandomization****
        indx=over2$filter1
        indm=over2$filter2
        overlap_index=intersect(indx,indm)#overlapping snps (for covariance estimation)
        overlap_m=which(indm%in%overlap_index) # the order in selected snps Sm
        overlap_x=which(indx%in%overlap_index) #new beta and se

        betanew_x=over2$gamma_exp1[indx]
        betanew_m=over2$gamma_exp2[indm]
        senew_x=over2$se1[indx]
        senew_m=over2$se2[indm]
        ele1=sum((betanew_x^2-senew_x^2)/(se_y[indx]^2))
        ele4=sum((betanew_m^2-senew_m^2)/(se_y[indm]^2))
        ele5=sum((betanew_x^2-senew_x^2)/(se_m[indx]^2))
        ele2=sum(over2$gamma_exp1[overlap_index]*over2$gamma_exp2[overlap_index]/(se_y[overlap_index]^2))+sum(over2$gamma_exp1[setdiff(indx,overlap_index)]*over2$gamma_exp2.carve[setdiff(indx,overlap_index)]/(se_y[setdiff(indx,overlap_index)]^2))
        ele3=sum(over2$gamma_exp1[overlap_index]*over2$gamma_exp2[overlap_index]/(se_y[overlap_index]^2))+sum(over2$gamma_exp1.carve[setdiff(indm,overlap_index)]*over2$gamma_exp2[setdiff(indm,overlap_index)]/(se_y[setdiff(indm,overlap_index)]^2))

        eleY1=sum(betanew_x*betahat_y[indx]/(se_y[indx]^2))
        eleY2=sum(betanew_m*betahat_y[indm]/(se_y[indm]^2))
        eleY3=sum(betanew_x[overlap_x]*over2$gamma_exp2[overlap_index]/(se_m[overlap_index]^2))+sum(over2$gamma_exp1[setdiff(indx,overlap_index)]*over2$gamma_exp2.carve[setdiff(indx,overlap_index)]/(se_m[setdiff(indx,overlap_index)]^2))
        #matrix
        ele2_original=c(ele2_original,sum(betanew_x[overlap_x]*betahat_m[overlap_index]/(se_y[overlap_index]^2)))
        ele2_rb=c(ele2_rb,ele2)
        sum_mat=matrix(c(ele1,ele2,0,ele3,ele4,0,0,0,ele5),nrow=3,ncol=3,byrow=TRUE)
        sum_vector=c(eleY1,eleY2,eleY3)
        estimator_rerand=solve(sum_mat)%*%sum_vector
        thetaRIVW=estimator_rerand[1]
        tauyRIVW=estimator_rerand[2]
        tauxRIVW=estimator_rerand[3]
        tauRIVW=tauyRIVW*tauxRIVW
        #covariance estimation
        regression_residuals1= sum( (betahat_y[overlap_index] *over2$gamma_exp1[overlap_index]-thetaRIVW * (over2$gamma_exp1[overlap_index]^2 - over2$se1[overlap_index]^2) -tauyRIVW*over2$gamma_exp1[overlap_index]*over2$gamma_exp2[overlap_index])^2/ se_y[overlap_index]^4)+sum( (betahat_y[setdiff(indx,overlap_index)] * over2$gamma_exp1[setdiff(indx,overlap_index)]-thetaRIVW * ( over2$gamma_exp1[setdiff(indx,overlap_index)]^2 -  over2$se1[setdiff(indx,overlap_index)]^2) -tauyRIVW*over2$gamma_exp1[setdiff(indx,overlap_index)]*over2$gamma_exp2.carve[setdiff(indx,overlap_index)])^2/ se_y[setdiff(indx,overlap_index)]^4)
        regression_residuals2= sum( (betahat_y[overlap_index] *over2$gamma_exp2[overlap_index]-tauyRIVW * (over2$gamma_exp2[overlap_index]^2 - over2$se2[overlap_index]^2) -thetaRIVW*over2$gamma_exp2[overlap_index]*over2$gamma_exp1[overlap_index])^2/ se_y[overlap_index]^4)+sum( (betahat_y[setdiff(indm,overlap_index)] * over2$gamma_exp2[setdiff(indm,overlap_index)]-tauyRIVW * ( over2$gamma_exp2[setdiff(indm,overlap_index)]^2 -  over2$se2[setdiff(indm,overlap_index)]^2) -thetaRIVW*over2$gamma_exp2[setdiff(indm,overlap_index)]*over2$gamma_exp1.carve[setdiff(indm,overlap_index)])^2/ se_y[setdiff(indm,overlap_index)]^4)
        regression_residuals3= sum( (betanew_m[overlap_m] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2))^2 / se_m[overlap_index]^4)+sum( (over2$gamma_exp2.carve[setdiff(indx,overlap_index)] * over2$gamma_exp1[setdiff(indx,overlap_index)] - tauxRIVW * (over2$gamma_exp1[setdiff(indx,overlap_index)]^2 - over2$se1[setdiff(indx,overlap_index)]^2))^2 / se_m[setdiff(indx,overlap_index)]^4)
        cov_12= sum( (betahat_y[overlap_index] * betanew_x[overlap_x]-thetaRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2) -tauyRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(betahat_y[overlap_index] * betanew_m[overlap_m] -tauyRIVW * (betanew_m[overlap_m]^2 - senew_m[overlap_m]^2) -thetaRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])/ se_y[overlap_index]^4)
        cov_13=sum( (betahat_y[overlap_index] * betanew_x[overlap_x]-thetaRIVW * (betanew_x[overlap_x]^2 - senew_x[overlap_x]^2) -tauyRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(over2$gamma_exp2[overlap_index] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2)) / (se_y[overlap_index]^2*se_m[overlap_index]^2))+sum( (betahat_y[setdiff(indx,overlap_index)] * over2$gamma_exp1[setdiff(indx,overlap_index)]-thetaRIVW * (over2$gamma_exp1[setdiff(indx,overlap_index)]^2 - over2$se1[setdiff(indx,overlap_index)]^2) -tauyRIVW*over2$gamma_exp1[setdiff(indx,overlap_index)]*over2$gamma_exp2.carve[setdiff(indx,overlap_index)])*(over2$gamma_exp2.carve[setdiff(indx,overlap_index)] * over2$gamma_exp1[setdiff(indx,overlap_index)] - tauxRIVW * (over2$gamma_exp1[setdiff(indx,overlap_index)]^2 - over2$gamma_exp1[setdiff(indx,overlap_index)]^2))/ (se_y[setdiff(indx,overlap_index)]^2*se_m[setdiff(indx,overlap_index)]^2))
       cov_23=sum((betahat_y[overlap_index] * betanew_m[overlap_m] -tauyRIVW * (betanew_m[overlap_m]^2 - senew_m[overlap_m]^2) -thetaRIVW*betanew_x[overlap_x]*betanew_m[overlap_m])*(over2$gamma_exp2[overlap_index] * betanew_x[overlap_x] - tauxRIVW * (betanew_x[overlap_x] ^2 - senew_x[overlap_x]^2))/(se_y[overlap_index]^2*se_m[overlap_index]^2))
        var_matrix=matrix(c(regression_residuals1,cov_12,cov_13,
                            cov_12,regression_residuals2,cov_23,
                            cov_13,cov_23,regression_residuals3),nrow=3,ncol=3,byrow=TRUE)
        var_matrix_final=solve(sum_mat)%*%var_matrix%*%t(solve(sum_mat))
        theta_rerandlist=c(theta_rerandlist,thetaRIVW)
        tauY_rerandlist=c(tauY_rerandlist,tauyRIVW)
        tauX_rerandlist=c(tauX_rerandlist,tauxRIVW)
        tau_rerandlist=c(tau_rerandlist,tauRIVW)

        v1_theta=var_matrix_final[1,1]
        v2_tauy=var_matrix_final[2,2]
        v3_taux=var_matrix_final[3,3]
        vector_C=c(tauyRIVW,tauxRIVW)
        sub_variance_matrix=var_matrix_final[2:3,2:3]
        sd_tau=msm::deltamethod(~ x1*x2,vector_C,sub_variance_matrix)
        theta_rerandse=c(theta_rerandse,sqrt(v1_theta))
        tauY_rerandse=c(tauY_rerandse,sqrt(v2_tauy))
        tauX_rerandse=c(tauX_rerandse,sqrt(v3_taux))
        tau_rerandse=c(tau_rerandse,sd_tau)
        #print(i)
      }

      result=data.frame(cbind(theta_rerandlist,theta_rerandse,tauX_rerandlist,tauX_rerandse,tauY_rerandlist,tauY_rerandse,tau_rerandlist,tau_rerandse,mvmr_theta,mvmr_tauX,mvmr_tauY,mvmr_tau,mvmr_theta_ivw,ivw_theta_se,mvmr_tauX_ivw,ivw_tauX_se,mvmr_tauY_ivw,ivw_tauY_se,mvmr_tau_ivw,twostep_tauX,twostep_tauX_sd,twostep_tauY,twostep_tauY_sd,twostep_tau,twostep_tau_sd,total_effects_list,total_effects_sd,difference_tau_list))
      colnames(result)=c("thetarerand","thetase","tauxrerand","tauxse","tauyrerand","tauyse","taurerand","tause","mvmr_theta","mvmr_tauX","mvmr_tauY","mvmr_tau","mvmrIVW_theta","mvmrIVW_theta_se","mvmrIVW_tauX","mvmrIVW_tauX_se","mvmrIVW_tauY","mvmrIVW_tauY_se","mvmrIVW_tau","2step_tauX","2step_tauX_sd","2step_tauY","2step_tauY_sd","2step_tau","2step_tau_sd","total_effect","total_effect_sd","difference_tau")
      estimation=data.frame(cbind(theta_rerandlist,tauX_rerandlist,tauY_rerandlist,tau_rerandlist,mvmr_theta,mvmr_tauX,mvmr_tauY,mvmr_tau,mvmr_theta_ivw,mvmr_tauX_ivw,mvmr_tauY_ivw,mvmr_tau_ivw,twostep_tauX,twostep_tauY,twostep_tau,total_effects_list,difference_tau_list))
      mean_val<-round(apply(result,2,mean),3)
      mc_sd<-round(apply(estimation,2,sd),3)
      cover_theta_rerand=cover(theta_rerandlist,theta_rerandse,theta)
      cover_tauX_rerand=cover(tauX_rerandlist,tauX_rerandse,tauX)
      cover_tauY_rerand=cover(tauY_rerandlist,tauY_rerandse,tauY)
      cover_tau_rerand=cover(tau_rerandlist,tau_rerandse,tau)
      cover_theta_ivw=cover(mvmr_theta_ivw,ivw_theta_se,theta)
      cover_tauX_ivw=cover(mvmr_tauX_ivw,ivw_tauX_se,tauX)
      cover_tauY_ivw=cover(mvmr_tauY_ivw,ivw_tauY_se,tauY)
      cover_tau_2step=cover(twostep_tau,twostep_tau_sd,tau)
      cover_tauX_2step=cover(twostep_tauX,twostep_tauX_sd,tauX)
      cover_tauY_2step=cover(twostep_tauY,twostep_tauY_sd,tauY)

      cover_total=cover(total_effects_list,total_effects_sd,tau+theta)

      abias_theta_rerand=abias(theta_rerandlist,theta)
      abias_tauX_rerand=abias(tauX_rerandlist,tauX)
      abias_tauY_rerand=abias(tauY_rerandlist,tauY)
      abias_tau_rerand=abias(tau_rerandlist,tau)
      abias_theta_ivw_mod=abias(mvmr_theta,theta)
      abias_tauX_ivw_mod=abias(mvmr_tauX,tauX)
      abias_tauY_ivw_mod=abias(mvmr_tauY,tauY)
      abias_tau_ivw_mod=abias(mvmr_tau,tau)
      abias_theta_ivw=abias(mvmr_theta_ivw,theta)
      abias_tauX_ivw=abias(mvmr_tauX_ivw,tauX)
      abias_tauY_ivw=abias(mvmr_tauY_ivw,tauY)
      abias_tau_ivw=abias(mvmr_tau_ivw,tau)


      abias_tauX_2step=abias(twostep_tauX,tauX)
      abias_tauY_2step=abias(twostep_tauY,tauY)
      abias_tau_2step=abias(twostep_tau,tau)
      abias_total=abias(total_effects_list,tau+theta)
      abias_tau_difference=abias(difference_tau_list,tau)




      cover_value=data.frame(cbind(cover_theta_rerand,cover_tauX_rerand,cover_tauY_rerand,cover_tau_rerand,cover_theta_ivw,cover_tauX_ivw,cover_tauY_ivw,cover_tauX_2step,cover_tauY_2step,cover_tau_2step,cover_total))
      abias_prop=data.frame(cbind(abias_theta_rerand/abs(theta),abias_tauX_rerand/abs(tauX),abias_tauY_rerand/abs(tauY),abias_tau_rerand/abs(tau),abias_theta_ivw/abs(theta),abias_tauX_ivw/abs(tauX),abias_tauY_ivw/abs(tauY),abias_tau_ivw/abs(tau),abias_theta_ivw_mod/abs(theta),abias_tauX_ivw_mod/abs(tauX),abias_tauY_ivw_mod/abs(tauY),abias_tau_ivw_mod/abs(tau),abias_tauX_2step/abs(tauX),abias_tauY_2step/abs(tauY),abias_tau_2step/abs(tau),abias_total/abs(theta+tau),abias_tau_difference/abs(tau)))
      colnames(cover_value)=c("theta_rerand","tauX_rerand","tauY_rerand","tau_rerand","theta_ivw","tauX_ivw","tauY_ivw","tauX_2step","tauY_2step","tau_2step","cover_total")
      colnames(abias_prop)=c("theta_rerand","tauX_rerand","tauY_rerand","tau_rerand","theta_ivw","tauX_ivw","tauY_ivw","tau_ivw","theta_ivw_mod","tauX_ivw_mod","tauY_ivw_mod","tau_ivw_mod","tauX_2step","tauY_2step","tau_2step","total effects","diff tau")
     #Please change the directory in your local computer
       save.image(paste0("D:\\Research\\structural equation\\ResultOCToriginal\\simudesign",design,"pi1=",pi1,"pi3=",pi3,"0.5e10-41e-4theta=",theta,"tauY=",tauY,"tauX=",tauX,"M=",M,".Rdata"))

    }
  }}}
}}


