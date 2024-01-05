## oracle selection
##simulation setting 1??not overlapping 2. partially overlapping 3.the same region
rm(list = ls())
theta=0.2
tauY=0.2
tauX=0.6
M=100000
pi_overlap=0.002
pix_minus_m=0.002
pim_other=0.001
nx=ny=nm=100000
se_x=rep(sqrt(1/nx),M)
se_y=rep(sqrt(1/ny),M)
se_m=rep(sqrt(1/nm),M)
sigma2x=sigma2y=0.5e-4
sigma2m=sigma2x_td=sigma2y_td=3e-4
p_tilde=c(0.1,0.5,1,3,6)


for (p in p_tilde)
{

      mvmr_theta=c()
      mvmr_tauY=c()
      theta_rerandlist=c()
      tauY_rerandlist=c()
      print(paste0("theta=",theta,"tauX=",tauX,"tauY=",tauY,"ptilde=",p))
      for (i in 1:1000)
      {
        gamma=rep(0,M)
        gammam=rep(0,M)
        indx_minus_m=sample(M,round(M*pix_minus_m))#valid for betaXj
        causalsnps_x= indx_minus_m
        causalsnps= indx_minus_m
        number=M*pix_minus_m/p

        indm_minus_x=sample(setdiff(1:M,causalsnps),round(number))
        causalsnps=c(causalsnps,indm_minus_x)
        causalsnps_m= indm_minus_x



        number=M*(pi_overlap)
        overlap_index=sample(setdiff(1:M,causalsnps),round(number))
        causalsnps_x=c(causalsnps_x,overlap_index)
        causalsnps=c(causalsnps,overlap_index)
        causalsnps_m=c(causalsnps_m,overlap_index)


        gamma[causalsnps_x]=rnorm(length(causalsnps_x),0,sd=sqrt(sigma2x))

        gammam=tauX*gamma
        gammam[causalsnps_m]=gammam[causalsnps_m]+rnorm(length(causalsnps_m),0,sd=sqrt(sigma2m))
        gammam[setdiff(1:M,causalsnps_m)]=0
        gammay=theta*gamma+tauY*gammam#no pleiotropy
        ##perfect selection
        ##
        betax=gamma
        betam=gammam
        betay=gammay

        betax=betax+rnorm(M,mean=0,sd=sqrt(1/nx))
        betam=betam++rnorm(M,mean=0,sd=sqrt(1/nm))
        betay=betay+rnorm(M,mean=0,sd=sqrt(1/ny))
        ele1=sum((betax[causalsnps_x]^2-se_x[causalsnps_x]^2)/(se_y[causalsnps_x]^2))
        ele2=sum(betax[causalsnps_x]*betam[causalsnps_x]/(se_y[causalsnps_x]^2))
        ele3=sum(betax[causalsnps_m]*betam[causalsnps_m]/(se_y[causalsnps_m]^2))
        ele4=sum((betam[causalsnps_m]^2-se_m[causalsnps_m]^2)/(se_y[causalsnps_m]^2))
        eleY1=sum((betay[causalsnps_x]*betax[causalsnps_x])/(se_y[causalsnps_x]^2))
        eleY2=sum((betay[causalsnps_m]*betam[causalsnps_m])/(se_y[causalsnps_m]^2))
        sum_mat=matrix(c(ele1,ele2,ele3,ele4),nrow=2,ncol=2,byrow=TRUE)
        sum_vector=c(eleY1,eleY2)
        estimator_rerand=solve(sum_mat)%*%sum_vector

        ele1_mvmr=sum((betax[causalsnps]^2-se_x[causalsnps]^2)/(se_y[causalsnps]^2))
        ele2_mvmr=sum(betax[causalsnps]*betam[causalsnps]/(se_y[causalsnps]^2))
        ele3_mvmr=sum(betax[causalsnps]*betam[causalsnps]/(se_y[causalsnps]^2))
        ele4_mvmr=sum((betam[causalsnps]^2-se_m[causalsnps]^2)/(se_y[causalsnps]^2))
        eleY1_mvmr=sum((betay[causalsnps]*betax[causalsnps])/(se_y[causalsnps]^2))
        eleY2_mvmr=sum((betay[causalsnps]*betam[causalsnps])/(se_y[causalsnps]^2))
        sum_mat_mvmr=matrix(c(ele1_mvmr,ele2_mvmr,ele3_mvmr,ele4_mvmr),nrow=2,ncol=2,byrow=TRUE)
        sum_vector_mvmr=c(eleY1_mvmr,eleY2_mvmr)
        estimator_mvmr=solve(sum_mat_mvmr)%*%sum_vector_mvmr
        mvmr_theta=c(mvmr_theta,estimator_mvmr[1])
        mvmr_tauY=c(mvmr_tauY,estimator_mvmr[2])
        theta_rerandlist=c(theta_rerandlist,estimator_rerand[1])
        tauY_rerandlist=c(tauY_rerandlist,estimator_rerand[2])
        #print(i)
      }
      result=data.frame(cbind(theta-theta_rerandlist,tauY-tauY_rerandlist,theta-mvmr_theta,tauY-mvmr_tauY))
      colnames(result)=c("thetarerand","tauyrerand","thetamvmr","tauymvmr")
      #asybias<-round(apply(result*sqrt(M),2,mean),4)
      mc_sd<-round(apply(result,2,sd),4)
      #asyvar_mc<-round(apply(result,2,var)*M,4)
      bias=round(abs(apply(result,2,mean)),4)
      estimation=data.frame(cbind(theta_rerandlist,tauY_rerandlist,mvmr_theta,mvmr_tauY))
      colnames(estimation)=c("thetarerand","tauyrerand","thetamvmr","tauymvmr")
      value=round(apply(estimation,2,mean),4)
      final_result=t(rbind(value,bias,mc_sd))
      #Please update to you own directory
      save.image(paste0("D:\\Research\\structural equation\\power_analysisNOV\\mvmr_effeciency_ptilde=",p,"pi_overlap=",pi_overlap,"sigma_x=",sigma2x,"sigma_m=",sigma2m,".Rdata"))
      print(final_result)
    }


