#' Main function for MAGIC
#'
#' @param beta.exposure  SNP effect size's vector of the exposure vairable (GWASI)
#' @param beta.mediator  SNP effect size's vector of the mediator vairable (GWASIII)
#' @param beta.outcome   SNP effect size's vector of the outcome vairable (GWASII)
#' @param se.exposure    SNP effect size's standard errors of \code{beta.exposure}
#' @param se.mediator    SNP effect size's standard errors of \code{beta.mediator}
#' @param se.outcome     SNP effect size's standard errors of \code{beta.outcome}
#' @param Conflevel Confidence level. Default is 0.95.
#' @param lambda A vector of specified pre-screening threshold in the ordering of (exposure, mediator). Default is c(4.06,4.06). (corresponding p value is 5e-5)
#' @param eta A vector of rerandomized scale in the ordering of  (exposure, mediator). Default is c(0.5,0.5).
#' @return A list
#' \describe{
#' \item{theta.hat}{Estimated direct effect from exposure to outcome variable}
#' \item{tauy.hat}{Estimated direct effect from mediator to outcome variable}
#' \item{taux.hat}{Estimated indirect effect from exposure to mediator variable}
#' \item{tau.hat}{Estimated mediation}
#' \item{theta.se}{Standard error of \code{theta.hat}}
#' \item{tauy.se}{Standard error of \code{tauy.hat}}
#' \item{taux.se}{Standard error of \code{taux.hat}}
#' \item{tau.se}{Standard error of \code{tau.hat}}
#' \item{n.IV.exp}{Number of IVs used in exposure dataset}
#' \item{n.IV.med}{Number of IVs used in mediator dataset}
#' }
#'
#' @references Rita Qiuran Lyu, Chong Wu, Xinwei Ma, Jingshen Wang (2023). Mediation Analysis with Mendelian Randomization and Efficient Multiple GWAS Integration. \url{https://arxiv.org/abs/2312.10563}.
#'
#' @import stats
#' @export
#'
#' @examples
#'
