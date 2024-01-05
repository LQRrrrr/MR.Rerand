#' Main function for RIVW
#' @param beta.exposure  SNP effect size's vector of the exposure vairable (GWASI)
#' @param beta.outcome   SNP effect size's vector of the outcome vairable (GWASII)
#' @param se.exposure    SNP effect size's standard errors of \code{beta.exposure}
#' @param se.outcome     SNP effect size's standard errors of \code{beta.outcome}
#' @param Conflevel Confidence level. Default is 0.95.
#' @param smoothing  Whether to use smoothing to decrease variance Default is FALSE.
#' @param lambda the value of specified pre-screening threshold. Default is 4.06. (corresponding p value is 5e-5)
#' @param eta A vector of rerandomized scale. Default is 0.5.
#' @return A list
#' \describe{
#' \item{beta.hat}{Estimated direct effect from exposure to outcome variable}
#' \item{beta.se}{Standard error of \code{beta.hat}}
#' \item{n.IV}{Number of IVs used in exposure dataset}
#' \item{F}{The value of F-statistic}
#' }
#'
#' @references   Xinwei Ma, Jingshen Wang, Chong Wu. (2023). Breaking the Winner’s Curse in Mendelian Randomization:Rerandomized Inverse Variance Weighted Estimator \url{https://projecteuclid.org/journals/annals-of-statistics/volume-51/issue-1/Breaking-the-winners-curse-in-Mendelian-randomization--Rerandomized-inverse/10.1214/22-AOS2247.full}.
#'
#' @import stats
#'
#' @examples
#'
