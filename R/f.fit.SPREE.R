#' spree: Structure Preserving Estimation
#'
#' Function \code{spree} updates counts and proportions arranged by domains
#' (areas) and a categorical variable (at least 2 categories). There are three
#' alternatives based on \cite{Purcell, N. J. and L. Kish (1980)},
#' \cite{Zhang, L.-C. and R. L. Chambers (2004)}, and
#' \cite{Luna-Hernández, A. (2016)}
#'
#' @param proxy a data frame that with exactly same size as \code{sample_data}
#' (A rows and J columns). It contains also \code{population_domains},
#' \code{row_margins}, and \code{col_margins}. No NAs are allowed. If this file
#'  contains zeros, they are treated as structural zeros, so they will remain as
#'   zeros in the output.
#'
#' @param survey This is the AxJ table obtained from the survey data. It is
#' **not** required if estimator is "spree". If method is "ml", survey contains
#' the observed survey counts in each category. If method is "iwls", survey
#' contains the estimated population counts using the expansion factors. Please
#' notice that when iwls is used, the survey table is internally transformed to
#'  the within-area proportions, so the expansion factors will only affect the
#'  estimate if the sampling design is *informative*, in the sense that the
#'  estimated within-area proportions are different from the observed
#'   within-area proportions in the sample.
#'
#' @param estimator Can take the values "spree" (no update of structure),
#' "gspree" (one coefficient) or "mspree" (J-1xJ-1 independent coefficients).
#' @param method Can take the values "ml" for Maximum Likelihood (Poisson or
#'  Multinomial are equivalent) or  "iwls" for Iterative Weighted Least Squares.
#'   ML estimates are obtained using the observed sample counts and assuming a
#'   Poisson/Multinomial distribution (both are equivalent). IWLS estimates are
#'    obtained using estimated counts obtained with sampling weights, and allow
#'    for the specification of category-specific design effects.
#' @param iwls.opt This is a list containing the additional elements required
#' for method "iwls". Has components nsa. (vector of realized area sample sizes),
#'  ini.prop (AxJ matrix of within-area proportions which will be used to
#'  initialize the algorithm), and deff (vector of size J containing
#'  column-specific design effects. If missing, it is set to 1).

#'
#' @references
#' Luna-Hernandez (2016) Multivariate Structure Preserving Estimation
#' for Population Compositions. Ph.D. thesis, University of Southampton.\cr \cr
#' Purcell, N.J. & Kish, L. (1980) Postcensal estimates for local areas
#' (or domains). International Statistical Review/Revue Internationale de
#' Statistique, 48, 3– 18.\cr \cr
#' Zhang, L.C. & Chambers, R.L. (2004) Small area estimates for
#' cross-classifications. Journal of the Royal Statistical Society.
#' Series B: Statistical Methodology, 66, 479– 496. \cr \cr
#'@examples
#' \dontrun{
# Loading proxz and sample data
#' data("census")
#' data("survey_17")
#'
#' Substituting the 0's in the proxy by 3.
#' censusa <- census
#' censusa[censusa==0] <- 3
#' census <- censusa
#' rm(censusa)
#' Increasing the size of the counts in the last two columns to avoid issues
#' with MSE estimation later on.
#' survey_12$Two <- survey_12$Two*10
#' survey_12$Three_Four <- survey_12$Three_Four*10
#' Filtering a common set of areas in proxy and survey and eliminating the area identifier:
#' old.proxy <- census[1:70,-1]
#' survey.data <- survey_12[1:70,-1]
#' Fixing an issue with the margins, as they don't sum exactly the same amount.
#' data("P2012c")
#' column_tot_new<-function(data_true, domains, total_true){
#' data_true<-as.data.frame(select(data_true, -one_of(domains)))
#' props<-prop.table(colSums(data_true, na.rm = T))
#' column_total <- NULL
#' for (i in 1:ncol(data_true)) {
#'    column_total <- rbind(column_total,as.vector(c((total_true)*props[i])))
#'  }
#'  return(column_total=c(column_total))
#' }
#' row.mar <-P2012c[1:70]
#' col.mar <- column_tot_new(data_true= survey_12, domains = "Canton", total_true= sum(P2012c[1:70]))
#' ctotals_12 <- column_tot_new(data_true= survey_12, domains = "Canton", total_true= sum(P2012c))
#'
#' SPREE estimates
#' SPREE does not require the estimation of a coefficient. However, to make the
#' treatment analogous to that of GSPREE and MSPREE, we also require a f.fit.SPREE step.
#' method has no effect in this case.
#' res1 <- f.fit.SPREE(proxy=old.proxy, survey=survey.data,
#'                    estimator="spree", method="ml")
#' res1
#'
#' SPREE estimates for the same proxy as in res. Only row margins are used. }
#' }


f.fit.SPREE <- function(proxy, survey=NULL,
                        estimator=c("spree","gspree", "mspree"),
                            method=c("ml", "iwls"), iwls.opt = list(nsa.=NULL, ini.prop =NULL, deff=NULL))
#for ml, survey contains the observed sample counts,
#for iwls, survey contains the estimated counts (using survey weights)
{
  Xaj <- as.matrix(proxy)
  A <- nrow(Xaj); J <- ncol(Xaj); AJ <- A*J
  alpha.Xaj.X <- as.matrix(f.LLRep(Xaj,A,J)$alpha_aj)
  if (estimator == "spree"){
    coeff <- list(beta=1)
    z <- NULL
    var.coeff <- 0}
  else {
    if (!is.null(survey))
      Yaj <- as.matrix(survey)
    col <- kronecker(matrix(seq(1:J),1,J), matrix(1,1,A))
    row <- kronecker(matrix(1,1,J), matrix(seq(1:A),1,A))
    zg <- kronecker(matrix(1,J,1), diag(A))
    matrixT <-rbind(diag(J-1), rep(-1,J-1))
    zl <- kronecker(matrixT, matrix(1,A,1))
    if (estimator == "gspree")
      z <- cbind(zg, zl, c(alpha.Xaj.X))
    if (estimator == "mspree")
      z <- cbind(zg, zl, kronecker(matrixT,(alpha.Xaj.X %*% matrixT)))
    if (method == "ml"){
      model<- stats::glm(c(Yaj) ~ -1 + z, family = "poisson")
      model$varcov <- vcov(model)
    }
    if(method=="iwls")
        model<- f.fit.iwls(paj=Yaj/rowSums(Yaj), design.matrix=z, nsa.=iwls.opt$nsa.,
                           ini.prop=iwls.opt$ini.prop, deff=iwls.opt$deff )
    if (estimator == "gspree"){
        coeff <- list(beta=matrix(model$coefficients[(A+J)]))
        var.coeff <- model$varcov[(A+J),(A+J)]
    }
    if (estimator == "mspree"){
        betat <- matrix(model$coefficients[-c(1:(A+J-1))],J-1,J-1,byrow = F)
        beta <- cbind(betat, -apply(betat, 1, sum))
        beta <- rbind(beta, -apply(beta, 2, sum))
        Bp <- matrix(NA, J, J)
        for (k in 1:J) {
          for (l in 1:J)
            Bp[k,l] <- beta[k,l] + (1/(J-2))*(beta[k,k] +
                          beta[l, l]-(sum(beta[col(beta) == row(beta)]))/(J-1))
        }
        coeff <- list (beta=beta, Bp=Bp)
        var.coeff <- diag(model$varcov[-c(1:(A+J-1)),-c(1:(A+J-1))])
    }
  }
  list(coeff=coeff, var.coeff = var.coeff,proxy= Xaj, estimator=estimator, proxy.alpha.aj=alpha.Xaj.X)
}
