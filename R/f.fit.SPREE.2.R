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
