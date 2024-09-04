
f.predict.SPREE <- function(fitted, new.proxy=NULL, row.margin, col.margin=NULL)
{
  if (!is.null(new.proxy)) new.proxy <- as.matrix(new.proxy)
  if (!is.null(col.margin)){
    mar <- list(1,2)
    indep <- outer(row.margin, col.margin)/sum(row.margin)
  }
  else{
    mar <- list(1)
    if (!is.null(new.proxy)) indep <- new.proxy*(row.margin/rowSums(new.proxy))
   else 
       indep <- fitted$proxy*(row.margin/rowSums(fitted$proxy))
  } 
  if (fitted$estimator == "spree"){
    if (!is.null(new.proxy))
      YS <- loglin(table=indep, margin=mar, start=new.proxy, fit = TRUE, eps = 1e-05, iter = 100)$fit
    else YS <- loglin(table=indep, margin=mar, start=fitted$proxy, fit = TRUE, eps = 1e-05, iter = 100)$fit
  }  else{
    if (!is.null(new.proxy))
      alpha.Xaj <- as.matrix(f.LLRep(as.matrix(new.proxy))$alpha_aj)
    else  alpha.Xaj <- fitted$proxy.alpha.aj     
  if(fitted$estimator == "gspree") beta <- diag(rep(fitted$coeff$beta,ncol(fitted$proxy.alpha.aj)))
  if(fitted$estimator=="mspree") beta <- fitted$coeff$beta
  start.table <- exp(alpha.Xaj%*%beta)
  YS <- loglin(table=indep, margin=mar, start=start.table, fit = TRUE, eps = 1e-05, iter = 100)$fit 
  }
  YP <- YS/apply(YS, 1, sum)
  list(counts = YS, proportions = YP)
}
