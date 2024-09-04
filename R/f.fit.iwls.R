###############################
#IWLS fitting

f.fit.iwls <- function(paj, design.matrix, nsa.=NULL, ini.prop=NULL, deff=NULL, max.iter=500, eps=1e-5){
  #paj is the table of within-row proportions estimated with sampling weights
  A <- nrow(ini.prop)
  J <- ncol(ini.prop)
  if(is.null(ini.prop)) stop("IWLS is chosen but no initial table of within-area proportions is provided.\n 
                               Please run the SPREE first to obtain an initial estimate.\n" )
  if(is.null(nsa.)) stop("IWLS is chosen but area sample sizes are not provided.\n")
  if (is.null(deff)||(min(deff)<=0)) {
    warning("IWLS is chosen but design effect is missing or takes an invalid value. It will be set as 1.")
    deff <- c(rep(1, J))
  }
  yaj <- as.matrix(c(as.matrix(paj*nsa.)))
  mu.hat = ini.prop*nsa.
  
  conv <- 0
  iter <- 1
  beta.ini <- matrix(NA,max.iter,ncol(design.matrix));
  #initial point
  yajstar <- round(yaj)
  model<- stats::glm(c(yajstar) ~ -1 + design.matrix, family = "poisson")
  beta.ini[1,] <- model$coefficients
  while (conv == 0 & iter < max.iter) {
    eta <- design.matrix%*%beta.ini[iter,]
    mu <- exp(eta)
    #W <- solve(diag(1/c(mu))%*%diag(c(mu))%*%diag(1/c(mu)))
    W <- diag(c(matrix(1/deff,A,J,byrow=T))*c(mu))
    XtWX <- solve(t(design.matrix)%*%W%*%design.matrix)
    eta.tilde <- (yaj-mu)/mu+eta
    #print(dim(eta.tilde))
    beta.ini[(iter+1),] <- XtWX%*%t(design.matrix)%*%W%*%eta.tilde
    if (max(abs(beta.ini[iter] - beta.ini[(iter+1)])) < eps) 
      conv <- 1
    iter <- iter+1
  }
  if (conv==1){
    list(coefficients = beta.ini[(iter),], varcov = XtWX)
  }
}