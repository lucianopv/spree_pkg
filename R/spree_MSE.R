
mse_spree<- function(YS,
                     Ya.,
                     Y.j,
                     popul_data,
                     sampl_data,
                     A,
                     J,
                     AJ,
                     type,
                     method,
                     design_effect,
                     B = 100){

  YS <- as.matrix(YS)
  MSE_b <- array(dim=c(AJ,B))
  updated_YS_true <- matrix(YS, AJ, byrow = F)
  MSE<-rep(0,AJ)

  for(b in 1:B){
    set.seed(b)

    #1. Generate pseudo-populations from Multinomial distibution with

    theta_p <- loglin(outer(Ya.,Y.j)/sum(Ya.),margin=list(1,2),
                      start=YS, fit=TRUE, eps=1.e-05, iter=100)$fit

    pi_p <- theta_p/rowSums(theta_p)

    na<- round((sum(popul_data)/sum(sampl_data))*Ya.,0)

    #population
    X_p <- NULL
    for (a in 1:A){
      set.seed(32*a+b*41)
      X_p <- cbind(X_p,(rmultinom(1, Ya.[a], pi_p[a,])))
    }
    X_p <- t(X_p)
    X_p <- loglin(outer(Ya.,Y.j)/sum(Ya.),margin=list(1,2),
                  start=X_p, fit=TRUE, eps=1.e-05, iter=100)$fit
    theta_s <- X_p/rowSums(X_p)

    #sample
    y_s <- NULL
    for (a in 1:A){
      set.seed(8*a+b*37)
      y_s <- cbind(y_s,(rmultinom(1, na[a], theta_s[a,])))
    }
    #here sample scenario
    y_s <- t(y_s)


    ####fitting


    #Obtener la matriz con las interacciones del censo

    alpha.Xp <- as.matrix(f.LLRep(X_p, A, J)$alpha_aj)


    # ID columns y rows
    col <- kronecker(matrix(seq(1:J),1,J),matrix(1,1,A)) #For table mspree
    row <- kronecker(matrix(1,1,J),matrix(seq(1:A),1,A)) #For table mspree


    # Design matrix

    # Row effect (domains)
    zg <- kronecker(matrix(1,J,1),diag(A)) # for table mspree

    # Column effect (categories)

    zl <- kronecker(rbind(diag(J-1),rep(-1,J-1)),matrix(1,A,1)) #mspree

    des <- alpha.Xp[,1:(J-1)]-alpha.Xp[,J] #mspree
    des[des==0]<-NA
    zalpha <- rbind(kronecker(diag(J-1),des),kronecker(matrix(1,1,J-1),-des)) ###aqui hay que hacer lo mismo que con parte1.1
    z <- cbind(zg,zl,zalpha)
    z[is.na(z)]<-0

    z.gspree <- cbind(zg,zl,c(alpha.Xp))  #gspree


    #spree

    if(type=="SPREE"){

    YS_mse <- loglin(outer(Ya., Y.j)/sum(Ya.), margin= list(1,2), start=X_p,
                     fit = TRUE, eps = 1.e-05, iter = 100)$fit


    MSE_b[,b]<-c(YS_mse)

    } else if(type== "GSPREE") {


    #gspree
    model_g <- glm(c(y_s) ~ -1+z.gspree,family = "poisson")

    ayfit.g<- model_g$coefficients[A+J]*alpha.Xp

    YS_mse <- loglin(outer(Ya., Y.j)/sum(Ya.),margin=list(1,2),
                     start=exp(ayfit.g), fit=TRUE, eps=1.e-05, iter=100)$fit

    MSE_b[,b]<-c(YS_mse)

    } else if(type== "MSPREE") {
    #mspree

    if( method == "ML") {

      if(ncol(y_s)<3) {
        warning("GSPREE is applied because there are only two categories.")
      }

      model_m <- glm(c(y_s) ~ -1+z,family = "poisson")
      YS_mse  <- loglin(outer(Ya.,Y.j)/sum(Ya.),margin=list(1,2),
                       start=model_m$fitted, fit=TRUE, eps=1.e-05, iter=100)$fit

      MSE_b[,b]<-c(YS_mse)

    } else if(method == "IWLS") {

        if(is.null(design_effect)){
          design_effect<-c(rep(1,J))
          warning("Design effect is set as 1.")
        }

        mu_hat_p<- loglin(outer(Ya.,Y.j)/sum(Ya.),margin=list(1,2),
                          start=X_p, fit=TRUE, eps=1.e-05, iter=100)$fit
        theta_hat_p <- mu_hat_p/rowSums(mu_hat_p)


        matrixT_p <- rbind(diag(J-1),matrix(-1,1,J-1))
        #for the multinomial, the last column is the reference
        talpha_p <- alpha.Xp%*%matrixT_p  #auxiliary structure
        colJ <- kronecker(matrix(seq(1:J),1,J),matrix(1,1,A))
        x <- cbind(zl[colJ!=J,],kronecker(diag(J-1),talpha_p))


        theta_dir_p <- y_s/rowSums(y_s)
        theta_dir_p <- theta_dir_p[,-J]
        theta_p <-theta_hat_p[,-J]

        deffs<-design_effect[-J]

        Vars_p<-vcovs.srs.prob(theta_p, Ya., deffs)

        col <- kronecker(matrix(seq(1:(J-1)),1,J-1),matrix(1,1,A))
        row <- kronecker(matrix(1,1,(J-1)),matrix(seq(1:A),1,A))
        conv <- 0
        iter <- 1
        beta_o_p <- rep(0,ncol(x))
        #if (class(vcovs)=="matrix")
        V <- Vars_p

        while (conv ==0 & iter < 500  ){

          theta_ref_p <- 1-rowSums(theta_p)
          grad_p <- matrix(0,A*(J-1),A*(J-1))
          for (r in 1:A){
            temp_p <- matrix(1/theta_ref_p[r],J-1,J-1)+diag(1/theta_p[r,])
            grad_p[row==r,row==r] <- temp_p
          }
          z_iwls_p <- c(log(theta_p/theta_ref_p))+
            (grad_p%*%c(theta_dir_p-theta_p))

          # if (class(vcovs)=="function")
          #  V <- vcovs(theta, na)
          W_p <- solve(grad_p%*%V%*%grad_p)

          H <- solve(t(x)%*%W_p%*%x)
          beta_n_p <- H%*%t(x)%*%W_p%*%z_iwls_p
          exbeta_p <-matrix(exp(x%*%beta_n_p),A,J-1,byrow=F)
          theta_p <- exbeta_p*matrix(1/(1+rowSums(exbeta_p)),A,J-1,byrow=F)
          if (max(abs(beta_n_p-beta_o_p))<0.00001) conv <- 1 else {
            beta_o_p <- beta_n_p
            iter <- iter+1
          }
        }
        beta_logit_p <- beta_n_p[-(1:(J-1))]
        transf_back_p <- (diag((J-1)*(J-1))-(1/J)*kronecker(matrix(1,J-1,J-1),diag(J-1)))
        beta_MSPREE_p <- transf_back_p%*%beta_logit_p

        ayfit_IWLS_p<- z[, -c(1:(A+J-1))]%*%beta_MSPREE_p
        ayfit_IWLS_p<- matrix(ayfit_IWLS_p, A, J, byrow = F)
        YS_mse<-loglin(outer(Ya., Y.j)/sum(Ya.), margin=list(1,2),
                       start=exp(ayfit_IWLS_p), fit = TRUE, eps = 1.e-05, iter = 100)$fit


        MSE_b[,b]<-c(YS_mse)

      }
    }
  }
        MSE<-rep(0,AJ)
        for (i in 1:AJ){

          MSE[i]<-mean((X_p[i]-MSE_b[i,])^2) #uncond
         # MSE[i]<-mean((YS[i]-MSE_b[i,])^2) #fp (mal!)
        }
        MSE<- matrix(MSE, A, J, byrow = F)

        RMSE<-sqrt(MSE)

        CV<- (sqrt(MSE)/YS)*100

        list(MSE   = MSE,
             RMSE = RMSE,
             CV    = CV)
      }



