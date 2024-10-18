#' spree: Structure Preserving Estimation
#'
#' Function \code{spree} updates counts and proportions arranged by domains
#' (areas) and a categorical variable (at least 2 categories). There are three
#' alternatives based on \cite{Purcell, N. J. and L. Kish (1980)},
#' \cite{Zhang, L.-C. and R. L. Chambers (2004)}, and
#' \cite{Luna-Hernández, A. (2016)}
#'
#' @param population_data a data frame that with exactly same size as
#' \code{sample_data} (A rows and J columns). It contains also
#' \code{population_domains},\code{row_margins}, and \code{col_margins}.
#' @param sample_data a data frame that with exactly same size as
#' \code{population_data}. (A rows and J columns). It contains also
#' \code{population_domains},\code{row_margins}, and \code{col_margins}.
#' @param population_domains a character string containing the name of a
#' variable that indicates domains in the sample data. The variable can be
#' numeric or a factor but needs to be of the same class as the variable named
#' in \code{sample_domains}.
#' @param sample_domains a character string containing the name of a
#' variable that indicates domains in the population data. The variable can be
#' numeric or a factor but needs to be of the same class as the variable named
#' in \code{population_domains}.
#' @param row_margins a numeric vector with A elements containing the ‘true’
#' domain sizes in t1. NULL is set as default and the
#' sum by rows in sample_data is used.
#' @param col_margins numeric vector with J elements containing the ‘true’
#' totals for each category. The sum by columns in sample_data is used as default.
#' @param type a character string. SPREE version to implement: (i) ‘SPREE’;
#' (ii) ‘GSPREE’; (iii) ‘MSPREE’. If this argument remains empty, the first
#' option is selected. The option ‘MSPREE’ requires J >= 3. If J = 2, the
#' ‘GSPREE’ will be applied.
#' @param B number of bootstraps used in the MSE estimation.
#' @param method a character string. Specify the method to obtain estimates of
#' the beta coefficient. This  argument is only required when \code{type} =
#' ‘MSPREE’. The alternatives are (i) ‘ML’ and (ii) ‘IWLS’.
#' @param design_effect a character string. If type =‘MSPREE’ and method= ‘IWLS’
#'  are selected, a design effect for the column totals can be provided. A
#'  vector containing 1’s is used as default.
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
#' # Load census and sample data
#' data("census")
#' data("survey_17")
#'
#' survey_17f <- prep_sample_data(population_domains = "Canton",
#'                               sample_domains = "Canton", population_data =
#'                                 census, sample_data = survey_17)
#'
#' ctotals_17 <- column_tot(data_true= survey_17f, domains = "Canton", total_true= sum(P2017c))
#'
#' SPREE_17 <-spree(population_domains = "Canton",
#'                        sample_domains = "Canton", population_data = census,
#'                        sample_data = survey_17f, row_margins = P2017c,
#'                        col_margins = ctotals_17, type = "SPREE")
#'
#' MSPREE_IWLS_17 <-spree(population_domains = "Canton",
#'                        sample_domains = "Canton", population_data = census,
#'                        sample_data = survey_17f, row_margins = P2017c,
#'                        col_margins = ctotals_17, type = "MSPREE",
#'                        method = "IWLS", design_effect = c(1,1,1,1))
#' }
#' @importFrom stats loglin glm dplyr
#' @importFrom dplyr n_distict




spree<- function(population_domains,
                 sample_domains,
                 population_data,
                 sample_data,
                 row_margins = NULL,
                 col_margins = NULL,
                 type = c("SPREE", "GSPREE","MSPREE"),
                 B = 100,
                 method = c("ML", "IWLS"),
                 design_effect = NULL){


  spree_check1(population_data=population_data, sample_data=sample_data)

  spree_check2(population_domains= population_domains, sample_domains=sample_domains,population_data=population_data, sample_data=sample_data)


  domains_vec <- as.data.frame(population_data[, population_domains])
  names(domains_vec)<- "domain"

  domains_vec2 <- as.data.frame(sample_data[, sample_domains])
  names(domains_vec2)<- "domain"

  if (!identical(as.character(domains_vec2$domain), as.character(domains_vec$domain))){
    stop('Domains must be the same in population and sample data.')
  }

  popul_data<- population_data[,-1]
  sampl_data<- sample_data[,-1]

  #Number of rows (domains)
  A <- nrow(popul_data)

  #Number of columns (categories)
  J <- ncol(popul_data)

  #Dimension
  AJ <- A*J


  #Preparing census
  Xaj<-as.matrix(popul_data)
  # Xaj[Xaj==0]<-1   # avoiding zeros in census
  Xaj[is.na(Xaj)]<-0  #avoiding NAs in census

  #Preparing survey
  Yaj<- as.matrix(sampl_data)
  Yaj[is.na(Yaj)]<-1
  Yaj[Yaj==0]<-1

  #Totals

  if(is.null(row_margins)) {
    Ya. <- rowSums(Yaj)
  } else {
    Ya.<-row_margins
  }

  if(is.null(col_margins)) {
    Y.j <- colSums(Yaj)
  } else {
    Y.j<-as.vector(col_margins)
  }

  if(!is.null(col_margins) & is.null(row_margins)){
    stop('row_margins argument is missing.')
  }

  if(is.null(col_margins) & !is.null(row_margins)){
    stop('col_margins argument is missing.')
  }





  #Get population matrix

  alpha.Xaj.X <- as.matrix(f.LLRep(Xaj, A, J)$alpha_aj)


  # ID columns y rows
  col <- kronecker(matrix(seq(1:J),1,J),matrix(1,1,A)) #For table mspree
  row <- kronecker(matrix(1,1,J),matrix(seq(1:A),1,A)) #For table mspree


  # Design matrix

  # Row effect (domains)
  zg <- kronecker(matrix(1,J,1),diag(A)) # for table mspree

  # Column effect (categories)

  zl <- kronecker(rbind(diag(J-1),rep(-1,J-1)),matrix(1,A,1)) #mspree

  des <- alpha.Xaj.X[,1:(J-1)]-alpha.Xaj.X[,J] #mspree
  des[des==0]<-NA
  zalpha <- rbind(kronecker(diag(J-1),des),kronecker(matrix(1,1,J-1),-des)) ###aqui hay que hacer lo mismo que con parte1.1
  z <- cbind(zg,zl,zalpha)
  z[is.na(z)]<-0

  z.gspree <- cbind(zg,zl,c(alpha.Xaj.X))  #gspree



  if(type=="SPREE"){

    #Procedimientos equivalentes
    # YS <- loglin(outer(Ya., Y.j)/sum(Ya.), margin= list(1,2),
    #              start=exp(alpha.Xaj.X), fit = TRUE, eps = 1.e-05, iter = 100)$fit
    YS <- loglin(outer(Ya., Y.j)/sum(Ya.), margin= list(1,2), start=(Xaj),
                 fit = TRUE, eps = 1.e-05, iter = 100)$fit

    point<- as.data.frame(cbind(domains_vec,YS))
    colnames(point)<-colnames(population_data)
    Beta<-1

    mse_spree_estimates<- mse_spree(YS=YS,
                                    Ya.=Ya.,
                                    Y.j=Y.j,
                                    popul_data=popul_data,
                                    sampl_data=sampl_data,
                                    A = A,
                                    J = J,
                                    AJ =AJ,
                                    type,
                                    method,
                                    design_effect,
                                    B = 100)

    MSE      = as.data.frame(cbind(domains_vec,mse_spree_estimates$MSE))
    CV      = as.data.frame(cbind(domains_vec,mse_spree_estimates$CV))
    colnames(MSE) <-colnames(CV)<-colnames(population_data)

    return(list(updated_point   =point,
                MSE       = MSE,
                CV        = CV,
                Beta =Beta))


  }else if(type== "GSPREE"){
    model.g <- glm(c(Yaj) ~ -1+z.gspree,family = "poisson")

    ayfit.g<- model.g$coefficients[A+J]*alpha.Xaj.X
    beta_gspree<-matrix(model.g$coefficients[-c(1:(A+J-1))])

    YS <- loglin(outer(Ya., Y.j)/sum(Ya.),margin=list(1,2),
                 start=exp(ayfit.g), fit=TRUE, eps=1.e-05, iter=100)$fit

    point<- as.data.frame(cbind(domains_vec,YS))
    colnames(point)<-colnames(population_data)

    mse_spree_estimates<- mse_spree(YS=YS,
                                    Ya.=Ya.,
                                    Y.j=Y.j,
                                    popul_data=popul_data,
                                    sampl_data=sampl_data,
                                    A = A,
                                    J = J,
                                    AJ =AJ,
                                    type,
                                    method,
                                    design_effect,
                                    B = 100)



    MSE      = as.data.frame(cbind(domains_vec,mse_spree_estimates$MSE))
    CV      = as.data.frame(cbind(domains_vec,mse_spree_estimates$CV))
    colnames(MSE) <-colnames(CV)<-colnames(population_data)

    return(list(updated_point   =point,
                MSE       = MSE,
                CV        = CV,
                Beta =beta_gspree))


  } else if(type=="MSPREE"){

    if( method == "ML") {

      if(ncol(Yaj)<3){
        warning("GSPREE is applied because there are only two categories.")
      }

      model.m <- glm(c(Yaj) ~ -1+z,family = "poisson")

      beta <- matrix(model.m$coefficients[-c(1:(A+J-1))],J-1,J-1,byrow=F)
      beta1 <- cbind(beta, -apply(beta,1,sum))
      beta1 <- rbind(beta1, -apply(beta1,2,sum))

      Bp <- matrix(NA,J,J)
      for (k in 1:J){
        for (l in 1:J)
          Bp[k,l] <- beta1[k,l]+(1/(J-2))*(beta1[k,k]+beta1[l,l]-(sum(beta1[col(beta1)==row(beta1)]))/(J-1))
      }
      beta.X <- beta1
      bp.X <- Bp

      #  parte1<- z[,-c(1:(A+J-1))]
      # parte1<- as.matrix(parte1)


      #ayfit <- multiplic.na(parte1,as.matrix(model.m$coefficients[-c(1:(A+J-1))]))

      # ayfit <- matrix(ayfit,A,J,byrow=F)

      # YS <- loglin(outer(Ya., Y.j)/sum(Ya.),margin=list(1,2),
      #               start=exp(ayfit), fit=TRUE, eps=1.e-05, iter=100)$fit
      YS <- loglin(outer(Ya., Y.j)/sum(Ya.),margin=list(1,2),
                   start=model.m$fitted.values, fit=TRUE, eps=1.e-05, iter=100)$fit
      point<- as.data.frame(cbind(domains_vec,YS))
      colnames(point)<-colnames(population_data)

      mse_spree_estimates<- mse_spree(YS=YS,
                                      Ya.=Ya.,
                                      Y.j=Y.j,
                                      popul_data=popul_data,
                                      sampl_data=sampl_data,
                                      A = A,
                                      J = J,
                                      AJ =AJ,
                                      type,
                                      method,
                                      design_effect,
                                      B = 100)

      MSE      = as.data.frame(cbind(domains_vec,mse_spree_estimates$MSE))
      CV      = as.data.frame(cbind(domains_vec,mse_spree_estimates$CV))
      colnames(MSE) <-colnames(CV)<-colnames(population_data)

      return(list(updated_point   =point,
                  MSE       = MSE,
                  CV        = CV,
                  Beta =Bp))

      # return(YS)

    }else if(method == "IWLS"){

      Xaj[Xaj==0]<-1
      Yaj[Yaj==0]<-1
      if(is.null(design_effect)){
        design_effect<-c(rep(1,J))
        warning("Design effect is set as 1.")
      }



      #variance covariance for IWLS
      vcovs.srs.prob <- function(theta.hat, na,deffs){
        A1<-nrow(theta.hat)
        J1 <- ncol(theta.hat)
        row1 <- kronecker(matrix(1,1,J1),matrix(seq(1:A1),1,A1))
        dtheta <- theta.hat*matrix(deffs,A1,J1,byrow=T)
        V <- matrix(0,A1*J1,A1*J1)
        for (q in 1:A1){
          vr <- (1/na[q])*(-as.matrix(dtheta[q,])%*%(dtheta[q,])+diag(dtheta[q,]*deffs))
          V[row1==q,row1==q] <- vr
        }
        V
      }
      mu_hat<- loglin(outer(Ya.,Y.j)/sum(Ya.),margin=list(1,2),
                      start=Xaj, fit=TRUE, eps=1.e-05, iter=100)$fit
      theta_hat <- mu_hat/rowSums(mu_hat)


      matrixT <- rbind(diag(J-1),matrix(-1,1,J-1))
      #for the multinomial, the last column is the reference
      talpha <- alpha.Xaj.X%*%matrixT  #auxiliary structure
      colJ <- kronecker(matrix(seq(1:J),1,J),matrix(1,1,A))
      x <- cbind(zl[colJ!=J,],kronecker(diag(J-1),talpha))


      theta_dir <- Yaj/rowSums(Yaj)
      theta_dir <- theta_dir[,-J]
      theta <-theta_hat[,-J]

      deffs<-design_effect[-J]

      Vars<-vcovs.srs.prob(theta, Ya., deffs)

      col <- kronecker(matrix(seq(1:(J-1)),1,J-1),matrix(1,1,A))
      row <- kronecker(matrix(1,1,(J-1)),matrix(seq(1:A),1,A))
      conv <- 0
      iter <- 1
      beta_o <- rep(0,ncol(x))


      #### modif
      #if (class(vcovs.srs.prob)=="matrix")
      V <- Vars

      while (conv ==0 & iter < 500  ){

        theta_ref <- 1-rowSums(theta)
        grad <- matrix(0,A*(J-1),A*(J-1))
        for (r in 1:A){
          temp <- matrix(1/theta_ref[r],J-1,J-1)+diag(1/theta[r,])
          grad[row==r,row==r] <- temp
        }
        z_iwls <- c(log(theta/theta_ref))+
          (grad%*%c(theta_dir-theta))

        #### start modif
        #if (class(vcovs.srs.prob)=="function")
        # V <- vcovs.srs.prob(theta, na)

        #### end modif
        W <- solve(grad%*%V%*%grad)

        H <- solve(t(x)%*%W%*%x)
        beta_n <- H%*%t(x)%*%W%*%z_iwls
        exbeta <-matrix(exp(x%*%beta_n),A,J-1,byrow=F)
        theta <- exbeta*matrix(1/(1+rowSums(exbeta)),A,J-1,byrow=F)
        if (max(abs(beta_n-beta_o))<0.00001) conv <- 1 else {
          beta_o <- beta_n
          iter <- iter+1
        }
      }
      beta_logit <- beta_n[-(1:(J-1))]
      transf_back <- (diag((J-1)*(J-1))-(1/J)*kronecker(matrix(1,J-1,J-1),diag(J-1)))
      beta_MSPREE <- transf_back%*%beta_logit
      beta_m <- c(beta_MSPREE)
      beta_m2 <- matrix(beta_m,J-1,J-1,byrow=F)
      beta_m3 <- cbind(beta_m2, -apply(beta_m2,1,sum))
      beta_m4 <- rbind(beta_m3, -apply(beta_m3,2,sum))

      Bp <- matrix(NA,J,J)
      for (k in 1:J){
        for (l in 1:J)
          Bp[k,l] <- beta_m4[k,l]+(1/(J-2))*(beta_m4[k,k]+beta_m4[l,l]-(sum(beta_m4[col(beta_m4)==row(beta_m4)]))/(J-1))
      }
      beta.X <- beta_m4
      bp.X <- Bp

      ayfit_IWLS<- z[, -c(1:(A+J-1))]%*%beta_MSPREE
      ayfit_IWLS<- matrix(ayfit_IWLS, A, J, byrow = F)
      YS<-loglin(outer(Ya., Y.j)/sum(Ya.), margin=list(1,2),
                 start=exp(ayfit_IWLS), fit = TRUE, eps = 1.e-05, iter = 100)$fit

      #

      #

      point<- as.data.frame(cbind(domains_vec,YS))
      colnames(point)<-colnames(population_data)
      # beta_mspree<-as.data.frame(beta_m4)

      mse_spree_estimates<- mse_spree(YS=YS,
                                      Ya.=Ya.,
                                      Y.j=Y.j,
                                      popul_data=popul_data,
                                      sampl_data=sampl_data,
                                      A = A,
                                      J = J,
                                      AJ =AJ,
                                      type,
                                      method,
                                      design_effect,
                                      B = 100)


      MSE      = as.data.frame(cbind(domains_vec,mse_spree_estimates$MSE))
      CV      = as.data.frame(cbind(domains_vec,mse_spree_estimates$CV))
      colnames(MSE) <-colnames(CV)<-colnames(population_data)

      return(list(updated_point   =point,
                  MSE       = MSE,
                  CV        = CV,
                  Beta = Bp))
    }
  }

}


