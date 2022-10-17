#Auxiliary functions


#Interaction structures
f.LLRep<-function(y){
  A <- nrow(y); J <- ncol(y);
  Z <- ifelse(y>0, log(y), 0)
  alpha0 <- mean(Z)
  alphaa <- rowSums(Z)/J-alpha0
  alphaj <- colSums(Z)/A-alpha0
  alphaaj <- ifelse(y>0, Z- matrix(alphaa,nrow=A,ncol=J,byrow=F)-
                      matrix(alphaj,nrow=A,ncol=J,byrow=T)-
                      matrix(alpha0,nrow=A, ncol=J), 0)
  alphaaj
  list(alpha_0 = alpha0, alpha_a = alphaa, alpha_j = alphaj, alpha_aj = alphaaj)
}


#Multiplication of vectors with NA
multiplic.na <- function(matriz1, matriz2) {
  rows <- nrow(matriz1)
  cols <- ncol(matriz2)
  matriz <- matrix(nrow = rows, ncol = cols)
  for (i in 1:rows){
    for(j in 1:cols){
      if(is.na(i)){
        i<- i+1
      } else{
        vec1 <- matriz1[i,]
        vec2 <- matriz2[,j]
        mult <- vec1 * vec2
        matriz[i,j] <- sum(mult, na.rm=TRUE)
      }
    }
  }
  return(matriz)
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
