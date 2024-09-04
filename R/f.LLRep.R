f.LLRep <- function (data, A, J) 
{
  Z <- log(data)
  alpha0 <- sum(Z)/(A * J)
  alphaa <- apply(Z, 1, mean) - alpha0
  alphaj <- apply(Z, 2, mean) - alpha0
  alphaaj <- Z - matrix(alphaa, nrow = A, ncol = J, byrow = F) - 
    matrix(alphaj, nrow = A, ncol = J, byrow = T) - matrix(alpha0, 
                                                           nrow = A, ncol = J)
  list(alpha_0 = alpha0, alpha_a = alphaa, alpha_j = alphaj, 
       alpha_aj = alphaaj)
}