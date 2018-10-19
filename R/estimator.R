


#' Kernel estimation of a varying coefficients model of a model
#'
#' $$ Y_{ig} = W_{ig}'theta(C_g) + U_{ig}$$
#'
#' where $U_{ig} = W_{ig}'(Theta_{ig} - E[Theta_{ig}|C_g])
#'
#' At present the function implements the uniform kernel. This functions estimates
#' $$ theta(cbar) = E[theta(C_g) | C_g = cbar]$$.
#'
#' The function multiplies each row of the data times a kernel weight. With these
#' transformed variables we compute an OLS regression.
#'
#' @param cbar A specific (non-random) point at which to evaluate $theta(r)$.
#' @param h The selected bandwidth.
#' @param Y The outcome variable.
#' @param W The regressor variable.
#' @param C The variable that determines the coeffients. In our case this is the
#' share of compliers.
#'
#' @return
#' @export
#'
#' @examples
thetahat_cbar_fn <- function(cbar,h,Y,W,C){
  omega <- sqrt((1/(2*h))*(abs(C - cbar) < h))
  omega[C == cbar] <- 0
  # if(cbar >= (1-h)){
  #    omega[ cbar >= (1-h)] <- omega / ((((1-cbar)/h)-1)*0.5)
  # }
  W_omega <- sweep(as.matrix(W),MARGIN = 1,omega, FUN = '*')
  Y_omega <- sweep(as.matrix(Y),MARGIN = 1,omega, FUN = '*')
  # thetahat_r <- solve(crossprod(W_omega,W_omega),crossprod(W_omega,Y_omega))
  thetahat_r <- MASS::ginv(crossprod(W_omega,W_omega))%*%crossprod(W_omega,Y_omega)
  return(thetahat_r)
}

#' Compute average treatment effects of the model
#'
#' $$ Y_{ig} = W_{ig}'theta(C_g) + U_{ig}$$
#'
#' where $U_{ig} = W_{ig}'(Theta_{ig} - E[Theta_{ig}|C_g])$.
#'
#' This function computes an estimate for $theta = E[theta(R_g)]$.
#' The function procedes in three steps: (1) For each city compute the estimated
#' coefficients $hat{theta}(R_g)$. (2) Compute a weighted averagem of the estimated
#' coeffients, depending on population size. (3) Adjust coefficients for the direct
#' and interaction effects to account for the total share of compliers.
#'
#' Notice that the kernel function is constructed with the whole data that
#' has individuals in each row rather than cities. We could do step (1) and (2)
#' in a single step by averaging over individuals rather than computing a weighted
#' average over cities. However, doing it this way is faster because we only need
#' to compute a kernel regression for each city and weight by the population in
#' that city. Otherwise we would need to repeat a lot of computations for individuals
#' in the same city.
#'
#' @param h Bandwidth choice to compute $theta(R_g)$.
#' @param Y Outcome variable.
#' @param W Set of regressors.
#' @param C The variable that determines the coeffients. In our case this is the
#' share of compliers.
#' @param D An individual take-up variable.
#' @param Z An individual offer variable.
#' @param C_g A vector with the share of compliers in each city.
#' @param rho_g A vector with the ratio of city size to the average city size.
#'
#' @return
#' @export
#'
#' @examples
thetahat_fn <- function(h,Y,W,C,D,Z,C_g,rho_g){

  # Step 1: "Localize" Compute kernel functions for each city.
  thetahat_vec     <- sapply(C_g,function(cbar) thetahat_cbar_fn(cbar,h,Y,W,C))

  # Step 2: "Average" Compute weighted average depending on city size.
  thetahat_vec_rho <- sweep(thetahat_vec, MARGIN = 2, rho_g, FUN = "*")
  thetahat         <- apply(thetahat_vec_rho,1,mean)

  # Step 3: "Adjust" Divide coefficients by total share of compliers to recover effect.
  prob_complier <- coef(lm(D~Z))[2]
  thetahat[c(2,4)] <- thetahat[c(2,4)]/prob_complier

  return(thetahat)

}





