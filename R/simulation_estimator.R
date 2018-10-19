
#' Simulate a dataset and compute IV.
#'
#' This function (1) Simulates the data. (2) Computes the IV coefficient and
#' (3) Outputs the difference between the estimated IV and the truth.
#'
#' #'
#' @param G
#' @param N_target
#'
#' @return $hat{theta}_{IV}-theta_{IV}$
#' @export
#'
#' @examples iv_sim(235,30000)
iv_sim <- function(G,N_target){

  # Simulate the dataset
  sim_param <- list(G = G,N_target=N_target)
  data      <- do.call(draw_sims,sim_param)
  # data      <- data %>% mutate(complier = (kappa > 0))

  # Compare the coefficients in different principal strata, and by city
  theta_comp      <- data %>% group_by(complier) %>%
                              summarise_at(vars(alpha,beta,gamma,delta),mean)
  theta_all       <- data %>% summarise_at(vars(alpha,beta,gamma,delta),mean)

  # Calculate the naive overall IV estimate for the entire sample
  ones      <- rep(1L, nrow(data))
  Z         <- as.matrix(data %>% transmute(ones = 1, z,saturation,saturation*z))
  X         <- as.matrix(data %>% transmute(ones = 1, d,dbar,d*dbar))
  Y         <- as.matrix(data[,'y'])
  theta_IV    <- solve(crossprod(Z, X),crossprod(Z, Y))

  results <- rep(NA,4)
  results[1] <- theta_IV[1] - theta_all$alpha
  results[2] <- theta_IV[2] - theta_comp[2,]$beta
  results[3] <- theta_IV[3] - theta_all$gamma
  results[4] <- theta_IV[4] - theta_comp[2,]$delta

  return(results)

}


#' Simulate a dataset and compute the random coefficients estimator.
#'
#' This function (1) Simulates the data. (2) Computes the random coefficient
#' (RC) estimator and
#' (3) Outputs the difference between the RC estimate and the truth.
#'
#' @param G
#' @param N_target
#' @param h
#'
#' @return
#' @export
#'
#' @examples
rcestimator_sim <- function(G,N_target,h){

  # Simulate the dataset
  sim_param <- list(G = G,N_target=N_target)
  data      <- do.call(draw_sims,sim_param)
  # data      <- data %>% mutate(complier = (kappa > 0))

  # Compare the coefficients in different principal strata, and by city
  theta_comp      <- data %>% group_by(complier) %>%
    summarise_at(vars(alpha,beta,gamma,delta),mean)
  theta_all       <- data %>% summarise_at(vars(alpha,beta,gamma,delta),mean)

  # Compute relevant quantities
  data <- data %>% filter(saturation %in% c(0.25,0.5,0.75,1)) # Only look at cities with saturation in the interior of 0 and 1 - why?
  data <- data %>% mutate(Chat = dbar/saturation,
                          ones = 1L,
                          z_dbar = z*dbar)

  # City data
  city_data <- data %>% group_by(city) %>%
    summarise_at(vars(saturation,dbar,complier),mean) %>%
    mutate(n_g = n())
  city_data <- city_data %>% mutate(rho = n_g/mean(n_g),
                                    Chat = dbar / saturation)

  # Prepare inputs for the estimator
  W <- as.matrix(subset(data,select = c('ones','z','dbar','z_dbar')))
  Y <- data$y
  C <- data$Chat
  # h <- 0.05
  D <- data$d
  Z <- data$z
  C_g <- city_data$Chat
  rho_g <- city_data$rho

  thetahat <- thetahat_fn(h,Y,W,C,D,Z,C_g,rho_g)

  results <- rep(NA,4)
  results[1] <- thetahat[1] - theta_all$alpha
  results[2] <- thetahat[2] - theta_comp[2,]$beta
  results[3] <- thetahat[3] - theta_all$gamma
  results[4] <- thetahat[4] - theta_comp[2,]$delta

  return(results)

}
