library(MASS)
library(dplyr)
library(tibble)

#------------------------------------------------------------------------------
draw_sims <- function(G=235,N_target=30000) {
  # G = number of groups, N_target = target number of total observations
  N_g <- as.matrix(rep(round(N_target/G),G)) # rpois(G, N_target / G)#  random draw for size of group g
  N <- sum(N_g) # total number of observations
  saturation <- rep(c(0, 0.25, 0.5, 0.75, 1), each = G / 5)
  noffers <- floor(saturation * N_g)

  m <- c(kappa = 0, alpha = 0.2, beta = 1, gamma = -0.2, delta = -0.1)
  #r_theta <- c(-0.15, 0.3, 0.2, 0.6)
  #r_epsilon <- c(-0.3, 0.3, 0.2, -0.3)
  r_theta <- c(0, 0.3, 0.5, -0.1) # This is the correlation between the coefficients
  r_epsilon <- c(0, 0.6, 0.3, -0.1) # This is the correlation between the epsilons?

  # sd_theta <- diag(c(kappa = 0.1, alpha = 0.025, beta = 0.25, gamma = 0.05, delta = 0.05))
  # sd_epsilon <- diag(c(kappa = 2, alpha = 0.1, beta = 1, gamma = 0.2, delta = 0.1))
  sd_theta <- diag(c(kappa = 0.1, alpha = 0.025, beta = 0.5, gamma = 0.2, delta = 0.3))
  sd_epsilon <- diag(c(kappa = 0.3, alpha = 0.1, beta = 1, gamma = 0.1, delta = 0.05))

  R_theta <- diag(length(m))
  R_theta[1, 2:5] <- R_theta[2:5, 1] <- r_theta
  Sigma_theta <- sd_theta %*% R_theta %*% sd_theta
  colnames(Sigma_theta) <- rownames(Sigma_theta) <- names(m)
  theta_g <- MASS::mvrnorm(G, rep(0, ncol(Sigma_theta)), Sigma_theta) # This generates the group-level random coefficients
  # using a multivariate normal distribution with mean zero and Sigma_theta as the covariance matrix

  R_epsilon <- diag(length(m))
  R_epsilon[1, 2:5] <- R_epsilon[2:5, 1] <- r_epsilon
  Sigma_epsilon <- sd_epsilon %*% R_epsilon %*% sd_epsilon
  colnames(Sigma_epsilon) <- rownames(Sigma_epsilon) <- names(m)
  epsilon_ig <- MASS::mvrnorm(sum(N_g), rep(0, ncol(Sigma_epsilon)), Sigma_epsilon) # This generates the idiosyncractic elements of the
  # random coefficients using a multivariate normal distribution with mean zero and Sigma_epsilon as the covariance matrix

  m_mat <- matrix(m, sum(N_g), 5, byrow = TRUE) # Create matrix with the parameters in each row
  city <- rep(1:G, N_g)
  theta_g_mat <- theta_g[city,]
  theta_ig <- m_mat + theta_g_mat + epsilon_ig
  # the random coefficients are a sum of the average coefficient and two random parts (group-level, individual-level)

  get_z_vector <- function(noffers, n) {
    c(rep(1, noffers), rep(0, n - noffers))
    #rbinom(n,1,noffers/n)
  }
  z <- as.numeric(unlist(mapply(get_z_vector, noffers, N_g)))

  # Generate the dataset
  dat <- data.frame(city,saturation = saturation[city],theta_ig,z)
  dat <- tibble::as.tibble(dat) # wrapper on the dataframe
  dat <- dat %>% mutate(d = z * (kappa > 0))
  dat <- dat %>% mutate(complier = (kappa > 0))

  # The following adds a simple strategic interaction:
  #dat <- dat %>% mutate(d = z * (kappa + (saturation == 0.25) * 0.5 > 0))

  dbar <- dat %>% group_by(city) %>% dplyr::summarize(dbar = mean(d))
  sc <- dat %>% group_by(city) %>% dplyr::summarize(sc = mean(kappa > 0)) # Is this is the proportion of compliers?
  dat <- full_join(dat, dbar, by = 'city')
  dat <- full_join(dat, sc, by = 'city')
  dat <- dat %>% mutate(y = alpha + beta * d + gamma * dbar + delta * d * dbar) # Add the outcomes data
  return(dat)
}


# Code to test for strategic interactions
extrapolate_takeup <- function(x, sat, compliers100, never100) {
  Ksat <- sat * (compliers100 + never100)
  mean(phyper(x, m = compliers100, n = never100, k = Ksat))
}
vextrapolate_takeup <- Vectorize(extrapolate_takeup, 'x')

get_test_stat <- function(dat_not100, dat100) {
  compliers100 <- dat100$nTakeup
  never100 <- with(dat100, nI - nTakeup)
  maxCompliers <- max(compliers100)

  get_KS <- function(sat) {
    Hhat <- ecdf(subset(dat_not100, saturation == sat)$nTakeup)
    Htilde <- function(x) vextrapolate_takeup(x, sat, compliers100, never100)
    xseq <- 0:max(knots(Hhat), maxCompliers)
    D <- Hhat(xseq) - Htilde(xseq)
    return(max(abs(D)))
  }
  vget_KS <- Vectorize(get_KS)
  max(vget_KS(c(0.25, 0.5, 0.75)))
}

get_test_stat_boot <- function(dat100) {
  compliers <- data.frame(nI = dat100$nI, nCompliers = dat100$nTakeup)
  n <- nrow(compliers)
  saturations <- c(0.25, 0.5, 0.75, 1)
  nboot <- n * length(saturations)
  dat_boot <- compliers[sample(1:n, size = nboot, replace = TRUE),]
  dat_boot$saturation <- rep(saturations, each = n)
  dat_boot$nOffered <- with(dat_boot, floor(saturation * nI))
  dat_boot$nTakeup <- with(dat_boot, rhyper(nn = nboot,
                                            m = nCompliers,
                                            n = nI - nCompliers,
                                            k = nOffered))
  dat_boot$nCompliers <- NULL
  dat100_boot <- subset(dat_boot, saturation == 1)
  dat_not100_boot <- subset(dat_boot, saturation != 1)
  get_test_stat(dat_not100_boot, dat100_boot)
}

testplot <- function(dat) {
  dat100 <- subset(dat, saturation == 1)
  compliers100 <- dat100$nTakeup
  never100 <- with(dat100, nI - nTakeup)
  dat75 <- subset(dat, saturation == 0.75)
  dat50 <- subset(dat, saturation == 0.50)
  dat25 <- subset(dat, saturation == 0.25)

  xseq <- 0:max(compliers100 + never100)
  H <- function(sat) vextrapolate_takeup(xseq, sat, compliers100, never100)
  par(mfrow = c(1, 3))

  with(dat25, plot(ecdf(nTakeup), main = '25%'))
  points(xseq, H(0.25), type = 's', col = 'red')

  with(dat50, plot(ecdf(nTakeup), main = '50%'))
  points(xseq, H(0.5), type = 's', col = 'red')

  with(dat75, plot(ecdf(nTakeup), main = '75%'))
  points(xseq, H(0.75), type = 's', col = 'red')

  par(mfrow = c(1, 1))
}


bootplot <- function(dat100) {
  compliers <- data.frame(nI = dat100$nI, nCompliers = dat100$nTakeup)
  n <- nrow(compliers)
  saturations <- c(0.25, 0.5, 0.75, 1)
  nboot <- n * length(saturations)
  dat_boot <- compliers[sample(1:n, size = nboot, replace = TRUE),]
  dat_boot$saturation <- rep(saturations, each = n)
  dat_boot$nOffered <- with(dat_boot, floor(saturation * nI))
  dat_boot$nTakeup <- with(dat_boot, rhyper(nn = nboot,
                                            m = nCompliers,
                                            n = nI - nCompliers,
                                            k = nOffered))
  dat_boot$nCompliers <- NULL
  testplot(dat_boot)
}

#------------------------------------------------------------------------------
