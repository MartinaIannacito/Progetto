library(deSolve)
library(ggplot2)

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    dS_d <- A + lambda * R_d + sigma * (1 - gamma) * E_d - beta * S_d * I_d - 
      (m + k) * S_d
    dE_d <- beta * S_d * I_d - sigma * (1 - gamma) * E_d - sigma * gamma * E_d - 
      (m + k) * E_d
    dI_d <- sigma * gamma * E_d - (m + mu) * I_d
    dR_d <- k * (S_d + E_d) - (m + lambda) * R_d
    
    dS_h <- B + lambda_h * R_h + sigma_h * (1 - gamma_h) * E_h -
      beta_h * S_h * I_d - m_h  * S_h
    dE_h <- beta_h * S_h * I_d - sigma_h * (1 - gamma_h) * E_h -
      sigma_h * gamma_h * E_h - (m_h + k_h) * E_h
    dI_h <- sigma_h * gamma_h * E_h - (m_h + mu_h) * I_h
    dR_h <- k_h * E_h - (m_h + lambda_h) * R_h
    
    return(list(c(dS_d, dE_d, dI_d, dR_d, dS_h, dE_h, dI_h, dR_h)))
  })
}

pars <- c(
  A = 3 * 10^6, # month^(-1) <-
  lambda = 1, # month^(-1) <-
  sigma = 6, # monrt^(-1)
  gamma = 12*0.4, # pure number
  m = 0.08, # same
  mu = 1, # same
  k = 0.09, # same
  beta = 1.58 * 10^(-7),
  B = 1.54 * 10^7,
  gamma_h = 0.4*12, # pure number
  lambda_h = 1,
  sigma_h = 6,
  m_h = 0.0066,
  beta_h = 2.29 * 10^(-12),
  k_h = 0.54,
  mu_h = 1
) / 12

# pars <- c(
#   A = 2.34 * 10^5, # month^(-1)
#   lambda = 1/6, # month^(-1)
#   sigma = 1/1.045, # monrt^(-1)
#   gamma = 0.49, # same
#   m = 0.0064, # same
#   mu = 1, # same
#   k = 0.09, # same
#   beta = 
# )

init <- c(S_d = 3.5 * 10 ^ 7, E_d = 2 * 10 ^ 5, I_d = 1 * 10 ^ 5,
  R_d = 2 * 10 ^ 5, S_h = 1.29 * 10 ^ 9, E_h = 250, I_h = 89, R_h = 2 * 10 ^ 5)
#init <- c(S_d = 3.5 * 10^2, E_d = 0, I_d = 1, R_d = 0)
times <- seq(1, 50 * 12, by = 1)
SEIR_out <- ode(init, times, SEIR, pars)

myTheme <- theme(axis.text = element_text(size = 20),
  axis.title = element_text(size = 25),
  axis.title.y = element_text(margin = margin(r = 20)),
  axis.title.x = element_text(margin = margin(t = 20)),
  axis.ticks = element_line(size = .7),
  axis.ticks.length = unit(.3, "cm"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "grey90"),
  panel.grid.minor = element_line(colour = "grey90"),
  axis.line = element_line(color = "black", size = .7),
  plot.title = element_text(hjust = 0.5, size = 25))

# Dogs

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_d), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_d), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme

# Humans

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_h), color = "blue") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_h), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_h), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, R_h), color = "green") + myTheme

# Dogs, 4 in 1 plot

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") +
  geom_line(mapping = aes(time, E_d), color = "orange") +
  geom_line(mapping = aes(time, I_d), color = "red") + 
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme

p <- as.list(c(init[1], pars))

R0 <- with(p, (beta * p[[1]] * sigma * gamma) / ((m + k + sigma) * (m + mu)))

I_d_star <- with(p, (m + sigma + k) * (m + lambda + k) * m * (R0 - 1) / 
    (beta * (m * (m + lambda + k) + sigma * gamma * (m + lambda))))

denominatore <- with(p, (m_h + lambda_h) * (m_h * (m_h + k_h + sigma_h) +
    beta_h * I_d_star * (m + k + sigma * gamma)) - 
    beta_h * I_d_star * lambda_h * k_h)

denominatore2 <- with(p, (m_h + lambda_h) * (m_h * (m_h + k_h + sigma_h) +
    beta_h * I_d_star * (m_h + k_h + sigma_h * gamma_h)) -
    beta_h * I_d_star * lambda_h * k_h)

E_h_star <- with(p, (beta_h * B * (m_h + lambda_h) * I_d_star) / denominatore)

I_h_star <- with(p, sigma_h * gamma_h * E_h_star / (m_h + mu_h))

################################################################################
# SOLVE SYSTEM WITH DIFFERENT INITIAL CONDITIONS
################################################################################

S_d_vect <- c(4, 3, 2, 1, 0.5, 0.3) *10^7

ode_system <- function(S_d_init) {
  init <- c(S_d = S_d_init, E_d = 2 * 10 ^ 5, I_d = 1 * 10 ^ 5,
    R_d = 2 * 10 ^ 5, S_h = 1.29 * 10 ^ 9, E_h = 250, I_h = 89, R_h = 2 * 10 ^ 5)
  #init <- c(S_d = 3.5 * 10^2, E_d = 0, I_d = 1, R_d = 0)
  ode(init, times, SEIR, pars)
}

res <- lapply(S_d_vect, ode_system)
cols <- grDevices::rainbow(length(res))

infected_h <- lapply(seq_along(res), function(i)
  geom_line(data = as.data.frame(res[[i]]),
    mapping = aes(time, I_h), color = cols[i]))

add_geoms <- function(l) {
  g <- ggplot()
  for (geom in l) g <- g + geom
  g
}

add_geoms(infected_h) + guides(color = "colorbar") + myTheme
