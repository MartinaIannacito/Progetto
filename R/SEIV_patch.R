library(deSolve)
library(ggplot2)

n_patch <- 2

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    ind <- matrix(1:(n_patch*8), ncol = n_patch, byrow = TRUE)
    S_d <- state[ind[1, ]]
    E_d <- state[ind[2, ]]
    I_d <- state[ind[3, ]]
    V_d <- state[ind[4, ]]
    
    S_h <- state[ind[5, ]]
    E_h <- state[ind[6, ]]
    I_h <- state[ind[7, ]]
    V_h <- state[ind[8, ]]
    
    dS_d <- A + lambda_d * V_d + sigma_d * (1 - gamma_d) * E_d - beta_d * S_d * I_d - (m_d + k_d) * S_d + phi_S %*% S_d
    dE_d <- beta_d * S_d * I_d - sigma_d * (1 - gamma_d) * E_d - sigma_d * gamma_d * E_d - (m_d + k_d) * E_d + phi_E %*% E_d
    dI_d <- sigma_d * gamma_d * E_d - (m_d + mu_d) * I_d + phi_I %*% I_d
    dV_d <- k_d * (S_d + E_d) - (m_d + lambda_d) * V_d + phi_V %*% V_d
    
    dS_h <- B + lambda_h * V_h + sigma_h * (1 - gamma_h) * E_h - beta_h * S_h * I_d - m_h  * S_h + psi_S %*% S_h
    dE_h <- beta_h * S_h * I_d - sigma_h * (1 - gamma_h) * E_h - sigma_h * gamma_h * E_h - (m_h + k_h) * E_h + psi_E %*% E_h
    dI_h <- sigma_h * gamma_h * E_h - (m_h + mu_h) * I_h + psi_I %*% I_h
    dV_h <- k_h * E_h - (m_h + lambda_h) * V_h + psi_V %*% V_h
    
    return(list(c(dS_d, dE_d, dI_d, dV_d, dS_h, dE_h, dI_h, dV_h)))
  })
}

phi <- matrix(c(0, 0.05, 0.6, 0), nrow = 2)
diag(phi) <- -rep(sum(phi), 2)
psi <- matrix(c(0, 0.32, 02, 0), nrow = 2)
diag(psi) <- -rep(sum(psi), 2)

pars <- list(
  A = c(3, 5) * 10 ^ 5,
  B = c(8.797, 4.101) * 10 ^ 5,
  lambda_d = rep(0.33, 2),
  lambda_h = c(2, 2),
  sigma_d = c(10, 10),
  sigma_h = c(6, 6),
  gamma_d = rep(0.4, 2) * 12,
  gamma_h = rep(0.475, 2) * 12,
  m_d = rep(0.345, 2),
  m_h = rep(0.00662, 2),
  k_d = c(0.09, 0),
  k_h = rep(0.5, 2),
  mu_d = c(1, 1),
  mu_h = c(1, 1),
  beta_d = c(2.45, 2.2) * 10 ^ (-6),
  beta_h = c(2.2, 1.4) * 10 ^ (-11),
  phi_S = phi,
  phi_E = phi,
  phi_I = phi,
  phi_V = phi,
  psi_S = psi,
  psi_E = psi,
  psi_I = matrix(0, nrow = 2, ncol = 2),
  psi_V = psi
)

pars <- lapply(pars, "/", 12)

init <- c(S_d = 3.5 * 10 ^ 7, E_d = 2 * 10 ^ 5, I_d = 1 * 10 ^ 5,
  R_d = 2 * 10 ^ 5, S_h = 1.29 * 10 ^ 9, E_h = 250, I_h = 89, R_h = 2 * 10 ^ 5,
  S_d = 3.5 * 10 ^ 7, E_d = 2 * 10 ^ 5, I_d = 1 * 10 ^ 5,
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


