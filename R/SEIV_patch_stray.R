library(deSolve)
library(ggplot2)

n_patch <- 2

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    ind <- matrix(1:(n_patch*12), ncol = n_patch, byrow = TRUE)
    S_d <- state[ind[1, ]]
    E_d <- state[ind[2, ]]
    I_d <- state[ind[3, ]]
    V_d <- state[ind[4, ]]
    
    S_s <- state[ind[5, ]]
    E_s <- state[ind[6, ]]
    I_s <- state[ind[7, ]]
    V_s <- state[ind[8, ]]
    
    S_h <- state[ind[9, ]]
    E_h <- state[ind[10, ]]
    I_h <- state[ind[11, ]]
    V_h <- state[ind[12, ]]
    
    dS_d <- A_d + lambda_d * V_d + sigma_d * (1 - gamma_d) * E_d - beta_dd * S_d * I_d - (m_d + k_d) * S_d + phi_S %*% S_d - l * S_d - beta_ds * S_d * I_s
    dE_d <- beta_dd * S_d * I_d - sigma_d * (1 - gamma_d) * E_d - sigma_d * gamma_d * E_d - (m_d + k_d) * E_d + phi_E %*% E_d - l * E_d + beta_ds * S_d * I_s
    dI_d <- sigma_d * gamma_d * E_d - (m_d + mu_d) * I_d + phi_I %*% I_d - l_i * I_d
    dV_d <- k_d * (S_d + E_d) - (m_d + lambda_d) * V_d + phi_V %*% V_d - l * V_d
    
    dS_s <- A_s + lambda_s * V_s + sigma_s * (1 - gamma_s) * E_s - beta_ss * S_s * I_s - (m_s + k_s) * S_s + rho_S %*% S_s + l * S_d - beta_ds * S_s * I_d
    dE_s <- beta_ss * S_s * I_s - sigma_s * (1 - gamma_s) * E_s - sigma_s * gamma_s * E_s - (m_s + k_s) * E_s + rho_E %*% E_s + l * E_d + beta_ds * S_s * I_d
    dI_s <- sigma_s * gamma_s * E_s - (m_s + mu_s) * I_s + rho_I %*% I_s + l_i * I_d
    dV_s <- k_s * (S_s + E_s) - (m_s + lambda_s) * V_s + rho_V %*% V_s + l * V_d
    
    dS_h <- B + lambda_h * V_h + sigma_h * (1 - gamma_h) * E_h - beta_hd * S_h * I_d - m_h  * S_h + psi_S %*% S_h - beta_hs * S_h * I_d
    dE_h <- beta_hd * S_h * I_d - sigma_h * (1 - gamma_h) * E_h - sigma_h * gamma_h * E_h - (m_h + k_h) * E_h + psi_E %*% E_h + beta_hs * S_h * I_d
    dI_h <- sigma_h * gamma_h * E_h - (m_h + mu_h) * I_h + psi_I %*% I_h
    dV_h <- k_h * E_h - (m_h + lambda_h) * V_h + psi_V %*% V_h
    
    return(list(c(dS_d, dE_d, dI_d, dV_d, dS_s, dE_s, dI_s, dV_s, dS_h, dE_h, dI_h, dV_h)))
  })
}

phi <- matrix(c(0, 0, 0, 0), nrow = 2) # assume domestic dogs don't move
rho <- matrix(c(0, 0.05, 0.6, 0), nrow = 2)
diag(rho) <- -c(0.05, 0.6)
psi <- matrix(c(0, 0.32, 0.2, 0), nrow = 2)
diag(psi) <- -c(0.32, 0.2)

pars <- list(
  A_d = c(3, 5) * 10 ^ 5,
  A_s = c(3, 5) * 10 ^ 5, # stray dogs can reproduce
  B = c(8.797, 4.101) * 10 ^ 5,
  lambda_d = rep(0.33, 2),
  lambda_s = rep(0.33, 2),
  lambda_h = c(2, 2),
  sigma_d = c(10, 10),
  sigma_s = c(10, 10),
  sigma_h = c(6, 6),
  gamma_d = rep(0.4, 2) * 12,
  gamma_s = rep(0.4, 2) * 12,
  gamma_h = rep(0.475, 2) * 12,
  m_d = rep(0.1, 2), # reduced mortality of domestics
  m_s = rep(0.345, 2),
  m_h = rep(0.00662, 2),
  k_d = c(0.09, 0),
  k_s = c(0.01, 0), # low vaccination rate of dogs
  k_h = rep(0.5, 2),
  mu_d = c(1, 1),
  mu_s = c(1, 1),
  mu_h = c(1, 1),
  beta_dd = c(2.45, 2.2) * 10 ^ (-7), # lower transmission rate between domestics
  beta_ds = c(2.45, 2.2) * 10 ^ (-6),
  beta_ss = c(2.45, 2.2) * 10 ^ (-6),
  beta_hd = c(2.2, 1.4) * 10 ^ (-11),
  beta_hs = c(2.2, 1.4) * 10 ^ (-11),
  phi_S = phi,
  phi_E = phi,
  phi_I = phi,
  phi_V = phi,
  psi_S = psi,
  psi_E = psi,
  psi_I = matrix(0, nrow = 2, ncol = 2),
  psi_V = psi,
  rho_S = rho,
  rho_E = rho,
  rho_I = rho,
  rho_V = rho,
  l = rep(0.014, 2),
  l_i = c(0.1, 0.1)
)

pars <- lapply(pars, "/", 12)

init <- c(
  S_d_1 = 5 * 10 ^ 5,
  S_d_2 = 2.4 * 10 ^ 6,
  E_d_1 = 10 ^ 3,
  E_d_2 = 2.9 * 10 ^ 4,
  I_d_1 = 10 ^ 3,
  I_d_2 = 5 * 10 ^ 4,
  V_d_1 = 10 ^ 3,
  V_d_2 = 2 * 10 ^ 4,
  
  S_s_1 = 6 * 10 ^ 3,
  S_s_2 = 4 * 10 ^ 4,
  E_s_1 = 20,
  E_s_2 =  4 * 10 ^ 2,
  I_s_1 = 30,
  I_s_2 = 2.04 * 10 ^ 3,
  V_s_1 = 0,
  V_s_2 = 0,
  
  S_h_1 = 7.1 * 10 ^ 7,
  S_h_2 = 3.8 * 10 ^ 7,
  E_h_1 = 10,
  E_h_2 = 50,
  I_h_1 = 1,
  I_h_2 = 20,
  V_h_1 = 6 * 10 ^ 2,
  V_h_2 = 6 * 10 ^ 2
)

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
  plot.title = element_text(hjust = 0.5, size = 25),
  legend.text = element_text(size = 13),
  legend.title= element_text(size = 15))

# Dogs

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_d_1), color = "blue") + 
  geom_line(mapping = aes(time, I_d_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_s_1), color = "blue") + 
  geom_line(mapping = aes(time, I_s_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d_1), color = "blue") + 
  geom_line(mapping = aes(time, S_d_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_s_1), color = "blue") + 
  geom_line(mapping = aes(time, S_s_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_h_1), color = "blue") + 
  geom_line(mapping = aes(time, I_h_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_h_1), color = "blue") + 
  geom_line(mapping = aes(time, E_h_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_h_1), color = "blue") + 
  geom_line(mapping = aes(time, S_h_2), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, V_h_1), color = "blue") + 
  geom_line(mapping = aes(time, V_h_2), color = "red") + myTheme

################################################################################
# PLOTS WITH ALL HUMANS
################################################################################

data <- SEIR_out[, c("time", "E_h_1", "E_h_2", "I_h_1", "I_h_2")]
#data[, c("S_h", "V_h")] <- log10(data[, c("S_h", "V_h")])*100

#labels <- 10^(seq(0, 9, by = 1))
#breaks <- log10(labels)*100

SEIR_out_long <- reshape2::melt(as.data.frame(data), id = "time")

ggplot(data = SEIR_out_long, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  scale_color_discrete(name = "class", labels =
      c(latex2exp::TeX("$E_h^{Hebei}$"), "2", "3", "4")) +
  myTheme + ylab(latex2exp::TeX("$E_h$, $I_h$")) +
   theme(legend.position="bottom",
    axis.title.y.right = element_text(margin = margin(l = 20))) + 
  xlab("time (months)")
