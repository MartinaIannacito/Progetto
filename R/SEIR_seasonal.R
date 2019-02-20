library(deSolve)
library(ggplot2)

SEIR <- function(t, state, pars) {
  with(as.list(c(t, state, pars)), {
    
    beta <- function(t) a * (1 + b * sin(pi / 6 * t + 5.5))
    beta_h <- function(t) a_1 * (1 + b_1 * sin(pi / 6 * t + 5.5))
    
    dS_d <- A + lambda * R_d + sigma * (1 - gamma) * E_d - m * S_d - beta(t) * S_d * I_d - k * S_d
    dE_d <- beta(t) * S_d * I_d - sigma * (1 - gamma) * E_d - k * E_d - m * E_d - sigma * gamma * E_d
    dI_d <- sigma * gamma * E_d - (m + mu) * I_d
    dR_d <- k * (S_d + E_d) - (m + lambda) * R_d
    
    dS_h <- B + lambda_h * R_h + sigma_h * (1 - gamma_h) * E_h - beta_h(t) * S_h * I_d - m_h  * S_h
    dE_h <- beta_h(t) * S_h * I_d - sigma_h * (1 - gamma_h) * E_h - sigma_h * gamma_h * E_h - (m_h + k_h) * E_h
    dI_h <- sigma_h * gamma_h * E_h - (m_h + mu_h) * I_h
    dR_h <- k_h * E_h - (m_h + lambda_h) * R_h
    
    return(list(c(dS_d, dE_d, dI_d, dR_d, dS_h, dE_h, dI_h, dR_h)))
  })
}

# months ^ -1
pars <- c(
  A = 2.34 * 10 ^ 5,
  lambda = 1 / 6,
  sigma = 1 / 1.045,
  gamma = 0.49,
  m = 0.0064,
  a = 9.9 * 10 ^ (-8),
  b = 0.41,
  mu = 1,
  k = 0.09,
  B = 1.34 * 10^6,
  lambda_h = 1 / 6,
  sigma_h = 1 / 2,
  gamma_h = 0.5,
  m_h = 0.00057,
  a_1 = 2.41 * 10 ^ (-11),
  b_1 = 0.23,
  k_h = 0.54,
  mu_h = 1
)

init <- c(S_d = 3.3 * 10 ^ 7, E_d = 2.2 * 10 ^ 4, I_d = 1.1 * 10 ^ 4,
  R_d = 3.3 * 10 ^ 6, S_h = 1.29 * 10 ^ 9, E_h = 178, I_h = 89, R_h = 6 * 10 ^ 7)
times <- seq(1, 150, by = 0.1)
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

################################################################################
# LONG TIME MODEL
################################################################################

times_2 <- seq(1, 1500, by = 0.1)
SEIR_out_2 <- ode(init, times_2, SEIR, pars)

ggplot(data = as.data.frame(SEIR_out_2)) +
  geom_line(mapping = aes(time, I_h), color = "red") + myTheme

################################################################################
# DIFFERENT INITIAL CONDITIONS
################################################################################

S_d_0_vect <- c(2.5, 2.7, 3, 3.2, 3.5)*10^7

ode_system_1 <- function(S_d_0) {
  init["S_d"] <- S_d_0
  times <- seq(1, 300, by = 1)
  ode(init, times, SEIR, pars)
}

res_1 <- lapply(S_d_0_vect, ode_system_1)

I_h_df_1 <- as.data.frame(lapply(res_1, `[`, , "I_h"))
colnames(I_h_df_1) <- S_d_0_vect
I_h_df_1 <- cbind(I_h_df_1, time = res_1[[1]][,"time"])
I_h_df_long_1 <- reshape2::melt(I_h_df_1, id = "time")

ggplot(data = I_h_df_long_1, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  myTheme

#-------------------------------------------------------------------------------

I_d_0_vect <- c(5, 8, 10, 12, 15)*10^3

ode_system_2 <- function(I_d_0) {
  init["I_d"] <- I_d_0
  times <- seq(1, 300, by = 1)
  ode(init, times, SEIR, pars)
}

res_2 <- lapply(I_d_0_vect, ode_system_2)

I_h_df_2 <- as.data.frame(lapply(res_2, `[`, , "I_h"))
colnames(I_h_df_2) <- I_d_0_vect
I_h_df_2 <- cbind(I_h_df_2, time = res_2[[1]][,"time"])
I_h_df_long_2 <- reshape2::melt(I_h_df_2, id = "time")

ggplot(data = I_h_df_long_2, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  myTheme
