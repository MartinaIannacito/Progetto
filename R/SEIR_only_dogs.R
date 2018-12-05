library(deSolve)
library(ggplot2)

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    dS_d <- A + lambda * R_d + sigma * (1 - gamma) * E_d - beta * S_d * I_d - (m + k) * S_d
    dE_d <- beta * S_d * I_d - sigma * (1 - gamma) * E_d - sigma * gamma * E_d - (m + k) * E_d
    dI_d <- sigma * gamma * E_d - (m + mu) * I_d
    dR_d <- k * (S_d + E_d) - (m + lambda) * R_d
    
    return(list(c(dS_d, dE_d, dI_d, dR_d)))
  })
}

pars <- c(
  A = 3 * 10^(6), # month^(-1)
  lambda = 1, # month^(-1)
  sigma = 6, # monrt^(-1)
  gamma = 0.4, # same
  m = 0.08, # same
  mu = 1, # same
  k = 0.09, # same
  beta = 1.58 * 10^(-7)
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

init <- c(S_d = 3.5 * 10^7, E_d = 2 * 10^5, I_d = 1 * 10^5, R_d = 2 * 10^5)
#init <- c(S_d = 3.5 * 10^2, E_d = 0, I_d = 1, R_d = 0)
times <- seq(1, 200*12, by = 1)
SEIR_out <- ode(init, times, SEIR, pars)

myTheme <- theme(axis.text=element_text(size=20),
  axis.title=element_text(size=25),
  axis.title.y = element_text(margin = margin(r=20)),
  axis.title.x = element_text(margin = margin(t=20)),
  axis.ticks = element_line(size = .7),
  axis.ticks.length = unit(.3, "cm"),
  panel.background = element_rect(fill="white"),
  panel.grid.major = element_line(colour = "grey90"),
  panel.grid.minor = element_line(colour = "grey90"),
  axis.line = element_line(color = "black", size = .7),
  plot.title = element_text(hjust = 0.5, size = 25))

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") +
  myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_d), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_d), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme

############################

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") +
  geom_line(mapping = aes(time, E_d), color = "orange") +
  geom_line(mapping = aes(time, I_d), color = "red") + 
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme
