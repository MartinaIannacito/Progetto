library(deSolve)
library(ggplot2)

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    dS_0 <- l * S_1 - (mu + c) * S_0 - beta_0 * I_0 * S_0 + p_0 * E_0 + delta_1 * V_0 
    dE_0 <- l * E_1 + beta_0 * S_0 * I_0 - (mu + c + sigma_0 + p_0) * E_0
    dI_0 <- sigma_0 * E_0 + epsilon * I_1 - (mu + c + alpha) * I_0
    dV_0 <- l * V_1 - (mu + c + delta_1) * V_0
    
    dS_1 <- A - beta * S_1 * I_0 - beta_1 * S_1 * I_1 - (d + nu + l) * S_1 + p_1*E_1 + delta_2 * V_1
    dE_1 <- beta * S_1 * I_0 + beta_1 * S_1 * I_1 - (d + nu + sigma_1 + l + p_1) * E_1 
    dI_1 <- sigma_1 * E_1 - (d + epsilon + k) * I_1
    dV_1 <- nu * (S_1 + E_1) - (d + l + delta_2) * V_1
    
    dS_h <- H - mu_h * S_h - lambda_1 * S_h * I_0 - lambda_2 * S_h * I_1 + delta_1h * E_h + delta_2h * V_h
    dE_h <- lambda_1 * S_h * I_0 + lambda_2 * S_h * I_1 - (mu_h + sigma_h + delta_1h + nu_h) * E_h
    dI_h <- sigma_h * E_h - (mu_h + alpha_h) * I_h
    dV_h <- nu_h * E_h - (mu_h + delta_2h) * V_h
    
    return(list(c(dS_0, dE_0, dI_0, dV_0, dS_1, dE_1, dI_1, dV_1, dS_h, dE_h, dI_h, dV_h)))
  })
}

pars <- c(
  A = 7.7 * 10^(5), # month^(-1)
  l = 0.014,
  mu = 0.24,
  c = 0.06,
  delta_1 = 0.5,
  delta_2 = 0.5,
  sigma_0 = 0.35,
  sigma_1 = 0.37,
  p_0 = 0.35,
  p_1 = 0.37,
  epsilon = 0.1,
  alpha = 1,
  beta_0 = 8*10^(-6),
  beta = 4*10^(-6),
  beta_1 = 3.2*10^(-7),
  d = 0.11,
  nu = 0.133,
  k = 0.79,
  h = 10^6,
  mu_h = 4.6*10^(-3),
  lambda_1 = 3.6*10^(-9), # month^(-1)
  lambda_2 = 4.8*10^(-10),
  delta_1h = 0.33,
  delta_2h = 1,
  sigma_h = 0.33,
  nu_h = 0.328,
  alpha_h = 1,
  H = 10 ^ 6
)/ 12



init <- c(S_0 = 4 * 10 ^ 4, E_0 = 4 * 10 ^ 2, I_0 = 2.04 * 10 ^ 3, V_0 = 0,
          S_1 = 2.4 * 10 ^ 6, E_1 = 2.9 * 10 ^ 4, I_1 = 2 * 10 ^ 4, V_1 = 6 * 10 ^ 5,
          S_h = 7.988 * 10 ^ 7, E_h = 7.13 * 10 ^ 2, I_h = 3.87 * 10 ^ 2, V_h = 6 * 10 ^ 5)

times <- seq(1, 100*12, by = 1)
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

# Dogs

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") +
  myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_d), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_d), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme

# Humans

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_h), color = "blue") +
  myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E_h), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I_h), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, V_h), color = "green") + myTheme

############################

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S_d), color = "blue") +
  geom_line(mapping = aes(time, E_d), color = "orange") +
  geom_line(mapping = aes(time, I_d), color = "red") + 
  geom_line(mapping = aes(time, R_d), color = "green") + myTheme

p<-as.list(c(init[1],12*pars))

R0<-with(p,(beta*p[[1]]*sigma*gamma)/((m+k+sigma)*(m+mu)))

I_d_star<-with(p,(m+sigma+k)*(m+lambda+k)*m*(R0-1)/(beta*(m*(m+lambda+k)+sigma*gamma*(m+lambda))))

denominatore<-with(p,(m_h+lambda_h)*(m_h*(m_h+k_h+sigma_h)+beta_h*I_d_star*(m+k+sigma*gamma))-beta_h*I_d_star*lambda_h*k_h)

denominatore2<-with(p,(m_h+lambda_h)*(m_h*(m_h+k_h+sigma_h)+beta_h*I_d_star*(m_h+k_h+sigma_h*gamma_h))-beta_h*I_d_star*lambda_h*k_h)


E_h_star<-with(p,(beta_h*B*(m_h+lambda_h)*I_d_star)/(denominatore2))

I_h_star<-with(p,sigma_h*gamma_h*E_h_star/((m_h+mu_h)))
