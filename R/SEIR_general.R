#install.packages("deSolve")
#install.packages("ggplot2")
library(deSolve)
library(ggplot2)

SEIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    dS <- b * (1 - S) - R_0 * S * I
    dE <- R_0 * S * I - (p + b) * E
    dI <- p * E - I
    dR <- I - b * R 
    
    return(list(c(dS, dE, dI, dR)))
  })
}

## only endemic equilibrium

pars <- c(
  R_0 <- 10,
  p <- 0.5,
  b <- 1
) 

init <- c(S = 9/10, E = 0, I = 1/10, R = 0)
times <- seq(1, 12*1.5, by = 0.1)
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


ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S), color = "blue") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, E), color = "orange") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, I), color = "red") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, R), color = "green") + myTheme

ggplot(data = as.data.frame(SEIR_out)) +
  geom_line(mapping = aes(time, S/10), color = "blue") +
  geom_line(mapping = aes(time, E), color = "orange") +
  geom_line(mapping = aes(time, I), color = "red") +
  geom_line(mapping = aes(time, R), color = "green") + myTheme


################################################################################
# PLOT WITH ALL CLASSES
################################################################################

data <- SEIR_out[, c("time", "S", "E", "I", "R")]
data[, "S"] <- data[, "S"]/2

SEIR_out_long <- reshape2::melt(as.data.frame(data), id = "time")

ggplot(data = SEIR_out_long, aes(x = time, y = value, colour = variable)) +
  geom_line() + scale_color_discrete(name = "class") +
  myTheme + ylab("E, I, R") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*2, name = "S")) + 
  theme(legend.position="bottom", axis.title.y.right =
      element_text(margin = margin(l = 20))) + 
  xlab("time")

################################################################################
# ONLY DISEASE-FREE EQUILIBRIUM
################################################################################

pars2 <- c(
  R_0 <- 0.5,
  p <- 0.5,
  b <- 1
)

SEIR_out2 <- ode(init, times, SEIR, pars2)

data2 <- SEIR_out2[, c("time", "S", "E", "I", "R")]
data2[, "S"] <- data2[, "S"]/10

SEIR_out_long2 <- reshape2::melt(as.data.frame(data2), id = "time")

ggplot(data = SEIR_out_long2, aes(x = time, y = value, colour = variable)) +
  geom_line() + scale_color_discrete(name = "class") +
  myTheme + ylab("E, I, R") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*10, name = "S")) + 
  theme(legend.position="bottom", axis.title.y.right =
      element_text(margin = margin(l = 20))) + 
  xlab("time")
