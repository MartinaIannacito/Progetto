<<<<<<< HEAD
rm(list=ls())
install.packages("deSolve")
=======
rm(list = ls())
>>>>>>> 8d707e1a3b91b43d8d42d947bed503bf4ac83720
library(deSolve)
SEIR <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        N <- S + E + I + R
        dS <- -b * S * I / N #susceptible
        dE <- b * S * I / N - k * E #exposed
        dI <- k * E - g * I  #infectious
        dR <- g * I
        res <- c(dS, dE, dI, dR)
        list(res)
    })
}

## Parameters
<<<<<<< HEAD
transmission.rate=0.25
duration.of.infection=1
recovery.rate=1/duration.of.infection
duration.of.NON.infection=4
infection.rate=1/duration.of.NON.infection
parms  <- c(b = transmission.rate, g = recovery.rate, k = infection.rate)
=======
transmission_rate <- 0.25
duration_of_infection <- 10
recovery_rate <- 1 / duration_of_infection
duration_of_NON_infection <- 6
infection_rate <- 1/duration_of_NON_infection
parms <- c(b = transmission_rate, g = recovery_rate, k = infection_rate)
>>>>>>> 8d707e1a3b91b43d8d42d947bed503bf4ac83720

## vector of timesteps
ndays <- 365
times  <- seq(0, ndays, length = 1001)

## initial conditions
pop_size <- 100000
S0 <- pop_size-1
I0 <- 1
R0 <- 0
E0 <- 0
## Start values for steady state
xstart <- c(S = S0, E = E0, I = I0, R = R0)

## Solving
<<<<<<< HEAD
out <-  lsoda(xstart, times, SEIR, parms) 
=======
out <- lsoda(xstart, times, SEIR, parms) 
>>>>>>> 8d707e1a3b91b43d8d42d947bed503bf4ac83720

## Translate the output into a data.frame
out2 <- as.data.frame(out)

## Plotting
<<<<<<< HEAD
plot(out2$time, out2$S, col='blue', ylim=c(0, max(out[,-1])), type='l')
lines(out2$time, out2$I, col='red')
lines(out2$time, out2$R, col='orange')
lines(out2$time, out2$E, col="green")
=======
plot(out2$time, out2$S, col = 'blue', ylim = c(0, max(out[, -1])), type = 'l')
lines(out2$time, out2$I, col = 'red')
lines(out2$time, out2$R, col = 'orange')
lines(out2$time, out2$E, col = "green")
>>>>>>> 8d707e1a3b91b43d8d42d947bed503bf4ac83720

