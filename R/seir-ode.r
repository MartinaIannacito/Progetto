rm(list=ls())
library(deSolve)
SEIR <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        N <- S+E+I+R
        dS <- - b*S*I/N  #susceptible
        dE <- b*S*I/N - k*E #exposed
        dI <- k*E - g*I  #infectious
        dR <- g*I
        res <- c(dS, dE, dI, dR)
        list(res)
    })
}

## Parameters
transmission.rate=0.25
duration.of.infection=10
recovery.rate=1/duration.of.infection
duration.of.NON.infection=6
infection.rate=1/duration.of.NON.infection
parms  <- c(b = transmission.rate, g = recovery.rate, k = infection.rate)

## vector of timesteps
ndays=365
times  <- seq(0, ndays, length = 1001)

## initial conditions
pop.size=100000
S0=pop.size-1
I0=1
R0=0
E0=0
## Start values for steady state
xstart <- c(S = S0, E=E0, I= I0, R = R0)

## Solving
out <-  lsoda(xstart, times, SIR, parms) 

## Translate the output into a data.frame
out2 <- as.data.frame(out)

## Plotting
plot(out2$time, out2$S, col='blue', ylim=c(0, max(out[,-1])), type='l')
lines(out2$time, out2$I, col='red')
lines(out2$time, out2$R, col='orange')
lines(out2$time, out2$E, col="green")
