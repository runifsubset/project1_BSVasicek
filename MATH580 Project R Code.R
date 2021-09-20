
## --------------------------------------------

## Course: MATH580 Stochastic Processes
## Title: Implementing Black-Scholes and Vasiciek
## Author: John Inston
## Date: November 2020
## Institution: Lancaster University

## --------------------------------------------

## Section 1: Libraries & Data ----------------
install.packages(ggplot2,reshape,tidyr,gridExtra)

library(ggplot2)
library(reshape)
library(tidyr)
library(gridExtra)

## Section 2: Black-Scholes Model -----------------

## 2.1 Price series plot

# Define initial values
S0 <-  1 
sigma <- sqrt(0.02)
rho = 0.03
c <- 1
t <- seq(0, 10, by=0.001)

# Define price
P <- S0*pnorm((log(S0/c)+(rho+(sigma^2)/2)*t)/(sigma*sqrt(t)))-(c*exp(-rho*t))*pnorm((log(S0/c)+(rho-(sigma^2)/2)*t)/(sigma*sqrt(t)))

# Produce plot for figure 1
fig1 <- data.frame(P, t)
plot1 <- ggplot(data = fig1, aes(t,P))+geom_line() + 
  ggtitle("Option Price over time.") +
  xlab("Time (t)") + ylab("Option Price (P)")

## 2.2 Investigation

# Set t=10
t10<-10

## 2.2.1 Varying volatility (sigma)

# Define varying sigma 
sigma1=seq(0,2, by=0.01)

# Calculate new price series
P1<-S0*pnorm((log(S0/c)+(rho+(sigma1^2)/2)*t10)/(sigma1*sqrt(t10)))-(c*exp(-rho*t10))*pnorm((log(S0/c)+(rho-(sigma1^2)/2)*t10)/(sigma1*sqrt(t10)))

# Produce plot for figure 2
fig2 <- data.frame(P1, sigma1)
plot2 <- ggplot(data = fig2, aes(sigma1,P1))+geom_line() + 
  ggtitle("Option price at time t=10 for various volatility levels.") +
  xlab("Volatility (sigma)") + ylab("Option Price (P) at t=10.")

## 2.2.2 Varying interest rate (rho)

# Define varying rho
rho2 <-seq(0, 0.8, by=0.005) 

# Calculate new price series
P2 <- S0*pnorm((log(S0/c)+(rho2+(sigma^2)/2)*t10)/(sigma*sqrt(t10)))-(c*exp(-rho2*t10))*pnorm((log(S0/c)+(rho2-(sigma^2)/2)*t10)/(sigma*sqrt(t10)))

# Produce plot for figure 3
fig3 <- data.frame(P2, rho2)
plot3 <- ggplot(data = fig3, aes(rho2,P2))+geom_line() + 
  ggtitle("Option price at time t=10 for various interest rates.") +
  xlab("Interest Rate (rho)") + ylab("Option Price (P) at t=10.")

## 2.2.3 Varying strike price (c)

# Define varying c
c3 <- seq(0,5, by=0.01)

# Calculate new price series
P3 <- S0*pnorm((log(S0/c3)+(rho+(sigma^2)/2)*t10)/(sigma*sqrt(t10)))-(c3*exp(-rho*t10))*pnorm((log(S0/c3)+(rho-(sigma^2)/2)*t10)/(sigma*sqrt(t10)))

# Produce plot for figure 4
fig4 <- data.frame(P3, c3)
plot4 <- ggplot(data = fig4, aes(c3,P3))+geom_line() + 
  ggtitle("Option price at time t=10 for various strike prices.") +
  xlab("Strike Price (c)") + ylab("Option price (P) at t=10.")

# Combining figures 2, 3 and 4
grid.arrange(plot2, plot3, plot4, ncol=3)



## Section 3: OU Process -----------------

## Define the OU function

rOU <- function(n,N,Delta,theta,sigma){ 
  times<-(0:n)*Delta ##vector of t_0,t_1,..,t_n 
  X<-matrix(0,nrow=N,ncol=n+1)
  for(i in 1:n){
    x<-X[,i]#current value
    m<-x*exp(-theta*Delta) #mean of new value 
    v<-sigma^2*(1-exp(-2*theta*Delta))/(2*theta) ##variance of new value 
    X[,i+1]=rnorm(N,m,sqrt(v)) ##simulate new value
  }
  return(list(X=X,times=times)) 
}

# Define initial conditions 

n=10000
N=10
Delta=10/n
theta=0.5
sigma=0.02

# Calculate 10 realisations of the OU process using our function
OU=rOU(n, N, Delta, theta, sigma)

# Produce plot for Figure 5
tX = data.frame(t(OU$X))
fig5 = data.frame(x=seq_along(tX[,1]), tX)
fig5 = melt(fig5, id.vars = "x")
cols = 2:11
plot5 = ggplot(fig5, aes(x = x, y = value, color = variable)) +
  geom_line()+ggtitle("OU process realisations over time.") +
  xlab("Time (t)") + ylab("OU process (X)")+
  guides(color = guide_legend(title = "Simulation"))+
  scale_color_manual(values = cols)


# Define transformation conditions
R0=0.1
mu=0.05

# Transform OU process

R=matrix(0, nrow=N, ncol=n+1)
for(i in 1:N){
  R[i,]=exp(-theta*OU$times)*R0+(1-exp(-theta*OU$times))*mu+OU$X[i,]
}

# Produce plot for Figure 6
tR = data.frame(t(R))
fig6 <- data.frame(x = seq_along(tR[, 1]), tR)
fig6 <- melt(fig6, id.vars = "x")
cols = 2:11
plot6 = ggplot(fig6, aes(x = x, y = value, color = variable)) +
  geom_line()+ggtitle("Spot rate simulations over time.") +
  xlab("Time (s)") + ylab("Spot Rate (R_s)")+
  guides(color = guide_legend(title = "Simulation"))+
  scale_color_discrete(labels = paste("S", 1:10))+
  scale_color_manual(values = cols)
  
# Combining figures 5 and 6
grid.arrange(plot5, plot6, ncol=2)

## Plot the Expectation and Variance

expRt=exp(-theta*OU$times)*R0+(1-exp(-theta*OU$times))*mu

fig7 = data.frame(expRt, OU$times)
plot7 = ggplot(data = fig7, aes(expRt,OU.times))+geom_line() + 
  ggtitle("OU process expected value over time.") +
  xlab("Time") + ylab("Expectated Value")

varRt=(sigma^2/(2*theta))*(1-exp(-2*theta*OU$times))

fig8 = data.frame(varRt, OU$times)
plot8 = ggplot(data = fig8, aes(varRt,OU.times))+geom_line() + 
  ggtitle("OU process variance over time.") +
  xlab("Time") + ylab("Variance")

# Combine figures 7 and 8
grid.arrange(plot7, plot8, ncol=2)


## Section 4: Vasicek Model -----------------

# Matrix to store Q
Q=matrix(0,nrow=N,ncol=n+1) 

# Loop to generate values of Q
for(i in 1:N){
  Q[i,]=exp(-Delta*cumsum(R[i,]))
}

# Plot bond price over time
tQ = data.frame(t(Q))
fig9 <- data.frame(x = seq_along(tQ[,1]),tQ)
fig9 <- melt(fig9, id.vars = "x")
cols = 2:11
plot9 = ggplot(fig9, aes(x = x, y = value, color = variable)) +
  geom_line()+ggtitle("Bond price Qt over time t.") +
  xlab("Time (t)") + ylab("Bond Price (Qt)")+
  guides(color = guide_legend(title = "Simulation"))+
  scale_color_manual(values = cols)


## Distribution of Qt

#Simulate 1000 realisations of Rt at t=10 (n=1000, delta=10/n so t=n*delta=10)

N1=1000
OU1=rOU(n, N1, Delta, theta, sigma)

R1=matrix(0, nrow=N1, ncol=n+1)

for(i in 1:N1){
  R1[i,]=exp(-theta*OU1$times)*R0+(1-exp(-theta*OU1$times))*mu+OU1$X[i,]
}

Q=rep(0, N1)
for(i in 1:N1){
  Q[i]=exp(-Delta*sum(R1[i, 2:(n+1)]))
}

hist(log(Q))
qqnorm(log(Q))
qqline(log(Q))


