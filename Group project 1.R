S0 <- 150        
r <- 0.04       
sigma <- 0.3    
T <- 1           
m <- 12 
dt <- T / m      
K <- 150         
n_sim <- 10000   

payoffs <- numeric(n_sim) 

for (j in 1:n_sim) {
  Z <- rnorm(m, mean=0, sd=1)  
  S_Ti <- numeric(m)
  S_Ti[1] <- S0
  

for (i in 2:m) {
    S_Ti[i] <- S_Ti[i-1] * exp((r * (sigma^2) / 2) * dt + sigma * sqrt(dt) * Z[i-1])
}
A_T <- prod(S_Ti)^(1/m)
payoffs[j] <- max(A_T - K, 0)
}

V_t <- mean(payoffs) * exp(-r * T)

