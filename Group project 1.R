#1.1
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

#1.2
s <- seq(0, T, length.out = m)
d1_values <- numeric(m)
d2_values <- numeric(m)

for (i in 1:m) {
  d1_values[i] <- (log(S_Ti[i] / K) + (r + (sigma^2) / 2) * (T - s[i])) / (sigma * sqrt(T - s[i]))
  d2_values[i] <- d1_values[i] - sigma * sqrt(T - s[i])
}

v_call_values <- numeric(m)

for (i in 1:m) {
  v_call_values[i] <- max(S_Ti[i] * pnorm(d1_values[i]) - K * exp(-r * (T - s[i])) * pnorm(d2_values[i]), 0)
}
v_call_values

#1.3
K_star <- 225 
for (j in 1:n_sim) {
  Z_star <- rnorm(m, mean=0, sd=1)  
  S_Ti_star <- numeric(m)
  S_Ti_Star[1] <- S0
  barrier_hit <- FALSE   #barrier not passed
    for (i in 2:m) {
    S_Ti_star[i] <- S_Ti_star[i - 1] * exp((r - (sigma^2) / 2) * dt + sigma * sqrt(dt) * Z_star[i - 1])
    if (S_Ti_star[i] > K_star) {
      barrier_hit <- TRUE  # barrier passed
    }
  }
  
A_T_star <- prod(S_Ti_star)^(1/m)
payoffs_up_in <- numeric(n_sim)
payoffs_up_out <- numeric(n_sim)
  if (barrier_hit) {
    payoffs_up_in[j] <- max(A_T - K, 0)
  } else {
    payoffs_up_in[j] <- 0
  }  #barrier up-in
  
  if (!barrier_hit) {
    payoffs_up_out[j] <- max(A_T - K, 0)
  } else {
    payoffs_up_out[j] <- 0
  }
} #barrier up-out

V_t_up_in <- mean(payoffs_up_in) * exp(-r * T)
V_t_up_out <- mean(payoffs_up_out) * exp(-r * T)
print(V_t_up_in)
print(V_t_up_out)
