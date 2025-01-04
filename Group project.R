S0 <- 150   # Initial price of the stock
K = 150
sigma <- 0.3        # volatility
r <- 0.04        # risk-free rate
T <- 1             # TTM
t=0.5
tau= T-t
m <- 12            # intervals number
delta_t <- T / m   # lenght of each interval
W <- numeric(m + 1)
W[1] <- 0      # We assumed as initial condition of our brownian motion W[1]=0

Z <- rnorm(m, mean = 0, sd = 1)  
for (i in 2:(m + 1)) {
  W[i] <- W[i - 1] + sqrt(delta_t) * Z[i - 1]
}  
print(W)
# I put the square of delta t because Ti-Ti-1 is delta T=T/m and so I made the square

S_T <- S0 * exp((r - (sigma^2) / 2) * tau + sigma * sqrt(tau) * Z)

i_t <- ceiling(t / delta_t)  #approximation to the higher value of t

alpha_t <- function(t) {
  return(1)  #return(1) because we don't have a value so i choose 1
} 

sigma_A_t <- (sigma * sqrt(delta_t) / m) * (m - i_t + 1) *
  sqrt(
    alpha_t(t) + ((m - i_t) * (2 * (m - i_t) + 1)) / (6 * (m - i_t + 1))
  )
#compute the value of sigma A(t)* (T-t)

mu_a_t <- ((r - sigma^2 / 2) * delta_t / m) * (
  (m - i_t + 1) * alpha_t(t) + ((m - i_t) * (m - i_t + 1)) / 2
) - (r - sigma_A_t^2 / 2) * (T - t)
#compute the value of mu A(t)* (T-t)

Y_t <- prod(S_Ti[1:(i_t - 1)]) * S_Ti[i_t]
X_t <- Y_t * exp(mu_a_t)
print(X_t)
S_T_tilde <- S_T * exp(
  (r - (sigma_A_t^2) / 2) * tau + sigma_A_t * sqrt(tau) * Z
)
payoff <- pmax((S_T_tilde - K) / X_t, 0)

expected_payoff <- mean(payoff_normalized)
V_t <- X_t * expected_payoff * exp(-r * tau)
