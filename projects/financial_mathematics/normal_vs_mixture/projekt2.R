
### read and transform stock price data
#######################################

data =read.csv('stock_data.csv', header=1, row.names=1)

n_stocks = length(data)

N = length(data[,1])

log_returns = log(data[-1, ]) - log(data[-nrow(data), ])

mixture_prob_density = function(x, params){
  params = as.matrix(params)
  sigma = params[1]
  beta  = params[2]
  p     = params[3]
  return((1 - p) * (1 - p) * dnorm(x, sd=beta) + p * dnorm(x, sd=sigma))
}

mixture_CDF = function(x, params){
  params = as.matrix(params)
  sigma = params[1]
  beta  = params[2]
  p     = params[3]
  return((1 - p) * pnorm(x, sd=beta) + p * pnorm(x, sd=sigma))
}

### method of moments visualization
####################################

# colors
if (n_stocks <= 12) {
  colors <- RColorBrewer::brewer.pal(n_stocks, "Set3")
} else {
  colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_stocks)
}

# save old graphics settings
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

# side-by-side layout
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))


### ---------------------------
### RIGHT PANEL
### ---------------------------

x <- seq(-0.2, 0.05, by = 0.001)

plot(
  NA, NA,
  xlim = range(x),
  ylim = c(1e-4, 1),
  log = "y",
  xlab = "Log return",
  ylab = expression(hat(F)(x)),
  main = "Left tail"
)

grid(col = "grey85", lty = "dotted")

for (s in 1:n_stocks) {
  r <- log_returns[, s]
  r <- r[is.finite(r)]
  
  ecdf_func <- ecdf(r)
  emp_cdf  <- ecdf_func(x)
  
  mu  <- mean(r)
  sig <- sd(r)
  norm_cdf <- pnorm(x, mean = results_normal[s, 1],
                    sd   = results_normal[s, 2])
  
  # empirical tail
  lines(x, pmax(emp_cdf, 1e-8), col = colors[s], lwd = 2)
  
  # fitted normal tail
  lines(x, pmax(norm_cdf, 1e-8), col = colors[s], lwd = 2, lty = 2)
}

legend(
  "topleft",
  legend = c("Empirical CDF", "Fitted normal CDF"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)


### ---------------------------
### right PANEL
### ---------------------------
x <- seq(-0.05, 0.2, by = 0.001)

plot(
  NA, NA,
  xlim = range(x),
  ylim = c(1e-4, 1),
  log = "y",
  xlab = "Log return",
  ylab = expression(1 - hat(F)(x)),
  main = "Right tail"
)

grid(col = "grey85", lty = "dotted")

for (s in 1:n_stocks) {
  r <- log_returns[, s]
  r <- r[is.finite(r)]
  
  ecdf_func <- ecdf(r)
  emp_tail  <- 1 - ecdf_func(x)
  
  mu  <- mean(r)
  sig <- sd(r)
  norm_tail <- 1 - pnorm(x, mean = results_normal[s, 1],
                         sd   = results_normal[s, 2])
  
  # empirical tail
  lines(x, pmax(emp_tail, 1e-8), col = colors[s], lwd = 2)
  
  # fitted normal tail
  lines(x, pmax(norm_tail, 1e-8), col = colors[s], lwd = 2, lty = 2)
}

legend(
  "topright",
  legend = c("Empirical tail", "Fitted normal tail"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

mtext("Empirical vs Gaussian CDF", outer = TRUE, cex = 1.1)


###
####################################################


### method of moments visualization
####################################

# colors
if (n_stocks <= 12) {
  colors <- RColorBrewer::brewer.pal(n_stocks, "Set3")
} else {
  colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_stocks)
}

# save old graphics settings
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

# side-by-side layout
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))


### ---------------------------
### RIGHT PANEL
### ---------------------------

x <- seq(-0.2, 0.05, by = 0.001)

plot(
  NA, NA,
  xlim = range(x),
  ylim = c(1e-4, 1),
  log = "y",
  xlab = "Log return",
  ylab = expression(hat(F)(x)),
  main = "Left tail"
)

grid(col = "grey85", lty = "dotted")

for (s in 1:n_stocks) {
  r <- log_returns[, s]
  r <- r[is.finite(r)]
  
  ecdf_func <- ecdf(r)
  emp_cdf  <- ecdf_func(x)

  mix_cdf <- mixture_CDF(x, results_mixture[s, 1:3])
  
  # empirical tail
  lines(x, pmax(emp_cdf, 1e-8), col = colors[s], lwd = 2)
  
  # fitted mixture tail
  lines(x, pmax(mix_cdf, 1e-8), col = colors[s], lwd = 2, lty = 2)
}

legend(
  "topleft",
  legend = c("Empirical CDF", "Fitted normal CDF"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)


### ---------------------------
### right PANEL
### ---------------------------
x <- seq(-0.05, 0.2, by = 0.001)

plot(
  NA, NA,
  xlim = range(x),
  ylim = c(1e-4, 1),
  log = "y",
  xlab = "Log return",
  ylab = expression(1 - hat(F)(x)),
  main = "Right tail"
)

grid(col = "grey85", lty = "dotted")

for (s in 1:n_stocks) {
  r <- log_returns[, s]
  r <- r[is.finite(r)]
  
  ecdf_func <- ecdf(r)
  emp_tail  <- 1 - ecdf_func(x)
  
  mu  <- mean(r)
  sig <- sd(r)
  mix_tail <- 1 - mixture_CDF(x, results_mixture[s, 1:3])
  
  # empirical tail
  lines(x, pmax(emp_tail, 1e-8), col = colors[s], lwd = 2)
  
  # fitted normal tail
  lines(x, pmax(mix_tail, 1e-8), col = colors[s], lwd = 2, lty = 2)
}

legend(
  "topright",
  legend = c("Empirical tail", "Fitted normal tail"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

mtext("Empirical vs Mixed Gaussian CDF", outer = TRUE, cex = 1.1)


###
####################################################

for(s in 1:n_stocks){
  pdf(file=paste("plots/plot", s, ".pdf", sep=""),
      width=8, height=6,)
  
  x = log_returns[,s]
  
  hist(x, breaks = "Scott", freq = FALSE, col = colors[s], border = "white",
       main = paste("Ticker:",colnames(log_returns)[s]),
       xlab = "log return", ylab="Density", xlim=c(-0.2, 0.2))
  
  curve_x = seq(min(x), max(x), length.out = 1000)
  
  
  lines(curve_x, mixture_prob_density(curve_x, results_mixture[s, 1:3]), 
        col = "black", lwd = 2.5, lty = 2)
  
  
  
  lines(curve_x, dnorm(curve_x, mean = results_normal[s, 1],
                      sd   = results_normal[s, 2]), 
        col = "black", lwd = 2.5, lty = 3)
  
  
  
  legend("topright", 
         legend = c("MLE - gaussian mixture", "MLE - normal"),
         col = c(colors[s], colors[s]), 
         lty = c(2,3), lwd = 2)
  
  dev.off()

}






















library(maxLik)

results_normal = data.frame(
  mu = numeric(),
  sigma = numeric(),
  SE_mu = numeric(),
  SE_sigma = numeric()
)

results_mixture = data.frame(
  sigma = numeric(),
  beta = numeric(),
  p = numeric(),
  SE_sigma = numeric(),
  SE_beta = numeric(),
  SE_p = numeric()
)


# metoda MLE dla rozkładu normalnego

for(s in 1:n_stocks){
  x = log_returns[,s]
  
  loglik = function(params){
    sigma = params[2]
    mu = params[1]
    n = length(x)
    return(-n/2*log(2*pi)-n*log(sigma)-1/2/sigma/sigma*sum((x-mu)^2))
  }
  
  
  grad = function(params){
    mu = params[1]
    sigma = params[2]
    n = length(x)
    
    d_mu = sum(x - mu) / sigma^2
    d_sigma = -n / sigma + sum((x - mu)^2) / sigma^3
    
    c(d_mu, d_sigma)
  }
  
  hess = function(params){
    mu = params[1]
    sigma = params[2]
    n = length(x)
    
    h11 = -n / sigma^2
    h12 = -2 * sum(x - mu) / sigma^3
    h22 = n / sigma^2 - 3 * sum((x - mu)^2) / sigma^4
    
    matrix(c(h11, h12,
             h12, h22), nrow = 2)
  }
  
  wynik = maxNR(
    loglik,
    grad   = grad,
    hess   = hess,
    start  = c(mean(x), sd(x))
  )
  
  
  SE_vector = sqrt(diag(solve(-wynik$hessian)))
  
  results_normal[s, ] = list(
    wynik$estimate[1], wynik$estimate[2],
    SE_vector[1], SE_vector[2]
  )
}

# test mu = mu_sp i sigma= sd_sp

mu_0 = mu_sp
sigma_0 = sd_sp

lnL_unrestricted = wynik$maximum

lnL_restricted = loglik(c(mu_0, sigma_0))

LR = 2*(lnL_unrestricted - lnL_restricted)
g = 2
alpha=0.05

qchisq(1-alpha, g)

p_value = 1-pchisq(LR, g)
p_value

(p_value)


# test mu = 0

mu_0 = 0 
Z = (wynik$estimate[1]-mu_0)/SE_vector[1]
p_value = 2 * pnorm(-abs(Z))

(p_value)


### MLE dla mieszanki rozkładów normalnych
###########################################


for(s in 1:n_stocks){
  x = log_returns[,s]
  loglik = function(params){
    sigma = params[1]
    beta  = params[2]
    p     = params[3]
    
    if(sigma <= 0 | beta <= 0 | p < 0 | p > 1) return(NA)
    
    return(sum(log(p/sqrt(2*pi)/sigma*exp(-x^2/2/sigma^2) + (1-p)/sqrt(2*pi)/beta*exp(-x^2/2/beta^2))))
  }
  
  
  grad = function(params){
    sigma = params[1]
    beta  = params[2]
    p     = params[3]
    
    num_sigma = 
      (p * x^2 * exp(-x^2 / (2 * sigma^2))) / (sqrt(2*pi) * sigma^4) - 
      (p * exp(-x^2 / (2 * sigma^2))) / (sqrt(2*pi) * sigma^2)
    
    num_beta = 
      ((1 - p) * x^2 * exp(-x^2 / (2 * beta^2))) / (sqrt(2*pi) * beta^4) - 
      ((1 - p) * exp(-x^2 / (2 * beta^2))) / (sqrt(2*pi) * beta^2)
    
    num_p = 
      exp(-x^2 / (2 * sigma^2)) / (sqrt(2*pi) * sigma) - 
      exp(-x^2 / (2 * beta^2)) / (sqrt(2*pi) * beta)
    
    denom = 
      ((1 - p) * exp(-x^2 / (2 * beta^2))) / (sqrt(2*pi) * beta) + 
      (p * exp(-x^2 / (2 * sigma^2))) / (sqrt(2*pi) * sigma)
    
    d_sigma = sum(num_sigma / denom)
    d_beta  = sum(num_beta / denom)
    d_p     = sum(num_p / denom)
    
    return(c(d_sigma, d_beta, d_p))
  }
  
  
  wynik = maxNR(
    loglik,
    grad = grad,
    start  = c(sd(x), sd(x) , 0.95)
  )
  
  summary(wynik)
  
  SE_vector = sqrt(diag(solve(-wynik$hessian)))
  
  results_mixture[s, ] = list(
    wynik$estimate[1], wynik$estimate[2],wynik$estimate[3],
    SE_vector[1], SE_vector[2], SE_vector[3]
  )
}

# test sigma = beta = sqrt(p*sigma^2 + (1-p)*beta^2)

sigma_0 = sd_sp
beta_0 = sd_sp

lnL_unrestricted = wynik$maximum

lnL_restricted = loglik(c(sigma_0, beta_0, wynik$estimate[3]))

LR = 2*(lnL_unrestricted - lnL_restricted)
g = 2
alpha=0.05

qchisq(1-alpha, g)

p_value = 1-pchisq(LR, g)

(p_value)


# test p=0

p_0 = 0 
Z = (wynik$estimate[3]-p_0)/SE_vector[3]
p_value = 2 * pnorm(-abs(Z))

(p_value)

