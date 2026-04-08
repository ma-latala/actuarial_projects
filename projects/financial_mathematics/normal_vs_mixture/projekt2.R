
### read and transform stock price data
#######################################

data =read.csv('stock_data.csv', header=1, row.names=1)

n_stocks = length(data)

N = length(data[,1])

log_returns = log(data[-1, ]) - log(data[-nrow(data), ])


dgnorm = function(x, params){
  params = as.matrix(params)
  beta = params[3]
  alpha = params[2]
  mu = params[1]
  return(beta/2/alpha/gamma(1/beta)*exp(-(abs(x-mu)/alpha)^beta))
}

pgnorm = function(x, params) {
  params = as.matrix(params)
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  
  if (alpha <= 0 || beta <= 0) {
    stop("alpha and beta must be positive")
  }
  
  z <- (abs(x - mu) / alpha)^beta
  G <- pgamma(z, shape = 1 / beta, scale = 1)
  
  ifelse(x < mu,
         0.5 * (1 - G),
         0.5 * (1 + G))
}

dmixture = function(x, params){
  params = as.matrix(params)
  sigma = params[1]
  beta  = params[2]
  p     = params[3]
  return(p*dnorm(x, sd=sigma)+(1-p)*dnorm(x, sd=beta)) 
}

pmixture = function(x, params){
  params = as.matrix(params)
  sigma = params[1]
  beta  = params[2]
  p     = params[3]
  return(p*pnorm(x, sd=sigma)+(1-p)*pnorm(x, sd=beta)) 
}

#######################################################
# ML ESTIMATION
#######################################################


library(maxLik)

# metoda MLE dla rozkładu normalnego

mu_0_p_values = rep(0, n_stocks)
results_normal = data.frame(
  mu = numeric(),
  sigma = numeric(),
  SE_mu = numeric(),
  SE_sigma = numeric()
)



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
  
  # test mu = 0
  
  mu_0 = 0
  Z = (wynik$estimate[1]-mu_0)/SE_vector[1]
  p_value = 2 * pnorm(-abs(Z))
  mu_0_p_values[s] = p_value
}



### MLE dla mieszanki rozkładów normalnych
###########################################

mixture_p_p_vals = rep(0, n_stocks)

results_mixture = data.frame(
  sigma = numeric(),
  beta = numeric(),
  p = numeric(),
  SE_sigma = numeric(),
  SE_beta = numeric(),
  SE_p = numeric()
)


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
  
  ###
  # test sigma=beta
  
  
  loglik = function(params){
    sigma = params[1]
    beta  = params[1]
    p     = params[2]
    
    if(sigma <= 0 | beta <= 0 | p < 0 | p > 1) return(NA)
    
    return(sum(log(p/sqrt(2*pi)/sigma*exp(-x^2/2/sigma^2) + (1-p)/sqrt(2*pi)/beta*exp(-x^2/2/beta^2))))
  }
  
  
  grad = function(params){
    sigma = params[1]
    beta  = params[1]
    p     = params[2]
    
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
    
    return(c(d_sigma, d_p))
  }
  
  
  wynik_restricted = maxNR(
    loglik,
    grad = grad,
    start  = c(.1, .5)
  )
  
  
  lnL_unrestricted = wynik$maximum
  lnL_restricted = wynik_restricted$maximum
  LR = 2*(lnL_unrestricted - lnL_restricted)
  g = 2
  alpha=0.05
  qchisq(1-alpha, g)
  p_value = 1-pchisq(LR, g)
  
  mixture_p_p_vals[s] = p_value

}



# metoda MLE dla generalized normal 
############################################

results_generalized = data.frame(
  mu = numeric(),
  alpha = numeric(),
  beta = numeric(),
  SE_mu = numeric(),
  SE_alpha = numeric(),
  SE_beta = numeric()
)

generalized_beta_2_p_values = rep(0, n_stocks)

g_beta <- function(beta, x, mu) {
  if (beta <= 0) return(NA_real_)
  
  y <- abs(x - mu)
  
  # problematic when any y == 0 because log(0) appears
  if (any(y == 0)) return(NA_real_)
  
  S0 <- sum(y^beta)
  S1 <- sum((y^beta) * log(y))
  
  1 +
    digamma(1 / beta) / beta -
    S1 / S0 +
    log((beta / length(x)) * S0) / beta
}


gprime_beta <- function(beta, x, mu) {
  if (beta <= 0) return(NA_real_)
  
  y <- abs(x - mu)
  
  # problematic when any y == 0 because log(0) appears
  if (any(y == 0)) return(NA_real_)
  
  S0 <- sum(y^beta)
  L1 <- sum((y^beta) * log(y))
  L2 <- sum((y^beta) * (log(y)^2))
  
  -digamma(1 / beta) / beta^2 -
    trigamma(1 / beta) / beta^3 +
    1 / beta^2 -
    L2 / S0 +
    (L1^2) / (S0^2) +
    L1 / (beta * S0) -
    log((beta / length(x)) * S0) / beta^2
}

estimate_params = function(x){
  
  perc_diff = 100
  
  mu = mean(x)
  beta = mean(x)/ sd(x)
  
  while(perc_diff>0.0001){
    new_beta = beta - g_beta(beta, x, mu)/gprime_beta(beta, x, mu)
    perc_diff = abs(beta-new_beta)/beta
    beta = new_beta
  }
  
  target_mu = function(y) return(-sum((abs(x-y))^beta))
  wynik = maxNR(
    target_mu,
    start  = c(mean(x))
  )
  mu = wynik$estimate[1]
  
  alpha = (beta/length(x)*sum((abs(x-mu))^beta))^(1/beta)
  
  return(c(mu, alpha, beta))

}


for(s in 1:n_stocks){
  x = log_returns[,s]

  res = estimate_params(x)
  
  results_generalized[s, ] = c(res, c(NA, NA, NA))
  
  
  # test beta=2
  
    loglik_unrestricted = sum(log(dgnorm(x, res)))
  
    mu = mean(x)
    beta = 2
    
    target_mu = function(y){return(-sum((abs(x-y))^beta))}
    wynik = maxNR(target_mu, start  = c(mean(x)))

    mu = wynik$estimate[1]
    
    alpha = (beta/length(x)*sum((abs(x-mu))^beta))^(1/beta)
    
    loglik_restricted = sum(log(dgnorm(x, c(mu, alpha, beta))))
  
    LR = 2*(loglik_unrestricted - loglik_restricted)
    g = 1
    alpha=0.05
    qchisq(1-alpha, g)
    p_value = 1-pchisq(LR, g)
    
    generalized_beta_2_p_values[s] = p_value
}







### histogram plots + save to svg
####################################################

colors[2] = "#F0F086"

for(s in 1:n_stocks){
  svg(file=paste("plots/plot", s, ".svg", sep=""),
       width=6, height=5,pointsize = 10)
  
  x = log_returns[,s]
  
  curve_x = seq(min(x), max(x), length.out = 1000)
  
  hist(x, breaks = "Scott", freq = FALSE, col = colors[s], border = "white",
       main = paste("Ticker:",colnames(log_returns)[s]),
       xlab = "log return", ylab="Density", xlim=c(-0.15, 0.15), ylim = c(0, max(dgnorm(curve_x, results_generalized[s,]))))
  
  
  
  lines(curve_x, dmixture(curve_x, results_mixture[s, 1:3]), 
        col = "black", lwd = 2, lty = 2)
  
  
  
  lines(curve_x, dnorm(curve_x, mean = results_normal[s, 1],
                       sd   = results_normal[s, 2]), 
        col = "black", lwd = 2, lty = 1)
  
  
  lines(curve_x, dgnorm(curve_x, results_generalized[s,]),
        col = "black", lwd = 2, lty = 3)
  
  
  legend("topright",
         legend = c("MLE: Gaussian", "MLE: Gaussian mixture",  "MLE: generalized Gaussian"),
         col = c("black", "black", "black"),
         lty = c(1, 2, 3), lwd = 2)
  
  dev.off()
  
}



## -------------------------------------------------
## helper functions
## -------------------------------------------------

# normal
cdf_normal <- function(x, pars) {
  pars = as.matrix(pars)
  pnorm(x, mean = pars[1], sd = pars[2])
}

pdf_normal <- function(x, pars) {
  pars = as.matrix(pars)
  dnorm(x, mean = pars[1], sd = pars[2])
}

# mixture
# assumes results_mixture[s, 1:3] are exactly the parameters expected
# by your existing functions dmixture() and mixture_CDF()
cdf_mixture <- function(x, pars) {
  pmixture(x, pars)
}

pdf_mixture <- function(x, pars) {
  dmixture(x, pars)
}

# generalized normal
# replace pgnorm() below by your actual CDF function if needed
cdf_generalized <- function(x, pars) {
  pgnorm(x, pars)
}

pdf_generalized <- function(x, pars) {
  dgnorm(x, pars)
}

## -------------------------------------------------
## colors
## -------------------------------------------------

if (n_stocks <= 12) {
  colors <- RColorBrewer::brewer.pal(n_stocks, "Set3")
} else {
  colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_stocks)
}
colors[2] = "#F0F086"

## -------------------------------------------------
## plotting loop: one figure per stock
## -------------------------------------------------

old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

for (s in 1:n_stocks) {
  svg(file=paste("plots/tail_plot", s, ".svg", sep=""),
      width=8, height=3.5,pointsize = 12)
  
  r <- log_returns[, s]
  r <- r[is.finite(r)]
  
  ecdf_func <- ecdf(r)
  
  # grids
  x_left  <- seq(-0.20, 0.00, by = 0.001)
  x_mid   <- seq(-0.025, 0.03, by = 0.001)
  x_right <- seq( 0.00, 0.20, by = 0.001)
  
  # empirical objects
  emp_left  <- ecdf_func(x_left)
  emp_mid  <- ecdf_func(x_mid)
  emp_right <- 1 - ecdf_func(x_right)
  
  # fitted CDFs
  norm_left  <- cdf_normal(x_left,  results_normal[s, ])
  mix_left   <- cdf_mixture(x_left, results_mixture[s, 1:3])
  gnorm_left <- cdf_generalized(x_left, results_generalized[s, ])
  
  norm_right  <- 1 - cdf_normal(x_right,  results_normal[s, ])
  mix_right   <- 1 - cdf_mixture(x_right, results_mixture[s, 1:3])
  gnorm_right <- 1 - cdf_generalized(x_right, results_generalized[s, ])
  
  # fitted PDFs
  norm_mid  <- cdf_normal(x_mid,  results_normal[s, ])
  mix_mid   <- cdf_mixture(x_mid, results_mixture[s, 1:3])
  gnorm_mid <- cdf_generalized(x_mid, results_generalized[s, ])
  
  ymax_mid <- max(
    hist(r, breaks = "Scott", plot = FALSE)$density,
    norm_mid, mix_mid, gnorm_mid,
    na.rm = TRUE
  )
  
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1), oma = c(0, 0, 2, 0))
  
  ## ---------------------------------------------
  ## 1) LEFT PANEL: left tail
  ## ---------------------------------------------
  plot(
    NA, NA,
    xlim = range(x_left),
    ylim = c(1e-6, 1),
    log  = "y",
    xlab = "Log return",
    ylab = expression(hat(F)(x)),
    main = "Left tail"
  )
  
  grid(col = "grey85", lty = "dotted")
  
  # empirical in stock color
  lines(x_left,  pmax(emp_left,  1e-12), col = colors[s], lwd = 2.2)
  
  # fitted in black with your line styles
  lines(x_left, pmax(norm_left,  1e-12), col = "black", lwd = 2, lty = 1)
  lines(x_left, pmax(mix_left,   1e-12), col = "black", lwd = 2, lty = 2)
  lines(x_left, pmax(gnorm_left, 1e-12), col = "black", lwd = 2, lty = 3)
  
  legend(
    "topleft",
    legend = c("Empirical", "Normal", "Gaussian mixture", "Generalized normal"),
    col    = c(colors[s], "black", "black", "black"),
    lwd    = 2,
    lty    = c(1, 1, 2, 3),
    bty    = "n"
  )
  
  ## ---------------------------------------------
  ## 2) MIDDLE PANEL: PDF around 0
  ## ---------------------------------------------
  plot(
    NA, NA,
    xlim = range(x_mid),
    ylim = c(1e-1, 1),
    log  = "y",
    xlab = "Log return",
    ylab = expression(hat(F)(x)),
    main = "Around 0"
  )
  
  grid(col = "grey85", lty = "dotted")
  
  # empirical in stock color
  lines(x_mid,  pmax(emp_mid,  1e-12), col = colors[s], lwd = 2.2)
  
  # fitted in black with your line styles
  lines(x_mid, pmax(norm_mid,  1e-12), col = "black", lwd = 2, lty = 1)
  lines(x_mid, pmax(mix_mid,   1e-12), col = "black", lwd = 2, lty = 2)
  lines(x_mid, pmax(gnorm_mid, 1e-12), col = "black", lwd = 2, lty = 3)
  
  ## ---------------------------------------------
  ## 3) RIGHT PANEL: right tail
  ## ---------------------------------------------
  plot(
    NA, NA,
    xlim = range(x_right),
    ylim = c(1e-6, 1),
    log  = "y",
    xlab = "Log return",
    ylab = expression(1 - hat(F)(x)),
    main = "Right tail"
  )
  
  grid(col = "grey85", lty = "dotted")
  
  # empirical in stock color
  lines(x_right, pmax(emp_right,  1e-12), col = colors[s], lwd = 2.2)
  
  # fitted in black
  lines(x_right, pmax(norm_right,  1e-12), col = "black", lwd = 2, lty = 1)
  lines(x_right, pmax(mix_right,   1e-12), col = "black", lwd = 2, lty = 2)
  lines(x_right, pmax(gnorm_right, 1e-12), col = "black", lwd = 2, lty = 3)
  
  mtext(
    paste("Ticker:", colnames(log_returns)[s]),
    outer = TRUE,
    cex = 1.2
  )
  dev.off()
}

