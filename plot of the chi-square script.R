
setwd("C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\generalized chi-square")
# ==============================================================================
# R Script: Generalized Chi-Squared Power Distribution (Fixed Symbol Rendering)
# Manuscript Reference: Obulezi et al.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Distribution Functions
# ------------------------------------------------------------------------------

# Probability Density Function
dgchi_power <- function(q, n, alpha, sigma) {
  pdf_val <- numeric(length(q))
  pos_idx <- q > 0
  q_pos   <- q[pos_idx]
  
  const <- (2^(1 - n/2)) / (gamma(n/2) * (1 + alpha) * sigma^(n / (1 + alpha)))
  pdf_val[pos_idx] <- const * (q_pos^(n / (1 + alpha) - 1)) * 
    exp(-0.5 * (q_pos / sigma)^(2 / (1 + alpha)))
  
  return(pdf_val)
}

# Survival Function
sgchi_power <- function(q, n, alpha, sigma) {
  sf_val  <- numeric(length(q))
  pos_idx <- q > 0
  q_pos   <- q[pos_idx]
  
  z <- 0.5 * (q_pos / sigma)^(2 / (1 + alpha))
  sf_val[pos_idx]  <- pgamma(z, shape = n/2, lower.tail = FALSE)
  sf_val[!pos_idx] <- 1
  
  return(sf_val)
}

# Hazard Rate Function
hgchi_power <- function(q, n, alpha, sigma) {
  pdf_vals <- dgchi_power(q, n, alpha, sigma)
  sf_vals  <- sgchi_power(q, n, alpha, sigma)
  
  hrf_vals <- pdf_vals / sf_vals
  hrf_vals[sf_vals <= 1e-14] <- NA
  return(hrf_vals)
}

# ------------------------------------------------------------------------------
# 2. Synchronized Parameter Configurations
# ------------------------------------------------------------------------------
params <- list(
  list(n = 1, alpha =  1.5, sigma = 1.0, color = "#E41A1C", lty = 1),
  list(n = 1, alpha = -0.2, sigma = 1.0, color = "#377EB8", lty = 2),
  list(n = 3, alpha =  0.75, sigma = 0.25, color = "#4DAF4A", lty = 3),
  list(n = 2, alpha =  1.0, sigma = 1.0, color = "#984EA3", lty = 4),
  list(n = 2, alpha =  0.0, sigma = 1.5, color = "#FF7F00", lty = 5),
  list(n = 4, alpha =  0.0, sigma = 1.0, color = "#A65628", lty = 6)
)

# Robust plot math expression builder using bquote
legend_exprs <- as.expression(lapply(params, function(p) {
  bquote(n == .(p$n) * "," ~ alpha == .(p$alpha) * "," ~ sigma == .(p$sigma))
}))

q_seq <- seq(0.01, 4.5, length.out = 1000)

# ------------------------------------------------------------------------------
# 3. Plot 1: PDF Window
# ------------------------------------------------------------------------------
windows(width = 9.5, height = 7)

pdf_matrix <- sapply(params, function(p) dgchi_power(q_seq, p$n, p$alpha, p$sigma))

plot(q_seq, pdf_matrix[, 1], type = "l", col = params[[1]]$color, lty = params[[1]]$lty,
     lwd = 2.2, xlim = c(0, 4.5), ylim = c(0, 1.2),
     xlab = expression(q), ylab = expression(f[Q](q)),
     main = expression(bold("Gen-" * chi^2 * "(n," ~ alpha * "," ~ sigma * ")")),
     cex.lab = 1.2, cex.axis = 1.1)
grid(col = "gray85", lty = "dotted")

for (i in 2:6) {
  lines(q_seq, pdf_matrix[, i], col = params[[i]]$color, lty = params[[i]]$lty, lwd = 2.2)
}

legend("topright", 
       legend = legend_exprs,
       col    = sapply(params, function(p) p$color),
       lty    = sapply(params, function(p) p$lty),
       lwd    = 2.2, bty = "o", bg = "white", box.col = "gray70", cex = 0.85)

# ------------------------------------------------------------------------------
# 4. Plot 2: HRF Window
# ------------------------------------------------------------------------------
windows(width = 9.5, height = 7)

hrf_matrix <- sapply(params, function(p) hgchi_power(q_seq, p$n, p$alpha, p$sigma))

plot(q_seq, hrf_matrix[, 1], type = "l", col = params[[1]]$color, lty = params[[1]]$lty,
     lwd = 2.2, xlim = c(0, 4.5), ylim = c(0, 3),
     xlab = expression(q), ylab = expression(h[Q](q)),
     main = expression(bold(" Gen-" * chi^2 * "(n," ~ alpha * "," ~ sigma * ")")),
     cex.lab = 1.2, cex.axis = 1.1)
grid(col = "gray85", lty = "dotted")

for (i in 2:6) {
  lines(q_seq, hrf_matrix[, i], col = params[[i]]$color, lty = params[[i]]$lty, lwd = 2.2)
}

legend("topleft", 
       legend = legend_exprs,
       col    = sapply(params, function(p) p$color),
       lty    = sapply(params, function(p) p$lty),
       lwd    = 2.2, bty = "o", bg = "white", box.col = "gray70", cex = 0.85)