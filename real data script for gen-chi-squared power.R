# ==============================================================================
# Real Data Application of Generalized Chi-Square Power Distribution
# Dataset: Ball Bearing Fatigue Life Data (Lieblein & Zelen, 1956)
# Competing Models: Generalized Chi-Square Power (Proposed), Weibull, Chi
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package Dependencies & Data Input
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "gridExtra", "fitdistrplus", "stats", "goftest")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Real Dataset: 23 Deep-Groove Ball Bearing Fatigue Life Times (Millions of Revolutions)
ball_bearing_data <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 
                       51.96, 54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 
                       93.12, 98.64, 105.84, 127.92, 128.04, 173.40, 173.40)

x <- sort(ball_bearing_data)
n_obs <- length(x)

# ------------------------------------------------------------------------------
# 2. Density, Cumulative, and Quantile Functions for Gen-Chi^2 Power
# ------------------------------------------------------------------------------

# PDF: Gen-Chi^2 Power (n, alpha, sigma)
dgchipower <- function(q, n, alpha, sigma, log = FALSE) {
  if (n <= 0 || alpha <= -1 || sigma <= 0) return(rep(0, length(q)))
  log_c1 <- (1 - n/2)*log(2) - lgamma(n/2) - log(1 + alpha) - (n/(1 + alpha))*log(sigma)
  log_pdf <- log_c1 + (n/(1 + alpha) - 1)*log(q) - 0.5 * (q/sigma)^(2/(1 + alpha))
  if (log) return(log_pdf)
  return(exp(log_pdf))
}

# CDF: Gen-Chi^2 Power (n, alpha, sigma)
pgchipower <- function(q, n, alpha, sigma) {
  if (n <= 0 || alpha <= -1 || sigma <= 0) return(rep(0, length(q)))
  z <- 0.5 * (q / sigma)^(2 / (1 + alpha))
  return(pgamma(z, shape = n/2, scale = 1))
}

# Quantile Function: Gen-Chi^2 Power (n, alpha, sigma)
qgchipower <- function(p, n, alpha, sigma) {
  z <- qgamma(p, shape = n/2, scale = 1)
  q <- sigma * (2 * z)^((1 + alpha) / 2)
  return(q)
}

# ------------------------------------------------------------------------------
# 3. Model Fitting & Standard Errors via Hessian
# ------------------------------------------------------------------------------

# 3.1 Proposed Gen-Chi^2 Power Model
nll_gchipower <- function(params, data) {
  n <- params[1]; alpha <- params[2]; sigma <- params[3]
  if (n <= 0.01 || alpha <= -0.99 || sigma <= 0.001) return(1e10)
  log_lik <- sum(dgchipower(data, n, alpha, sigma, log = TRUE))
  if (is.na(log_lik) || is.nan(log_lik)) return(1e10)
  return(-log_lik)
}

fit_genchi <- optim(par = c(n = 3.0, alpha = 0.5, sigma = 10.0),
                    fn = nll_gchipower, data = x, method = "Nelder-Mead", hessian = TRUE)
par_genchi <- fit_genchi$par
se_genchi  <- sqrt(diag(solve(fit_genchi$hessian)))

# 3.2 Weibull Distribution
fit_weibull <- fitdist(x, "weibull")
par_weibull <- fit_weibull$estimate
se_weibull  <- fit_weibull$sd

# 3.3 Standard Chi Distribution (Scale-adjusted)
nll_chi <- function(params, data) {
  df <- params[1]; scale <- params[2]
  if (df <= 0.01 || scale <= 0.001) return(1e10)
  log_lik <- sum(log(2^(1-df/2) / (gamma(df/2) * scale^df)) + (df-1)*log(data) - 0.5*(data/scale)^2)
  return(-log_lik)
}
fit_chi <- optim(par = c(df = 2.0, scale = 50.0), fn = nll_chi, data = x, hessian = TRUE)
par_chi <- fit_chi$par
se_chi  <- sqrt(diag(solve(fit_chi$hessian)))

# ------------------------------------------------------------------------------
# 4. Computation of Goodness-of-Fit Metrics & Export CSV Table
# ------------------------------------------------------------------------------

compute_all_gof <- function(nll, k, sample_size, data, cdf_closure, mles, ses) {
  loglik <- -nll
  aic    <- 2 * k - 2 * loglik
  caic   <- aic + (2 * k^2 + 2 * k) / (sample_size - k - 1)
  bic    <- k * log(sample_size) - 2 * loglik
  hqic   <- 2 * k * log(log(sample_size)) - 2 * loglik
  
  # EDF Goodness-of-Fit Statistics using closures
  ks_res  <- ks.test(data, cdf_closure)
  ks_stat <- as.numeric(ks_res$statistic)
  p_val   <- as.numeric(ks_res$p.value)
  
  ad_stat  <- as.numeric(ad.test(data, null = cdf_closure)$statistic)
  cvm_stat <- as.numeric(cvm.test(data, null = cdf_closure)$statistic)
  
  # Format MLEs (Standard Errors)
  mle_se_str <- paste(sapply(1:length(mles), function(i) {
    sprintf("%.4f (%.4f)", mles[i], ses[i])
  }), collapse = "; ")
  
  return(data.frame(
    LogLik  = round(loglik, 4),
    AIC     = round(aic, 4),
    CAIC    = round(caic, 4),
    BIC     = round(bic, 4),
    HQIC    = round(hqic, 4),
    A       = round(ad_stat, 4),
    W       = round(cvm_stat, 4),
    KS      = round(ks_stat, 4),
    p_value = round(p_val, 4),
    MLEs_SE = mle_se_str
  ))
}

# Define closures for each model CDF
cdf_genchi_closure  <- function(v) pgchipower(v, n = par_genchi[1], alpha = par_genchi[2], sigma = par_genchi[3])
cdf_weibull_closure <- function(v) pweibull(v, shape = par_weibull[1], scale = par_weibull[2])
cdf_chi_closure     <- function(v) pgamma(0.5 * (v / par_chi[2])^2, shape = par_chi[1] / 2)

# Compile Table
gof_csv_table <- data.frame(
  Model = c("Gen-Chi^2 Power (Proposed)", "Weibull", "Chi"),
  rbind(
    compute_all_gof(fit_genchi$value, 3, n_obs, x, cdf_genchi_closure, par_genchi, se_genchi),
    compute_all_gof(-fit_weibull$loglik, 2, n_obs, x, cdf_weibull_closure, par_weibull, se_weibull),
    compute_all_gof(fit_chi$value, 2, n_obs, x, cdf_chi_closure, par_chi, se_chi)
  )
)

# Export to CSV
write.csv(gof_csv_table, "Goodness_of_Fit_Comparison.csv", row.names = FALSE)
cat(" Goodness-of-Fit Summary CSV exported successfully to: Goodness_of_Fit_Comparison.csv\n\n")
print(gof_csv_table)

# ------------------------------------------------------------------------------
# 5. Diagnostic Plots (Histogram with Legend Included)
# ------------------------------------------------------------------------------

theme_set(theme_bw(base_size = 12))

# --- Plot 1: Boxplot Superimposed on Violin Plot ---
df_data <- data.frame(Value = x)
p1 <- ggplot(df_data, aes(x = "Ball Bearings", y = Value)) +
  geom_violin(fill = "lightskyblue", alpha = 0.5, color = "dodgerblue4", linewidth = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.color = "red", outlier.size = 2) +
  geom_jitter(width = 0.05, color = "navy", alpha = 0.7, size = 2) +
  labs(title = "A: Data Distribution (Boxplot & Violin)", x = "", y = "Fatigue Life (Millions of Revolutions)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- Plot 2: Total Time on Test (TTT) Plot ---
i <- 1:n_obs
sum_x <- sum(x)
r_stat <- sapply(i, function(k) (sum(x[1:k]) + (n_obs - k) * x[k]) / sum_x)
df_ttt <- data.frame(i_n = i / n_obs, TTT = r_stat)

p2 <- ggplot(df_ttt, aes(x = i_n, y = TTT)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "darkblue", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "B: Total Time on Test (TTT) Plot", x = "i / n", y = "T(i/n)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- Plot 3: Empirical Density & Model Fits WITH HISTOGRAM LEGEND ---
x_grid <- seq(min(x) - 5, max(x) + 15, length.out = 300)
df_grid <- data.frame(
  x = x_grid,
  GenChi = dgchipower(x_grid, par_genchi[1], par_genchi[2], par_genchi[3]),
  Weibull = dweibull(x_grid, shape = par_weibull[1], scale = par_weibull[2]),
  Chi = dgchipower(x_grid, par_chi[1], 0, par_chi[2])
)

p3 <- ggplot(df_data) +
  geom_histogram(aes(x = Value, y = ..density.., fill = "Empirical Data"), bins = 8, color = "black", alpha = 0.4) +
  geom_density(aes(x = Value, color = "Kernel Density"), linewidth = 0.8, linetype = "dotted") +
  geom_line(data = df_grid, aes(x = x, y = GenChi, color = "Gen-Chi^2 Power"), linewidth = 1.2) +
  geom_line(data = df_grid, aes(x = x, y = Weibull, color = "Weibull"), linewidth = 1.0, linetype = "dashed") +
  geom_line(data = df_grid, aes(x = x, y = Chi, color = "Chi"), linewidth = 1.0, linetype = "dotdash") +
  scale_fill_manual(name = "Histogram", values = c("Empirical Data" = "grey80")) +
  scale_color_manual(name = "Model Density", values = c("Kernel Density" = "black", "Gen-Chi^2 Power" = "red", 
                                                        "Weibull" = "blue", "Chi" = "darkorange")) +
  labs(title = "C: Fitted Densities Superimposed on Data Histogram", x = "Fatigue Life", y = "Density") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")

# --- Plot 4: Survival Function Fits ---
emp_surv <- 1 - (1:n_obs - 0.5) / n_obs
df_surv <- data.frame(
  x = x_grid,
  GenChi = 1 - pgchipower(x_grid, par_genchi[1], par_genchi[2], par_genchi[3]),
  Weibull = 1 - pweibull(x_grid, shape = par_weibull[1], scale = par_weibull[2]),
  Chi = 1 - pgamma(0.5 * (x_grid / par_chi[2])^2, shape = par_chi[1] / 2)
)

p4 <- ggplot() +
  geom_point(data = data.frame(x = x, Surv = emp_surv), aes(x = x, y = Surv), size = 2) +
  geom_line(data = df_surv, aes(x = x, y = GenChi, color = "Gen-Chi^2 Power"), linewidth = 1.2) +
  geom_line(data = df_surv, aes(x = x, y = Weibull, color = "Weibull"), linewidth = 1.0, linetype = "dashed") +
  geom_line(data = df_surv, aes(x = x, y = Chi, color = "Chi"), linewidth = 1.0, linetype = "dotdash") +
  scale_color_manual(name = "Models", values = c("Gen-Chi^2 Power" = "red", "Weibull" = "blue", "Chi" = "darkorange")) +
  labs(title = "D: Empirical vs Fitted Survival Functions", x = "Fatigue Life", y = "Survival Probability") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")

# --- Plot 5: Probability-Probability (P-P) Plot ---
p_emp <- (1:n_obs - 0.5) / n_obs
p_genchi <- pgchipower(x, par_genchi[1], par_genchi[2], par_genchi[3])

p5 <- ggplot(data.frame(Empirical = p_emp, Model = p_genchi), aes(x = Empirical, y = Model)) +
  geom_point(color = "red", size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  labs(title = "E: P-P Plot (Gen-Chi^2 Power)", x = "Empirical Cumulative Probability", y = "Theoretical Probability") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# --- Plot 6: Quantile-Quantile (Q-Q) Plot ---
q_genchi <- qgchipower(p_emp, par_genchi[1], par_genchi[2], par_genchi[3])

p6 <- ggplot(data.frame(Observed = x, Theoretical = q_genchi), aes(x = Theoretical, y = Observed)) +
  geom_point(color = "red", size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  labs(title = "F: Q-Q Plot (Gen-Chi^2 Power)", x = "Theoretical Quantiles", y = "Observed Quantiles") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ------------------------------------------------------------------------------
# 6. Display & Save Combined Diagnostic Panel
# ------------------------------------------------------------------------------
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)

ggsave("Gen_ChiSquare_Power_RealData_Diagnostics.png", 
       plot = marrangeGrob(list(p1, p2, p3, p4, p5, p6), nrow = 3, ncol = 2), 
       width = 14, height = 14, dpi = 300)