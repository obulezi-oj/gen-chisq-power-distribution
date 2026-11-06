# ==============================================================================
# Monte Carlo & MCMC Simulation Script for Gen-Chi^2(n, alpha, sigma)
# Scenario 2 Parameter Configuration: n = 3.0, alpha = 0.75, sigma = 0.80
# ==============================================================================

suppressPackageStartupMessages({
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("gridExtra")) install.packages("gridExtra")
  library(ggplot2)
  library(gridExtra)
})

set.seed(2026)

# ------------------------------------------------------------------------------
# 1. Random Variate Generation & Log-Likelihood Functions
# ------------------------------------------------------------------------------
rgchi_power <- function(m, n, alpha, sigma) {
  r_sq <- rchisq(m, df = n)
  r <- sqrt(r_sq)
  q <- sigma * (r^(1 + alpha))
  return(q)
}

log_lik_gchi <- function(params, q) {
  n     <- params[1]
  alpha <- params[2]
  sigma <- params[3]
  
  if (n <= 0.01 || alpha <= -0.99 || sigma <= 0.001) return(-1e10)
  
  m <- length(q)
  term1 <- m * (1 - n/2) * log(2) - m * lgamma(n/2) - m * log(1 + alpha) - (m * n / (1 + alpha)) * log(sigma)
  term2 <- (n / (1 + alpha) - 1) * sum(log(q))
  term3 <- -0.5 * (1 / (sigma^(2 / (1 + alpha)))) * sum(q^(2 / (1 + alpha)))
  
  val <- term1 + term2 + term3
  if (is.na(val) || is.nan(val)) return(-1e10)
  return(val)
}

# ------------------------------------------------------------------------------
# 2. Estimation Routines (MLE & MCMC Bayesian)
# ------------------------------------------------------------------------------
fit_mle <- function(q, init_params = c(n = 3, alpha = 0.75, sigma = 0.8)) {
  neg_log_lik <- function(p) -log_lik_gchi(p, q)
  
  fit <- tryCatch({
    optim(par = init_params, fn = neg_log_lik, method = "L-BFGS-B",
          lower = c(0.1, -0.95, 0.01), upper = c(50, 10, 50), hessian = TRUE)
  }, error = function(e) { NULL })
  
  if (is.null(fit) || fit$convergence != 0) {
    fit <- tryCatch({
      optim(par = init_params, fn = neg_log_lik, method = "Nelder-Mead", hessian = TRUE)
    }, error = function(e) { NULL })
  }
  
  if (is.null(fit)) return(list(est = rep(NA, 3), se = rep(NA, 3)))
  
  est <- fit$par
  se  <- rep(NA, 3)
  if (!is.null(fit$hessian)) {
    inv_h <- tryCatch(solve(fit$hessian), error = function(e) NULL)
    if (!is.null(inv_h) && all(diag(inv_h) > 0)) {
      se <- sqrt(diag(inv_h))
    }
  }
  return(list(est = est, se = se))
}

log_prior <- function(params) {
  n     <- params[1]
  alpha <- params[2]
  sigma <- params[3]
  if (n <= 0.01 || alpha <= -0.99 || sigma <= 0.001) return(-Inf)
  
  lp_n     <- dgamma(n, shape = 3, rate = 1, log = TRUE)
  lp_alpha <- dgamma(alpha + 1, shape = 2, rate = 1, log = TRUE)
  lp_sigma <- dgamma(sigma, shape = 2, rate = 2.5, log = TRUE)
  return(lp_n + lp_alpha + lp_sigma)
}

log_posterior <- function(params, q) {
  lp <- log_prior(params)
  if (is.infinite(lp)) return(-Inf)
  return(log_lik_gchi(params, q) + lp)
}

fit_bayes_mcmc <- function(q, n_samples = 1500, burnin = 500, init_params = c(3, 0.75, 0.8)) {
  chain <- matrix(NA, nrow = n_samples, ncol = 3)
  colnames(chain) <- c("n", "alpha", "sigma")
  
  curr_p <- init_params
  curr_lp <- log_posterior(curr_p, q)
  prop_sd <- c(0.12, 0.08, 0.06)
  
  for (i in 1:n_samples) {
    prop_p <- curr_p + rnorm(3, mean = 0, sd = prop_sd)
    prop_lp <- log_posterior(prop_p, q)
    if (log(runif(1)) < (prop_lp - curr_lp)) {
      curr_p <- prop_p
      curr_lp <- prop_lp
    }
    chain[i, ] <- curr_p
  }
  post_chain <- chain[(burnin + 1):n_samples, ]
  
  est <- colMeans(post_chain)
  ci_lower <- apply(post_chain, 2, quantile, probs = 0.025)
  ci_upper <- apply(post_chain, 2, quantile, probs = 0.975)
  return(list(est = est, ci_lower = ci_lower, ci_upper = ci_upper))
}

# ------------------------------------------------------------------------------
# 3. Execution of Scenario 2 Simulation
# ------------------------------------------------------------------------------
true_n     <- 3.0
true_alpha <- 0.75
true_sigma <- 0.80
true_params <- c(n = true_n, alpha = true_alpha, sigma = true_sigma)

sample_sizes <- c(10, 20, 50, 100, 200, 500, 1000)
N_reps       <- 1000

results_list <- list()
res_counter  <- 1

cat("======================================================================\n")
cat(" Running Scenario 2: Gen-Chi^2(n =", true_n, ", alpha =", true_alpha, ", sigma =", true_sigma, ")\n")
cat("======================================================================\n\n")

for (m in sample_sizes) {
  cat(sprintf("Processing Sample Size m = %-4d ... ", m))
  start_time <- Sys.time()
  
  mle_est_mat  <- matrix(NA, nrow = N_reps, ncol = 3)
  mle_ll_mat   <- matrix(NA, nrow = N_reps, ncol = 3)
  mle_ul_mat   <- matrix(NA, nrow = N_reps, ncol = 3)
  
  bayes_est_mat <- matrix(NA, nrow = N_reps, ncol = 3)
  bayes_ll_mat  <- matrix(NA, nrow = N_reps, ncol = 3)
  bayes_ul_mat  <- matrix(NA, nrow = N_reps, ncol = 3)
  
  for (r in 1:N_reps) {
    q_data <- rgchi_power(m, true_n, true_alpha, true_sigma)
    
    mle_res <- fit_mle(q_data, init_params = true_params)
    if (!any(is.na(mle_res$est)) && !any(is.na(mle_res$se))) {
      mle_est_mat[r, ] <- mle_res$est
      mle_ll_mat[r, ]  <- mle_res$est - 1.96 * mle_res$se
      mle_ul_mat[r, ]  <- mle_res$est + 1.96 * mle_res$se
    }
    
    bayes_res <- fit_bayes_mcmc(q_data, n_samples = 1500, burnin = 500, init_params = true_params)
    bayes_est_mat[r, ] <- bayes_res$est
    bayes_ll_mat[r, ]  <- bayes_res$ci_lower
    bayes_ul_mat[r, ]  <- bayes_res$ci_upper
  }
  
  param_names <- c("n", "alpha", "sigma")
  for (p_idx in 1:3) {
    p_name <- param_names[p_idx]
    p_true <- true_params[p_idx]
    
    # MLE Metrics
    valid_mle  <- !is.na(mle_est_mat[, p_idx])
    mle_vals   <- mle_est_mat[valid_mle, p_idx]
    mle_ll     <- mle_ll_mat[valid_mle, p_idx]
    mle_ul     <- mle_ul_mat[valid_mle, p_idx]
    
    results_list[[res_counter]] <- data.frame(
      Sample_Size = m, Method = "MLE", Parameter = p_name, True_Value = p_true,
      Bias = mean(mle_vals) - p_true, MSE = mean((mle_vals - p_true)^2),
      Coverage_Probability = mean(mle_ll <= p_true & mle_ul >= p_true),
      Lower_Interval = mean(mle_ll), Upper_Interval = mean(mle_ul)
    )
    res_counter <- res_counter + 1
    
    # Bayesian Metrics
    bayes_vals <- bayes_est_mat[, p_idx]
    bayes_ll   <- bayes_ll_mat[, p_idx]
    bayes_ul   <- bayes_ul_mat[, p_idx]
    
    results_list[[res_counter]] <- data.frame(
      Sample_Size = m, Method = "Bayesian", Parameter = p_name, True_Value = p_true,
      Bias = mean(bayes_vals) - p_true, MSE = mean((bayes_vals - p_true)^2),
      Coverage_Probability = mean(bayes_ll <= p_true & bayes_ul >= p_true),
      Lower_Interval = mean(bayes_ll), Upper_Interval = mean(bayes_ul)
    )
    res_counter <- res_counter + 1
  }
  end_time <- Sys.time()
  cat(sprintf("Done in %.2f mins\n", difftime(end_time, start_time, units = "mins")))
}

# ------------------------------------------------------------------------------
# 4. Save Scenario 2 Output Files
# ------------------------------------------------------------------------------
final_results_df2 <- do.call(rbind, results_list)
csv_filename2 <- "Gen_ChiSquare_Power_Simulation_Scenario2_Results.csv"
write.csv(final_results_df2, file = csv_filename2, row.names = FALSE)

theme_set(theme_bw(base_size = 12))

p_mse2 <- ggplot(final_results_df2, aes(x = factor(Sample_Size), y = MSE, color = Method, group = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ Parameter, scales = "free_y", labeller = label_bquote(cols = .(Parameter))) +
  scale_color_manual(values = c("MLE" = "#E41A1C", "Bayesian" = "#377EB8")) +
  labs(title = expression(bold("Mean Squared Error (MSE) vs. Sample Size")),
       x = "Sample Size (m)", y = "MSE") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

p_cp2 <- ggplot(final_results_df2, aes(x = factor(Sample_Size), y = Coverage_Probability, color = Method, group = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray30") +
  facet_wrap(~ Parameter, labeller = label_bquote(cols = .(Parameter))) +
  scale_color_manual(values = c("MLE" = "#E41A1C", "Bayesian" = "#377EB8")) +
  ylim(0.85, 1.0) +
  labs(title = expression(bold("95% Coverage Probability vs. Sample Size")),
       x = "Sample Size (m)", y = "Coverage Probability") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

ggsave("Gen_ChiSquare_Scenario2_MSE_Comparison.png", plot = p_mse2, width = 10, height = 5, dpi = 300)
ggsave("Gen_ChiSquare_Scenario2_CP_Comparison.png", plot = p_cp2, width = 10, height = 5, dpi = 300)