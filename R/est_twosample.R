#' @title Sample Size Estimation for Two Independent Proportions
#'
#' @description
#' This function estimates the sample sizes required for comparing two independent proportions
#' (\code{p1} and \code{p2}) to ensure the lower limit of a confidence interval for their
#' difference is no less than a specified threshold (\code{delta}), with a desired assurance probability.
#' The sample sizes are calculated using the Wald and Wilson score methods, and simulations
#' are performed to evaluate empirical assurance and coverage probabilities.
#'
#' @param p1 Numeric. True proportion in group 1 (e.g., experimental group).
#' @param p2 Numeric. True proportion in group 2 (e.g., control group).
#' @param delta Numeric. The threshold for the lower limit of the confidence interval for \code{p1 - p2}.
#' @param r Numeric. Sample size allocation ratio (\code{n1 / n2}) between group 1 and group 2.
#' @param alpha Numeric. Significance level for the two-sided confidence interval. Default is 0.05.
#' @param assurance Numeric. Desired assurance probability (e.g., 0.8 for 80% assurance).
#' @param rep Integer. Number of simulation replications to evaluate assurance and coverage probabilities.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{\code{n1_w, n2_w}}{Required sample sizes for group 1 and group 2 using the Wald method.}
#'   \item{\code{n1_s, n2_s}}{Required sample sizes for group 1 and group 2 using the Wilson score method.}
#'   \item{\code{EAPw}}{Empirical assurance probability (lower limit ≥ \code{delta}) using the Wald method (%).}
#'   \item{\code{ECPw}}{Empirical coverage probability (true difference within CI) using the Wald method (%).}
#'   \item{\code{EAPs}}{Empirical assurance probability using the Wilson score method (%).}
#'   \item{\code{ECPs}}{Empirical coverage probability using the Wilson score method (%).}
#' }
#'
#' @details
#' The function calculates sample sizes for two independent groups to achieve a specified
#' assurance probability that the lower bound of the confidence interval for the difference
#' in proportions (\code{p1 - p2}) is no less than \code{delta}.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item Wald Method: A standard normal approximation is used for the confidence interval.
#'   \item Wilson Score Method: The equation is solved numerically using the \code{uniroot} function for greater accuracy, particularly for small sample sizes or extreme proportions.
#'   \item Simulation: Parallel computing is used to evaluate the empirical assurance and coverage probabilities of the calculated sample sizes.
#' }
#'
#' \strong{Simulation Details:}
#' Simulations generate repeated confidence intervals to evaluate:
#' \itemize{
#'   \item Empirical Assurance Probability (EAP): The proportion of simulated CIs where the lower limit exceeds \code{delta}.
#'   \item Empirical Coverage Probability (ECP): The proportion of simulated CIs containing the true difference (\code{p1 - p2}).
#' }
#'
#' @examples
#' # Example 1: Equal group sizes (r = 1)
#' result <- est_twosample(p1 = 0.97, p2 = 0.98, delta = -0.06, r = 1,
#'                         alpha = 0.05, assurance = 0.97, rep = 10000)
#' print(result)
#'
#' # Example 2: Group size ratio r = 2
#' result <- est_twosample(p1 = 0.97, p2 = 0.98, delta = -0.06, r = 2,
#'                         alpha = 0.05, assurance = 0.97, rep = 10000)
#' print(result)
#'
#' @note
#' \itemize{
#'   \item The Wilson score method is recommended for extreme proportions (close to 0 or 1),
#'         where the Wald method may be less reliable.
#'   \item Parallel computing is employed for large \code{rep} values, so ensure sufficient computational resources.
#' }
#'
#' @import doParallel
#' @import DescTools
#' @export



est_twosample <- function(p1, p2, delta, r, alpha, assurance, rep) {
  # Load necessary libraries
  library(doParallel)
  library(DescTools)

  # Compute sample size for two proportions
  compute_sample_size <- function(p1, p2, delta, alpha, assurance, r) {
    d <- p1 - p2 - delta
    z_alpha <- qnorm(1 - alpha / 2)
    z_assurance <- qnorm(assurance)

    # Wald method
    var_p1 <- p1 * (1 - p1)
    var_p2 <- p2 * (1 - p2)
    n2_w <- ceiling(((var_p1 / r + var_p2) * (z_alpha + z_assurance)^2) / d^2)
    n1_w <- ceiling(r * n2_w)

    # Wilson score method (using uniroot for solving the equation)
    wilson_eq <- function(n2, p1, p2, d, z_alpha, z_assurance, r) {
      term1 <- ((z_alpha * (2 * p1 - 1) + sqrt(z_alpha^2 + 4 * r * n2 * p1 * (1 - p1))) /
                  (2 * (r * n2 + z_alpha^2)))^2
      term2 <- ((z_alpha * (1 - 2 * p2) + sqrt(z_alpha^2 + 4 * n2 * p2 * (1 - p2))) /
                  (2 * (n2 + z_alpha^2)))^2
      term1 + term2 - d^2 / (z_alpha + z_assurance)^2
    }

    n2_s <- ceiling(uniroot(wilson_eq, c(1, 1e6), p1 = p1, p2 = p2, d = d,
                            z_alpha = z_alpha, z_assurance = z_assurance, r = r)$root)
    n1_s <- ceiling(r * n2_s)

    # Exact method
    exact_twosample <- function(p1, p2, delta, r,
                            alpha = 0.05,
                            assurance = 0.8,
                            n_range = 10:1000,
                            plot = TRUE) {
  # Quantiles for CI and assurance
  z1 <- qnorm(1 - alpha / 2)
  z2 <- qnorm(assurance)

  # Left-hand side: theoretical standardized squared distance
  lhs <- (p1 - p2 - delta)^2 / (z1 + z2)^2

  # Right-hand side: squared distance using lower/upper limits of CI
  rhs_values <- sapply(n_range, function(n) {
    n1 <- r * n
    n2 <- n
    k1 <- n1 * p1
    k2 <- n2 * p2

    ci_lower1 <- qbeta(alpha / 2, k1, n1 - k1 + 1)
    ci_upper2 <- qbeta(1 - alpha / 2, k2 + 1, n2 - k2)

    (p1 - ci_lower1)^2 / z1^2 + (ci_upper2 - p2)^2 / z1^2
  })

  # Identify the closest intersection point
  intersection_idx <- which.min(abs(rhs_values - lhs))
  n_intersect <- n_range[intersection_idx]
  rhs_intersect <- rhs_values[intersection_idx]

  # Optional plot
  if (plot) {
    plot(n_range, rhs_values, type = "l", lwd = 2, col = "blue",
         ylab = "RHS", xlab = "Sample Size (n)",
         main = "RHS vs. Sample Size with LHS Reference Line")
    abline(h = lhs, col = "red", lty = 2, lwd = 2)
    points(n_intersect, rhs_intersect, pch = 19, col = "darkgreen")
    text(n_intersect, rhs_intersect,
         labels = paste0("n = ", n_intersect),
         pos = 3, col = "darkgreen")
    legend("topright", legend = c("RHS", "LHS", "Intersection"),
           col = c("blue", "red", "darkgreen"),
           lty = c(1, 2, NA), pch = c(NA, NA, 19), lwd = 2)
  }

  return(list(
    estimated_sample_size = n_intersect,
    lhs = lhs,
    rhs_at_n = rhs_intersect
  ))
}

    h<-exact_twosample(p1 = p1,
                       p2 = p2,
                       delta = delta,
                       r = r,
                       alpha = alpha,
                       assurance = assurance,
                       n_range = 10:10000,
                       plot = FALSE)
    n2_e <- h$estimated_sample_size
    n1_e <- ceiling(r*n2_e)
 
    return(data.frame(n1_w, n2_w, n1_s, n2_s, n1_e, n2_e))
  }
  

  # Simulate confidence intervals
  simulate_ci <- function(p1, p2, alpha, n1_w, n2_w, n1_s, n2_s, n1_e, n2_e) {
    z_alpha <- qnorm(1 - alpha / 2)
    pd <- p1 - p2

    # Wald method
    x1 <- rbinom(n1_w, 1, p1)
    x2 <- rbinom(n2_w, 1, p2)
    p1_hat <- mean(x1)
    p2_hat <- mean(x2)
    ml_w <- p1_hat - p2_hat - z_alpha * sqrt(p1_hat * (1 - p1_hat) / n1_w + p2_hat * (1 - p2_hat) / n2_w)
    mu_w <- p1_hat - p2_hat + z_alpha * sqrt(p1_hat * (1 - p1_hat) / n1_w + p2_hat * (1 - p2_hat) / n2_w)

    # Wilson method
    x1_s <- rbinom(n1_s, 1, p1)
    x2_s <- rbinom(n2_s, 1, p2)
    p1_hat_s <- mean(x1_s)
    p2_hat_s <- mean(x2_s)

    # Wilson score method
    l1 <- (2 * n1_s * p1_hat_s + z_alpha^2 - z_alpha * sqrt(z_alpha^2 + 4 * n1_s * p1_hat_s * (1 - p1_hat_s))) /
      (2 * (n1_s + z_alpha^2))
    u1 <- (2 * n1_s * p1_hat_s + z_alpha^2 + z_alpha * sqrt(z_alpha^2 + 4 * n1_s * p1_hat_s * (1 - p1_hat_s))) /
      (2 * (n1_s + z_alpha^2))

    l2 <- (2 * n2_s * p2_hat_s + z_alpha^2 - z_alpha * sqrt(z_alpha^2 + 4 * n2_s * p2_hat_s * (1 - p2_hat_s))) /
      (2 * (n2_s + z_alpha^2))
    u2 <- (2 * n2_s * p2_hat_s + z_alpha^2 + z_alpha * sqrt(z_alpha^2 + 4 * n2_s * p2_hat_s * (1 - p2_hat_s))) /
      (2 * (n2_s + z_alpha^2))

    ml_s <- p1_hat_s - p2_hat_s - sqrt((p1_hat_s - l1)^2 + (u2 - p2_hat_s)^2)
    mu_s <- p1_hat_s - p2_hat_s + sqrt((u1 - p1_hat_s)^2 + (p2_hat_s - l2)^2)

    # Exact method
    x1_e <- rbinom(n1_e, 1, p1)
    x2_e <- rbinom(n2_e, 1, p2)
    p1_hat_e = mean(x1_e)
    p2_hat_e = mean(x2_e)
    k1_e = sum(x1_e)
    k2_e = sum(x2_e)
    l1_e = qbeta(alpha/2, k1_e, n1_e-k1_e+1)
    u1_e = qbeta(1-alpha/2, k1_e+1, n1_e-k1_e)
    l2_e = qbeta(alpha/2, k2_e, n2_e-k2_e+1)
    u2_e = qbeta(1-alpha/2, k2_e+1, n2_e-k2_e)

    ml_e <- p1_hat_e - p2_hat_e - sqrt((p1_hat_e - l1_e)^2 + (u2_e - p2_hat_e)^2)
    mu_e <- p1_hat_e - p2_hat_e + sqrt((u1_e - p1_hat_e)^2 + (p2_hat_e - l2_e)^2)
    
    return(data.frame(ml_w, mu_w, ml_s, mu_s, ml_e, mu_e))
  }

  # Summarize simulation results
  summarize_results <- function(p1, p2, delta, ci_data) {
    pd <- p1 - p2
    results <- data.frame(
      EAPw = mean(ci_data$ml_w > delta) * 100,
      ECPw = mean(ci_data$ml_w < pd & ci_data$mu_w > pd) * 100,
      EAPs = mean(ci_data$ml_s > delta) * 100,
      ECPs = mean(ci_data$ml_s < pd & ci_data$mu_s > pd) * 100，
      EAPe = mean(ci_data$ml_e > delta) * 100,
      ECPe = mean(ci_data$ml_e < pd & ci_data$mu_e > pd) * 100
    )
    return(results)
  }

  # Parameter grid
  params <- expand.grid(r = r, delta = delta, p1 = p1, p2 = p2, assurance = assurance)
  results <- list()

  for (i in seq_len(nrow(params))) {
    row <- params[i, ]
    sample_sizes <- compute_sample_size(row$p1, row$p2, row$delta, alpha, row$assurance, row$r)

    cl <- makeCluster(6)
    registerDoParallel(cl)
    ci_data <- foreach(j = seq_len(rep), .combine = rbind) %dopar% {
      simulate_ci(row$p1, row$p2, alpha,
                  sample_sizes$n1_w, sample_sizes$n2_w,
                  sample_sizes$n1_s, sample_sizes$n2_s,
                 sample_sizes$n1_e, sample_sizes$n2_e,)
    }
    stopCluster(cl)

    stats <- summarize_results(row$p1, row$p2, row$delta, ci_data)
    results[[i]] <- cbind(row, sample_sizes, stats)
  }

  final_results <- do.call(rbind, results)
  return(final_results)
}

