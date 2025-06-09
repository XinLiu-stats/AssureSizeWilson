#' @title Sample Size Estimation for a Single Proportion
#'
#' @description
#' This function estimates the sample size required to ensure the lower bound of a confidence interval for a single proportion is no less than a specified threshold, with a desired level of assurance, using the Wald, Wilson score, and Exact methods.
#' It evaluates empirical assurance probabilities and empirical coverage probabilities through simulations.
#'
#' @param p Numeric. The true population proportion to be estimated. Must be in the range (0, 1).
#' @param p0 Numeric. The threshold proportion for evaluating the lower bound of the CI. Must be in the range (0, 1).
#' @param alpha Numeric. The significance level for the confidence interval (e.g., 0.05 for 95% confidence). Must be in the range (0, 1).
#' @param assurance Numeric. The desired assurance probability (e.g., 0.8 for 80% assurance). Must be in the range (0, 1).
#' @param rep Integer. The number of simulation replications to evaluate assurance and coverage probabilities. Must be a positive integer.
#'
#' @details
#' This function estimates sample sizes using the Wald and Wilson score methods. It simulates binomial data to construct confidence intervals
#' and evaluates empirical assurance probabilities and coverage probabilities through the following steps:
#' \enumerate{
#'   \item Computes sample sizes using the Wald, Wilson score, and Exact methods.
#'   \item Simulates confidence intervals for the specified sample sizes.
#'   \item Calculates assurance probability (percentage of simulations where the lower CI bound exceeds \code{p0}) and
#'         coverage probability (percentage of simulations where the true proportion lies within the CI).
#' }
#'
#' The Wilson score method generally provides more accurate sample size than the Wald method, especially when \code{p} is near 0 or 1.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{p0}: The specified threshold proportion.
#'   \item \code{p}: The true population proportion.
#'   \item \code{assurance}: The desired assurance probability.
#'   \item \code{Nw}: Sample size estimated using the Wald method.
#'   \item \code{Ns}: Sample size estimated using the Wilson score method.
#'   \item \code{EAPw}: Empirical assurance probability for the Wald method (%).
#'   \item \code{EAPs}: Empirical assurance probability for the Wilson score method (%).
#'   \item \code{EAPe}: Empirical assurance probability for the Exact method (%).
#'   \item \code{ECPw}: Empirical coverage probability for the Wald method (%).
#'   \item \code{ECPs}: Empirical coverage probability for the Wilson score method (%).
#'   \item \code{ECPe}: Empirical coverage probability for the Exact method (%).
#' }
#'
#' @examples
#' # Example 1: High true proportion
#' result <- est_onesample(
#'   p = 0.937,
#'   p0 = 0.88,
#'   alpha = 0.05,
#'   assurance = 0.8,
#'   rep = 10000
#' )
#' print(result)
#'
#' # Example 2: Moderate true proportion
#' result <- est_onesample(
#'   p = 0.6,
#'   p0 = 0.5,
#'   alpha = 0.05,
#'   assurance = 0.9,
#'   rep = 10000
#' )
#' print(result)
#'
#' @note
#' \itemize{
#'   \item \code{p0} represents the minimum acceptable value for the lower CI bound.
#'   \item For extreme proportions, the Wilson score method is generally more reliable.
#'   \item Larger \code{rep} values yield more precise results but require greater computation time.
#' }
#'
#' @import doParallel
#' @import DescTools
#' @export



est_onesample <- function(p, p0, alpha, assurance, rep) {
  # Load necessary libraries
  library(doParallel)
  library(DescTools)
  # Function to calculate sample size
  calculate_sample_size <- function(p, p0, alpha, assurance) {
    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(assurance)
    d <- p - p0

    # Wald method
    Nw <- ceiling((z1 + z2)^2 * p * (1 - p) / d^2)

    # Wilson score method
    Ns <- ceiling((z1 + z2)^2 * p * (1 - p) / d^2 + z1 * (z1 + z2) * (2 * p - 1) / d - z1^2)

    # Exact method
    exact_sample <- function(p, p0, alpha = 0.05, assurance = 0.8,
                         n_range = 10:1000, plot = TRUE) {
  # Quantiles for CI and assurance
  z1 <- qnorm(1 - alpha / 2)
  z2 <- qnorm(assurance)

  # Left-hand side (target standardized squared distance)
  lhs <- (p - p0)^2 / (z1 + z2)^2

  # Right-hand side for each candidate sample size (using beta approximation)
  rhs_values <- sapply(n_range, function(n) {
    k <- n * p
    ci_lower <- qbeta(alpha / 2, k, n - k + 1)
    (p - ci_lower)^2 / z1^2
  })

  # Find the closest intersection point
  intersection_idx <- which.min(abs(rhs_values - lhs))
  n_intersect <- n_range[intersection_idx]
  rhs_intersect <- rhs_values[intersection_idx]

  # Plot results if required
  if (plot) {
    plot(n_range, rhs_values, type = "l", lwd = 2, col = "blue",
         ylab = "RHS", xlab = "Sample Size (n)",
         main = "RHS vs. Sample Size with LHS Reference")
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

    h<-exact_sample(p = p,
                    p0 = p0,
                    alpha = alpha,
                    assurance = assurance,
                    n_range = 10:10000,
                    plot = FALSE)
    
    Ne <- h$estimated_sample_size
    
    return(c(Nw = Nw, Ns = Ns, Ne = Ne))
  }

  # Function to simulate confidence intervals
  simulate_confidence_intervals <- function(p, N, alpha) {
    Nw <- N["Nw"]
    Ns <- N["Ns"]
    Ne <- N["Ne"]

    # Wald method
    x1 <- rbinom(Nw, 1, p)
    p1 <- mean(x1)
    l1 <- p1 - qnorm(1 - alpha / 2) * sqrt(p1 * (1 - p1) / Nw)
    u1 <- p1 + qnorm(1 - alpha / 2) * sqrt(p1 * (1 - p1) / Nw)

    # Wilson score method
    x2 <- rbinom(Ns, 1, p)
    p2 <- mean(x2)
    a <- Ns + qnorm(1 - alpha / 2)^2
    b <- -(2 * Ns * p2 + qnorm(1 - alpha / 2)^2)
    c <- Ns * p2^2
    l2 <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
    u2 <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)

    # Exact method
    x3 <- rbinom(Ne, 1, p)
    k <- sum(x3)
    l3 <- qbeta(alpha / 2, k, Ne - k + 1)
    u3 <- qbeta(1- alpha/2, k + 1, Ne - k )
    return(c(l1 = l1, l2 = l2, l3 = l3, u1 = u1, u2 = u2, u3 = u3))
  }

  # Function to calculate assurance probabilities
  calculate_assurance <- function(p, p0, ci) {
    data <- as.data.frame(ci)
    colnames(data) <- c("l1", "l2", "l3", "u1", "u2", "u3")

    EAPw <- round(mean(data$l1 > p0, na.rm = TRUE) * 100, 2)
    EAPs <- round(mean(data$l2 > p0, na.rm = TRUE) * 100, 2)
    EAPe <- round(mean(data$l3 > p0, na.rm = TRUE) * 100, 2)
    ECPw <- round(mean(data$l1 < p & data$u1 > p, na.rm = TRUE) * 100, 2)
    ECPs <- round(mean(data$l2 < p & data$u2 > p, na.rm = TRUE) * 100, 2)
    ECPe <- round(mean(data$l3 < p & data$u3 > p, na.rm = TRUE) * 100, 2)
    
    return(c(EAPw = EAPw, EAPs = EAPs, EAPe = EAPe, ECPw = ECPw, ECPs = ECPs, ECPe = ECPe))
  }

  # Generate parameter grid
  param_grid <- expand.grid(p0 = p0, p = p, assurance = assurance)
  set_num <- nrow(param_grid)

  # Initialize results
  results <- list()

  # Simulation loop
  for (i in 1:set_num) {
    p0 <- param_grid[i, "p0"]
    p <- param_grid[i, "p"]
    assurance <- param_grid[i, "assurance"]

    # Calculate sample sizes
    N <- calculate_sample_size(p, p0, alpha, assurance)

    # Parallel processing for confidence interval simulation
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    ci <- foreach(x = 1:rep, .combine = 'rbind') %dopar% simulate_confidence_intervals(p, N, alpha)
    stopCluster(cl)

    # Calculate assurance probabilities
    assurance_stats <- calculate_assurance(p, p0, ci)

    # Store results
    results[[i]] <- c(param_grid[i, ], N, assurance_stats)
  }

  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  colnames(results_df) <- c("p0", "p", "assurance", "Nw", "Ns", "Ne", "EAPw", "EAPs", "EAPe", "ECPw", "ECPs", "ECPe")

  return(as.data.frame(results_df))
}
