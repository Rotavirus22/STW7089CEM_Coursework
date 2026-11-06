# Task 3: Approximate Bayesian Computation (rejection ABC) for Model 2
# Needs the objects from Task 2 (df, y_vec, X_m2, theta_m2), so run Task 2 first if they are missing
if (!exists("theta_m2")) source("scripts/task2_regression.R")

# ---- Task 3.1: Choose the 2 parameters with the largest |theta| ----
theta_m2
order(abs(theta_m2), decreasing = TRUE)

# ---- Task 3.2: Uniform priors around the estimates ----
# standard errors of the least squares estimates, used to set the prior width
n <- length(y_vec)
k <- ncol(X_m2)

rss_m2    <- sum((y_vec - X_m2 %*% theta_m2)^2)
sigma2_m2 <- rss_m2 / (n - k)

se_m2 <- sqrt(diag(sigma2_m2 * solve(t(X_m2) %*% X_m2)))
se_m2

# ---- Task 3.3: Sample from the priors and run rejection ABC ----

# fixed parameters (x1^3 and x3^4 terms stay at their least squares values)
fixed_part <- theta_m2[2] * df$x1^3 + theta_m2[3] * df$x3^4

# uniform priors: estimate +/- 10 standard errors
prior_bias <- c(theta_m2[4] - 10 * se_m2[4], theta_m2[4] + 10 * se_m2[4])
prior_x4   <- c(theta_m2[1] - 10 * se_m2[1], theta_m2[1] + 10 * se_m2[1])
prior_bias
prior_x4

# draw samples from the priors
set.seed(456)
N <- 100000
bias_samples <- runif(N, prior_bias[1], prior_bias[2])
x4_samples   <- runif(N, prior_x4[1],   prior_x4[2])

# tolerance: accept fits within 0.05% of the least squares RSS
eps <- rss_m2 * 1.0005

# compute the RSS for every sample
rss_samples <- numeric(N)
for (i in 1:N) {
  y_sim <- bias_samples[i] + x4_samples[i] * df$x4 + fixed_part
  rss_samples[i] <- sum((y_vec - y_sim)^2)
}

# keep the samples that fit well enough
accept    <- rss_samples <= eps
post_bias <- bias_samples[accept]
post_x4   <- x4_samples[accept]

sum(accept)     # number accepted
mean(accept)    # acceptance rate

# ---- Task 3.5: Posterior summaries (used to explain the results) ----

# posterior mean and 95% credible interval for each parameter
mean(post_bias)
quantile(post_bias, c(0.025, 0.975))

mean(post_x4)
quantile(post_x4, c(0.025, 0.975))

# correlation between the two parameters in the posterior
cor(post_bias, post_x4)

# ---- Task 3.4: Plot joint and marginal posteriors ----
png("figures/fig7_abc_posterior.png", width = 1500, height = 500, res = 120)
par(mfrow = c(1, 3))

# joint posterior
plot(post_bias, post_x4, pch = ".", col = "steelblue",
     main = "Joint posterior", xlab = "theta_bias", ylab = "theta_x4")
points(theta_m2[4], theta_m2[1], pch = 4, col = "red", cex = 2, lwd = 2)

# marginal posterior of theta_bias
hist(post_bias, breaks = 40, freq = FALSE, col = "lightgrey",
     main = "Marginal posterior: theta_bias", xlab = "theta_bias")
abline(v = theta_m2[4], col = "red", lwd = 2)

# marginal posterior of theta_x4
hist(post_x4, breaks = 40, freq = FALSE, col = "lightgrey",
     main = "Marginal posterior: theta_x4", xlab = "theta_x4")
abline(v = theta_m2[1], col = "red", lwd = 2)

dev.off()