# Task 2: Polynomial regression models
# load the data files
X    <- read.csv("data/X.csv", header = TRUE)
y    <- read.csv("data/y.csv", header = TRUE)
time <- read.csv("data/time.csv", header = TRUE)

df <- data.frame(time = time$time_seconds,
                 x1 = X$cpu_utilization,
                 x2 = X$memory_usage,
                 x3 = X$disk_io,
                 x4 = X$network_traffic,
                 y  = y$power_consumption)

head(df)

# ---- Task 2.1: Least squares estimation ----

# output as a plain vector
y_vec <- df$y

# column of ones for the bias term
ones <- rep(1, nrow(df))

# Model 1: x4, x1^2, x1^3, x2^4, x1^4, bias
X_m1 <- cbind(x4 = df$x4, x1_2 = df$x1^2, x1_3 = df$x1^3,
              x2_4 = df$x2^4, x1_4 = df$x1^4, bias = ones)

theta_m1 <- solve(t(X_m1) %*% X_m1) %*% t(X_m1) %*% y_vec
theta_m1

# check Model 1 against lm() - estimates should match theta_m1
check_m1 <- lm(y ~ x4 + I(x1^2) + I(x1^3) + I(x2^4) + I(x1^4), data = df)
coef(check_m1)

# Model 2: x4, x1^3, x3^4, bias
X_m2 <- cbind(x4 = df$x4, x1_3 = df$x1^3, x3_4 = df$x3^4, bias = ones)
theta_m2 <- solve(t(X_m2) %*% X_m2) %*% t(X_m2) %*% y_vec
theta_m2

# Model 3: x3^3, x3^4, bias
X_m3 <- cbind(x3_3 = df$x3^3, x3_4 = df$x3^4, bias = ones)
theta_m3 <- solve(t(X_m3) %*% X_m3) %*% t(X_m3) %*% y_vec
theta_m3

# Model 4: x2, x1^3, x3^4, bias
X_m4 <- cbind(x2 = df$x2, x1_3 = df$x1^3, x3_4 = df$x3^4, bias = ones)
theta_m4 <- solve(t(X_m4) %*% X_m4) %*% t(X_m4) %*% y_vec
theta_m4

# Model 5: x4, x1^2, x1^3, x3^4, bias
X_m5 <- cbind(x4 = df$x4, x1_2 = df$x1^2, x1_3 = df$x1^3,
              x3_4 = df$x3^4, bias = ones)
theta_m5 <- solve(t(X_m5) %*% X_m5) %*% t(X_m5) %*% y_vec
theta_m5

# ---- Task 2.2: Residual sum of squares (RSS) ----

# predictions from each model
yhat_m1 <- X_m1 %*% theta_m1
yhat_m2 <- X_m2 %*% theta_m2
yhat_m3 <- X_m3 %*% theta_m3
yhat_m4 <- X_m4 %*% theta_m4
yhat_m5 <- X_m5 %*% theta_m5

# RSS = sum of squared errors
rss_m1 <- sum((y_vec - yhat_m1)^2)
rss_m2 <- sum((y_vec - yhat_m2)^2)
rss_m3 <- sum((y_vec - yhat_m3)^2)
rss_m4 <- sum((y_vec - yhat_m4)^2)
rss_m5 <- sum((y_vec - yhat_m5)^2)

rss_all <- c(rss_m1, rss_m2, rss_m3, rss_m4, rss_m5)
rss_all

# ---- Task 2.3: Log-likelihood ----

n <- length(y_vec)            # number of samples (10000)

# estimated noise variance for each model
sigma2_all <- rss_all / (n - 1)

# log-likelihood for each model
loglik_all <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2_all) - rss_all / (2 * sigma2_all)

sigma2_all
loglik_all

# ---- Task 2.4: AIC and BIC ----

# number of parameters in each model (including bias)
k_all <- c(ncol(X_m1), ncol(X_m2), ncol(X_m3), ncol(X_m4), ncol(X_m5))
k_all

aic_all <- 2 * k_all - 2 * loglik_all
bic_all <- k_all * log(n) - 2 * loglik_all

# summary table of all results so far
results <- data.frame(model  = 1:5,
                      k      = k_all,
                      RSS    = rss_all,
                      logLik = loglik_all,
                      AIC    = aic_all,
                      BIC    = bic_all)
results

# ---- Task 2.5: Residual distributions ----

# residuals (errors) for each model
res_m1 <- as.vector(y_vec - yhat_m1)
res_m2 <- as.vector(y_vec - yhat_m2)
res_m3 <- as.vector(y_vec - yhat_m3)
res_m4 <- as.vector(y_vec - yhat_m4)
res_m5 <- as.vector(y_vec - yhat_m5)

png("figures/fig5_residuals.png", width = 1000, height = 1800, res = 120)
par(mfrow = c(5, 2))

hist(res_m1, breaks = 50, main = "Model 1: residuals", xlab = "Residual")
qqnorm(res_m1, pch = ".", main = "Model 1: Q-Q plot"); qqline(res_m1, col = "red")

hist(res_m2, breaks = 50, main = "Model 2: residuals", xlab = "Residual")
qqnorm(res_m2, pch = ".", main = "Model 2: Q-Q plot"); qqline(res_m2, col = "red")

hist(res_m3, breaks = 50, main = "Model 3: residuals", xlab = "Residual")
qqnorm(res_m3, pch = ".", main = "Model 3: Q-Q plot"); qqline(res_m3, col = "red")

hist(res_m4, breaks = 50, main = "Model 4: residuals", xlab = "Residual")
qqnorm(res_m4, pch = ".", main = "Model 4: Q-Q plot"); qqline(res_m4, col = "red")

hist(res_m5, breaks = 50, main = "Model 5: residuals", xlab = "Residual")
qqnorm(res_m5, pch = ".", main = "Model 5: Q-Q plot"); qqline(res_m5, col = "red")

dev.off()

# Shapiro-Wilk normality test on a random 5000 residuals (the test allows at most 5000)
set.seed(1)
sw_idx <- sample(1:n, 5000)
sapply(list(m1 = res_m1, m2 = res_m2, m3 = res_m3, m4 = res_m4, m5 = res_m5),
       function(r) shapiro.test(r[sw_idx])$p.value)

# ---- Task 2.6: Model selection ----

# lowest AIC and lowest BIC
which.min(aic_all)
which.min(bic_all)

# difference from the best model (smaller = better, 0 = best)
results$dAIC <- aic_all - min(aic_all)
results$dBIC <- bic_all - min(bic_all)
results

# ---- Task 2.7: Train/test split and prediction ----

set.seed(123)   # makes the random split the same every time the script runs

# randomly choose 70% of the row numbers for training
train_idx <- sample(1:n, size = round(0.7 * n))

# Model 2 design matrix, split by the SAME rows for X and y
X_train <- X_m2[train_idx, ]
y_train <- y_vec[train_idx]
X_test  <- X_m2[-train_idx, ]
y_test  <- y_vec[-train_idx]

nrow(X_train)   # should be 7000
nrow(X_test)    # should be 3000

# 1. estimate parameters using the training data only
theta_train <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
theta_train

# 2. predict the test data
yhat_test <- as.vector(X_test %*% theta_train)

# how good are the test predictions?
rmse_test <- sqrt(mean((y_test - yhat_test)^2))
r2_test   <- 1 - sum((y_test - yhat_test)^2) / sum((y_test - mean(y_test))^2)
rmse_test
r2_test

# 3. 95% confidence intervals for the test predictions

k <- ncol(X_train)                      # number of parameters (4)
rss_train    <- sum((y_train - X_train %*% theta_train)^2)
sigma2_train <- rss_train / (nrow(X_train) - k)   # noise variance

cov_theta <- sigma2_train * solve(t(X_train) %*% X_train)   # parameter uncertainty

# variance of each test prediction: x*' Cov(theta) x*, done for all rows at once
var_fit <- rowSums((X_test %*% cov_theta) * X_test)
se_fit  <- sqrt(var_fit)

t_crit <- qt(0.975, df = nrow(X_train) - k)   # t-value for 95%

ci_lower <- yhat_test - t_crit * se_fit
ci_upper <- yhat_test + t_crit * se_fit

# prediction interval: also includes the noise in a single new observation
pi_lower <- yhat_test - t_crit * sqrt(var_fit + sigma2_train)
pi_upper <- yhat_test + t_crit * sqrt(var_fit + sigma2_train)

t_crit
summary(ci_upper - ci_lower)   # width of the confidence intervals
summary(pi_upper - pi_lower)   # width of the prediction intervals

# what share of actual test values fall inside each interval?
mean(y_test >= ci_lower & y_test <= ci_upper)
mean(y_test >= pi_lower & y_test <= pi_upper)

# 4. plot predictions, 95% CIs and test data (first 40 test points)
idx    <- 1:40
t_test <- df$time[-train_idx][idx]     # time of each test point

png("figures/fig6_confidence_intervals.png", width = 1200, height = 1000, res = 120)
par(mfrow = c(2, 1))

# panel 1: actual values, predictions and CI error bars
plot(t_test, y_test[idx], pch = 16, col = "grey40",
     main = "Model 2: test data, predictions and 95% confidence intervals",
     xlab = "Time (s)", ylab = "Power y")
points(t_test, yhat_test[idx], pch = 4, col = "red")
# the CI error bars are under 1 W wide, too short to draw on this scale,
# so they are shown in the zoomed panel below
legend("topright", legend = c("Test data", "Prediction"),
       pch = c(16, 4), col = c("grey40", "red"), cex = 0.8)

# panel 2: zoomed view - distance of each point from its prediction
dev_y  <- y_test[idx] - yhat_test[idx]
plot(t_test, dev_y, pch = 16, col = "grey40",
     ylim = range(c(dev_y, pi_lower[idx] - yhat_test[idx], pi_upper[idx] - yhat_test[idx])),
     main = "Zoomed view: test data relative to prediction",
     xlab = "Time (s)", ylab = "Actual - predicted (W)")
abline(h = 0, col = "red")
arrows(t_test, ci_lower[idx] - yhat_test[idx], t_test, ci_upper[idx] - yhat_test[idx],
       angle = 90, code = 3, length = 0.03, col = "blue")
lines(t_test, pi_lower[idx] - yhat_test[idx], lty = 2, col = "darkgreen")
lines(t_test, pi_upper[idx] - yhat_test[idx], lty = 2, col = "darkgreen")
legend("topright", legend = c("Test data", "Prediction (0)", "95% CI", "95% prediction interval"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 2),
       col = c("grey40", "red", "blue", "darkgreen"), cex = 0.8)

dev.off()

# predicted vs actual for all 3000 test points
png("figures/fig6b_pred_vs_actual.png", width = 800, height = 750, res = 120)
plot(y_test, yhat_test, pch = 20, cex = 0.4, col = "grey30",
     main = "Model 2: predicted vs actual (test set)",
     xlab = "Actual power y (W)", ylab = "Predicted power (W)")
abline(0, 1, col = "red", lwd = 2)
legend("topleft", legend = c("Test point", "y = x line"),
       pch = c(20, NA), lty = c(NA, 1), col = c("grey30", "red"), cex = 0.8)
dev.off()
