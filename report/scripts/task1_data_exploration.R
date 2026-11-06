# Task 1: Preliminary data analysis
# Run from the project folder (the one that contains data/, scripts/, figures/)

# load the data files
X    <- read.csv("data/X.csv", header = TRUE)
y    <- read.csv("data/y.csv", header = TRUE)
time <- read.csv("data/time.csv", header = TRUE)

# put everything in one table
df <- data.frame(time = time$time_seconds,
                 x1 = X$cpu_utilization,
                 x2 = X$memory_usage,
                 x3 = X$disk_io,
                 x4 = X$network_traffic,
                 y  = y$power_consumption)

# check for missing values and look at summary
colSums(is.na(df))
summary(df)

# standard deviation and skewness of each signal
skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
round(sapply(df[, 2:6], sd), 3)
round(sapply(df[, 2:6], skew), 3)

dir.create("figures", showWarnings = FALSE)

# ---- Task 1.1: Time series plots ----

# full 1000 s record (10000 samples)
png("figures/fig1_time_series.png", width = 1000, height = 1400, res = 120)
par(mfrow = c(5, 1))

plot(df$time, df$x1, type = "l", main = "x1: CPU utilization", xlab = "Time (s)", ylab = "x1")
plot(df$time, df$x2, type = "l", main = "x2: Memory usage", xlab = "Time (s)", ylab = "x2")
plot(df$time, df$x3, type = "l", main = "x3: Disk I/O", xlab = "Time (s)", ylab = "x3")
plot(df$time, df$x4, type = "l", main = "x4: Network traffic", xlab = "Time (s)", ylab = "x4")
plot(df$time, df$y, type = "l", main = "y: Power consumption", xlab = "Time (s)", ylab = "y")

dev.off()

# zoomed view of the first 20 seconds (200 samples) - the full plot is too dense to read
zoom <- df$time <= 20

png("figures/fig1b_time_series_zoom.png", width = 1000, height = 1400, res = 120)
par(mfrow = c(5, 1))

plot(df$time[zoom], df$x1[zoom], type = "o", pch = 20, cex = 0.5,
     main = "x1: CPU utilization (first 20 s)", xlab = "Time (s)", ylab = "x1")
plot(df$time[zoom], df$x2[zoom], type = "o", pch = 20, cex = 0.5,
     main = "x2: Memory usage (first 20 s)", xlab = "Time (s)", ylab = "x2")
plot(df$time[zoom], df$x3[zoom], type = "o", pch = 20, cex = 0.5,
     main = "x3: Disk I/O (first 20 s)", xlab = "Time (s)", ylab = "x3")
plot(df$time[zoom], df$x4[zoom], type = "o", pch = 20, cex = 0.5,
     main = "x4: Network traffic (first 20 s)", xlab = "Time (s)", ylab = "x4")
plot(df$time[zoom], df$y[zoom], type = "o", pch = 20, cex = 0.5,
     main = "y: Power consumption (first 20 s)", xlab = "Time (s)", ylab = "y")

dev.off()

# autocorrelation at lag 1 - checks whether samples depend on the previous sample
round(sapply(df[, 2:6], function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2]), 3)

# ---- Task 1.2: Distribution of each signal ----

png("figures/fig2_distributions.png", width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3))

hist(df$x1, breaks = 30, freq = FALSE, main = "x1: CPU utilization", xlab = "x1")
lines(density(df$x1), col = "red")

hist(df$x2, breaks = 30, freq = FALSE, main = "x2: Memory usage", xlab = "x2")
lines(density(df$x2), col = "red")

hist(df$x3, breaks = 30, freq = FALSE, main = "x3: Disk I/O", xlab = "x3")
lines(density(df$x3), col = "red")

hist(df$x4, breaks = 30, freq = FALSE, main = "x4: Network traffic", xlab = "x4")
lines(density(df$x4), col = "red")

# y has a long right tail, so the y-axis limit is set from the density curve
d_y <- density(df$y)
hist(df$y, breaks = 60, freq = FALSE, ylim = c(0, max(d_y$y) * 1.1),
     main = "y: Power consumption", xlab = "y")
lines(d_y, col = "red")

dev.off()

# ---- Task 1.3: Correlation and scatter plots ----

# correlation between all signals (time column left out)
cor_pearson <- cor(df[, 2:6], method = "pearson")
round(cor_pearson, 3)

cor_spearman <- cor(df[, 2:6], method = "spearman")
round(cor_spearman, 3)

# correlation heatmap (Pearson) with the values written in each cell
png("figures/fig3a_correlation.png", width = 800, height = 700, res = 120)
par(mar = c(4, 4, 3, 6))
vars <- colnames(cor_pearson)
cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(41)
image(1:5, 1:5, cor_pearson[, 5:1], zlim = c(-1, 1), col = cols,
      axes = FALSE, xlab = "", ylab = "", main = "Pearson correlation matrix")
axis(1, at = 1:5, labels = vars)
axis(2, at = 1:5, labels = rev(vars), las = 1)
for (i in 1:5) for (j in 1:5)
  text(i, 6 - j, sprintf("%.2f", cor_pearson[i, j]), cex = 0.9)
box()
dev.off()

# scatter plots of each input against power
png("figures/fig3_scatter.png", width = 1100, height = 900, res = 120)
par(mfrow = c(2, 2))

plot(df$x1, df$y, pch = ".", main = "y vs x1 (CPU)", xlab = "x1", ylab = "y")
lines(lowess(df$x1, df$y), col = "red", lwd = 2)

plot(df$x2, df$y, pch = ".", main = "y vs x2 (Memory)", xlab = "x2", ylab = "y")
lines(lowess(df$x2, df$y), col = "red", lwd = 2)

plot(df$x3, df$y, pch = ".", main = "y vs x3 (Disk I/O)", xlab = "x3", ylab = "y")
lines(lowess(df$x3, df$y), col = "red", lwd = 2)

plot(df$x4, df$y, pch = ".", main = "y vs x4 (Network)", xlab = "x4", ylab = "y")
lines(lowess(df$x4, df$y), col = "red", lwd = 2)

dev.off()

# baseline linear model with all four inputs
model_lin <- lm(y ~ x1 + x2 + x3 + x4, data = df)
summary(model_lin)

# residual plots for the linear model
png("figures/fig4_linear_residuals.png", width = 1100, height = 500, res = 120)
par(mfrow = c(1, 2))

plot(fitted(model_lin), resid(model_lin), pch = ".",
     main = "Residuals vs fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(resid(model_lin), pch = ".", main = "Q-Q plot of residuals")
qqline(resid(model_lin), col = "red")

dev.off()
