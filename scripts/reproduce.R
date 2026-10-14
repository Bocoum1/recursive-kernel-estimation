# Read only the original function definitions, not the costly notebook experiments.
path <- "analysis/coursework.Rmd"
if (!file.exists(path)) stop("Run from the repository root: Rscript scripts/reproduce.R")
lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
inside <- FALSE
code <- character()
for (line in lines) {
  if (grepl("^```\\s*\\{r", line)) {
    inside <- TRUE
  } else if (grepl("^```", line)) {
    inside <- FALSE
  } else if (inside) {
    code <- c(code, line)
  }
}
needed <- c("f", "fprime", "Epach", "kprime_Epach", "quadratique",
            "kprime_quadratique", "gaussian", "kprime_gaussian", "h_hat",
            "h_hatprime", "g_hat", "g_hatprime", "f_hatprime", "f_tildeprime", "f_qprime")
course <- new.env()
for (expr in parse(text = code)) {
  if (is.call(expr) && is.name(expr[[1]]) &&
      as.character(expr[[1]]) %in% c("=", "<-") && is.name(expr[[2]]) &&
      as.character(expr[[2]]) %in% needed && is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function"))) eval(expr, course)
}
stopifnot(all(needed %in% ls(course)))
kernels <- list(
  Epanechnikov = list(k = course$Epach, kp = course$kprime_Epach),
  Biweight = list(k = course$quadratique, kp = course$kprime_quadratique),
  Gaussian = list(k = course$gaussian, kp = course$kprime_gaussian)
)
grid <- seq(0.1, 0.9, length.out = 100)
eps <- 1e-6
stopifnot(max(abs((course$f(grid + eps) - course$f(grid - eps)) / (2 * eps) -
                  course$fprime(grid))) < 1e-5)
for (name in names(kernels)) {
  bound <- if (name == "Gaussian") Inf else 1
  stopifnot(abs(integrate(kernels[[name]]$k, -bound, bound)$value - 1) < 1e-6)
}
estimate <- function(x, X, Y, kernel, alpha) {
  c(eq9 = course$f_hatprime(x, X, Y, alpha, kernel$k, kernel$kp),
    eq10 = course$f_tildeprime(x, X, Y, alpha, function(z) 1,
                             function(z) 0, kernel$k, kernel$kp),
    eq11 = course$f_qprime(x, X, Y, alpha, function(z) 1, kernel$kp))
}
n <- 10000L
repetitions <- 200L
alpha <- 0.3
seed <- 20261014L
points <- c(0.5, 0.8)
set.seed(seed)
errors <- array(NA_real_, dim = c(repetitions, length(points), 3, 3),
                dimnames = list(NULL, points, names(kernels), c("eq9", "eq10", "eq11")))
for (b in seq_len(repetitions)) {
  X <- runif(n)
  Y <- course$f(X) + rnorm(n)
  for (j in seq_along(points)) for (k in seq_along(kernels)) {
    fit <- estimate(points[j], X, Y, kernels[[k]], alpha)
    stopifnot(all(is.finite(fit)), abs(fit[2] - fit[3]) < 1e-10)
    errors[b, j, k, ] <- fit - course$fprime(points[j])
  }
}
rows <- list()
for (j in seq_along(points)) for (k in seq_along(kernels)) for (m in 1:3) {
  error <- errors[, j, k, m]
  rows[[length(rows) + 1L]] <- data.frame(
    x = points[j], kernel = names(kernels)[k], estimator = dimnames(errors)[[4]][m],
    mse = mean(error^2), bias = mean(error),
    mse_mc_se = sd(error^2) / sqrt(repetitions)
  )
}
result <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(result, "results/monte-carlo.csv", row.names = FALSE)
writeLines(c(paste("n =", n), paste("repetitions =", repetitions),
             paste("alpha =", alpha), paste("seed =", seed), R.version.string),
           "results/config.txt")

# A fresh fixed draw for the illustration; do not recenter or hide boundary errors.
set.seed(seed + 1L)
X <- runif(n)
Y <- course$f(X) + rnorm(n)
xx <- seq(0, 1, length.out = 250)
curves <- vapply(xx, function(x) estimate(x, X, Y, kernels[[1]], alpha), numeric(3))
png("results/derivative.png", width = 1440, height = 900, res = 150)
par(mar = c(4.5, 4.5, 3, 1), bg = "#faf8f2", fg = "#102d2b")
plot(xx, course$fprime(xx), type = "l", lwd = 2.5,
     ylim = range(curves, course$fprime(xx)), xlab = "x", ylab = "Derivative",
     main = "Recursive derivative estimation | n = 10,000")
grid(col = "#e3e3dc")
lines(xx, curves[1, ], col = "#176076", lwd = 2)
lines(xx, curves[2, ], col = "#df5b2d", lwd = 2, lty = 2)
legend("bottomleft", c("True derivative", "Estimated density (9)", "Known density (10 = 11)"),
       col = c("#102d2b", "#176076", "#df5b2d"), lwd = 2, lty = c(1, 1, 2),
       bg = "#faf8f2", cex = 0.8)
dev.off()
print(result, row.names = FALSE, digits = 5)
cat("PASS: derivative, kernel integrals and uniform-density identity verified.\n")
