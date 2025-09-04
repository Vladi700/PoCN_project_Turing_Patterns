suppressPackageStartupMessages({
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

set.seed(8)

# --- Output dir ---------------------------------------------------------------
out_dir <- "figs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Graph --------------------------------------------------------------------
N <- 200
m <- 5
g <- sample_pa(n = N, m = m, directed = FALSE)

A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
D <- diag(degree(g))
L <- A - D

# Plot: raw network
png(file.path(out_dir, "01_network_raw.png"), width = 1600, height = 1200, res = 150)
plot(g, vertex.size = 5, vertex.label = NA, edge.arrow.size = 0)
dev.off()

# --- Eigen decomposition -------------------------------------------------------
ev <- eigen(L, symmetric = TRUE)
Lambda_vals <- ev$values       # Laplacian eigenvalues (≤ 0 for this convention)
Phi_vals    <- ev$vectors

# --- Growth-rate (λ+) function & params ---------------------------------------
lambda_plus <- function(Lambda, fu, gv, fv, gu, sigma, eps,
                        complex_ok = TRUE, tol = 1e-12) {
  A    <- fu + gv + (1 + sigma) * eps * Lambda
  disc <- 4 * fv * gu + (fu - gv + (1 - sigma) * eps * Lambda)^2
  B <- if (complex_ok) sqrt(as.complex(disc)) else sqrt(pmax(disc, 0))
  0.5 * (A + B)
}

fu <- 3.33; fv <- -5; gu <- 10; gv <- -4

sigma_c <- function(fu, gv, fv, gu) {
  disc <- fv * gu * (fv * gu - fu * gv)
  (fu * gv - 2 * fv * gu + 2 * sqrt(disc)) / (fu^2)
}
sigma_c_val <- sigma_c(fu, gv, fv, gu)
cat(sprintf("sigma_c = %.6f\n", Re(sigma_c_val)))

# Modes at (sigma, eps) used for curve (from your notebook cell)
sigma_curve <- 15.5
eps_curve   <- 0.425

# Dispersion curve over Λ
Lambda_seq <- seq(min(Lambda_vals), -1e-4, length.out = 10000)
y_curve <- Re(lambda_plus(Lambda_seq, fu, gv, fv, gu, sigma_curve, eps_curve))

png(file.path(out_dir, "02_dispersion_curve.png"), width = 1600, height = 1200, res = 150)
plot(log(-Lambda_seq), y_curve, type = "l",
     xlab = "log(-Λ)", ylab = "Re(λ+)")
abline(h = 0, lty = 2)
dev.off()

# Index of eigenmode with max growth (for the curve parameters)
ac_idx <- which.max(Re(lambda_plus(Lambda_vals, fu, gv, fv, gu,
                                   sigma = sigma_curve, eps = eps_curve)))
cat(sprintf("Critical mode index (arg max growth for curve params): %d\n", ac_idx))

# --- Color by eigenvector (thresholded) ---------------------------------------
color_by_phi <- function(phi_ac, t = 0.1) {
  ifelse(phi_ac >= t, "red",
         ifelse(phi_ac <= -t, "blue", "yellow"))
}

thr <- 0.1
phi <- Phi_vals[, ac_idx]   # use the critical mode index found above

cls <- cut(phi, breaks = c(-Inf, -thr, thr, Inf),
           labels = c("≤ -0.1", "(-0.1, 0.1)", "≥ 0.1"))

df_phi <- data.frame(node = seq_along(phi), phi = phi, cls = cls)

p_phi <- ggplot(df_phi, aes(node, phi, color = cls)) +
  geom_point(size = 2) +
  geom_hline(yintercept = c(-thr, thr), linetype = "dashed") +
  scale_color_manual(values = c("≤ -0.1" = "blue",
                                "(-0.1, 0.1)" = "yellow",
                                "≥ 0.1" = "red"),
                     na.value = "grey80") +
  labs(x = "Node index", y = expression(phi[alpha[c]]),
       color = "Class") +
  theme_minimal()

ggsave(file.path(out_dir, "03_phi_scatter.png"), p_phi, width = 8, height = 5, dpi = 150)

# Plot: network colored by phi
V(g)$color <- color_by_phi(phi, t = thr)
png(file.path(out_dir, "04_network_colored_by_phi.png"), width = 1600, height = 1200, res = 150)
plot(g, vertex.color = V(g)$color, vertex.size = 6, vertex.label = NA, edge.arrow.size = 0)
dev.off()

# --- Mimura–Murray dynamics on the network -----------------------------------
a <- 35; b <- 16; c <- 9; d <- 2/5
f  <- function(u, v) ((a + b * u - u * u) / c - v) * u
gF <- function(u, v) (u - (1 + d * v)) * v

eps_sim   <- 0.12
sigma_sim <- 15.6

set.seed(8)
u <- rep(5, N)  + 1e-2 * runif(N)
v <- rep(10, N) + 1e-2 * runif(N)

steps <- 200
dt    <- 0.005
U <- matrix(NA_real_, nrow = steps + 1, ncol = N); U[1, ] <- u
V <- matrix(NA_real_, nrow = steps + 1, ncol = N); V[1, ] <- v

for (step in 1:steps) {
  du <- f(u, v) + eps_sim        * (L %*% u)
  dv <- gF(u, v) + sigma_sim*eps_sim * (L %*% v)
  u  <- u + dt * as.vector(du)
  v  <- v + dt * as.vector(dv)
  U[step + 1, ] <- u
  V[step + 1, ] <- v
}

# Rank by degree (high -> low) and plot u at final time
k   <- igraph::degree(g)
ord <- order(-k, seq_along(k))
t_idx <- nrow(U)

df_u <- data.frame(
  rank   = seq_along(ord),
  node   = ord,
  degree = k[ord],
  u      = U[t_idx, ord]
)

p_u <- ggplot(df_u, aes(rank, u)) +
  geom_point(size = 2) +
  labs(x = "Rank by degree (high → low)", y = "Activator density u(t)") +
  theme_minimal()

ggsave(file.path(out_dir, "05_u_vs_degree_rank.png"), p_u, width = 8, height = 5, dpi = 150)

# Degree vs rank (base plot)
png(file.path(out_dir, "06_degree_by_rank.png"), width = 1600, height = 1200, res = 150)
plot(seq_along(k), k[ord], type = "b", pch = 19,
     xlab = "Rank (high → low)", ylab = "Degree")
dev.off()

cat("Saved figures to ", normalizePath(out_dir), "\n", sep = "")
