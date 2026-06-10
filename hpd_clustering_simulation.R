#!/usr/bin/env Rscript



# ==============================================================================
# Configuration & Dependencies
# ==============================================================================

# User configuration for Monte Carlo replications
# Set to 15 for full replication (paper results); set to 2 for a fast test run.
N_MC_RUNS <- 15

# Automatically install missing dependencies
required_packages <- c("astsa", "mclust", "cluster")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

# Load libraries
suppressPackageStartupMessages({
  library(astsa)
  library(mclust)
  library(cluster)
})

# Colors for plotting
col_groupes <- c("#E69F00", "#56B4E9")

# Create output folder for plots if it doesn't exist
plots_dir <- "plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

message("Environment configured successfully. Output plots will be saved to: ", plots_dir)

# ==============================================================================
# Mathematical Utility Functions (Manifold Tools & Distances)
# ==============================================================================

# Complex HPD -> Real SPD Block Isomorphism
iso_real <- function(H) {
  rbind(cbind(Re(H), -Im(H)), cbind(Im(H), Re(H)))
}

# Real SPD -> Complex HPD Inverse Isomorphism
iso_complex <- function(M) {
  n <- nrow(M) / 2
  M11 <- M[1:n, 1:n]
  M21 <- M[(n + 1):(2 * n), 1:n]
  matrix(complex(real = M11, imaginary = M21), n, n)
}

# Extracted Symmetric Square Root (Real SPD)
sym_sqrt <- function(M, inv = FALSE) {
  e <- eigen(M, symmetric = TRUE)
  val <- pmax(e$values, 1e-12)
  if (inv) val <- 1.0 / sqrt(val) else val <- sqrt(val)
  e$vectors %*% diag(val) %*% t(e$vectors)
}

# 1. Distances
riemannian_distance <- function(H1, H2) {
  M1 <- iso_real(H1)
  M2 <- iso_real(H2)
  M1_inv_half <- sym_sqrt(M1, inv = TRUE)
  core <- M1_inv_half %*% M2 %*% M1_inv_half
  core <- (core + t(core)) / 2
  evals_core <- pmax(eigen(core, symmetric = TRUE, only.values = TRUE)$values, 1e-12)
  sqrt(sum(log(evals_core)^2)) / sqrt(2)
}

airm_distance <- riemannian_distance # Alias

lerm_distance <- function(H1, H2) {
  M1 <- iso_real(H1)
  M2 <- iso_real(H2)
  e1 <- eigen(M1, symmetric = TRUE)
  log_M1 <- e1$vectors %*% diag(log(pmax(e1$values, 1e-12))) %*% t(e1$vectors)
  e2 <- eigen(M2, symmetric = TRUE)
  log_M2 <- e2$vectors %*% diag(log(pmax(e2$values, 1e-12))) %*% t(e2$vectors)
  sqrt(sum((log_M1 - log_M2)^2)) / sqrt(2)
}

euclidean_distance <- function(H1, H2) {
  sqrt(sum(Mod(H1 - H2)^2))
}

# 2. Differential Geometry (Log & Exp Maps on Real SPD manifold)
logmap <- function(Base, X) {
  B_inv_half <- sym_sqrt(Base, inv = TRUE)
  core <- B_inv_half %*% X %*% B_inv_half
  core <- (core + t(core)) / 2
  e <- eigen(core, symmetric = TRUE)
  log_core <- e$vectors %*% diag(log(pmax(e$values, 1e-12))) %*% t(e$vectors)
  B_half <- sym_sqrt(Base, inv = FALSE)
  B_half %*% log_core %*% B_half
}

expmap <- function(Base, V) {
  B_inv_half <- sym_sqrt(Base, inv = TRUE)
  core <- B_inv_half %*% V %*% B_inv_half
  core <- (core + t(core)) / 2
  e <- eigen(core, symmetric = TRUE)
  exp_core <- e$vectors %*% diag(exp(e$values)) %*% t(e$vectors)
  B_half <- sym_sqrt(Base, inv = FALSE)
  B_half %*% exp_core %*% B_half
}

# 3. Fréchet Mean Solver (Optimized Gradient Descent)
frechet_mean_hpd <- function(hpd_list, max_iter = 15, tol = 1e-5) {
  n_mats <- length(hpd_list)
  spd_list <- lapply(hpd_list, iso_real)
  M <- Reduce("+", spd_list) / n_mats # Initialization with Arithmetic Mean

  for (iter in 1:max_iter) {
    # Precomputation of M square roots to optimize speed (massive 20x speedup)
    M_half <- sym_sqrt(M, inv = FALSE)
    M_inv_half <- sym_sqrt(M, inv = TRUE)
    
    log_sum <- matrix(0, nrow(M), ncol(M))
    for (i in 1:n_mats) {
      X <- spd_list[[i]]
      core <- M_inv_half %*% X %*% M_inv_half
      core <- (core + t(core)) / 2
      e <- eigen(core, symmetric = TRUE)
      log_core <- e$vectors %*% diag(log(pmax(e$values, 1e-12))) %*% t(e$vectors)
      log_sum <- log_sum + M_half %*% log_core %*% M_half
    }
    
    grad <- log_sum / n_mats
    step_norm <- sqrt(sum(diag(crossprod(grad))))
    if (step_norm < tol) break
    
    M <- expmap(M, grad)
    M <- (M + t(M)) / 2 # Enforce exact symmetry
  }
  iso_complex(M)
}

# 4. Riemannian K-Means
kmeans_riemann_hpd <- function(hpd_list, k = 2, max_iter = 50, nstart = 1) {
  n <- length(hpd_list)
  best_inertia <- Inf
  best_labels <- rep(0, n)
  best_centroids <- list()
  best_iter <- 0

  for (start in 1:nstart) {
    centroids <- hpd_list[sample(1:n, k)]
    labels <- rep(0, n)

    for (iter in 1:max_iter) {
      old_labels <- labels

      # Assignment Step
      for (i in 1:n) {
        dists <- sapply(1:k, function(c) airm_distance(hpd_list[[i]], centroids[[c]]))
        labels[i] <- which.min(dists)
      }

      if (all(labels == old_labels)) break

      # Update Step (Fréchet Mean of each cluster)
      for (c in 1:k) {
        cluster_pts <- hpd_list[labels == c]
        if (length(cluster_pts) > 0) {
          centroids[[c]] <- frechet_mean_hpd(cluster_pts)
        }
      }
    }

    # Geodesic within-cluster inertia
    inertia <- 0
    for (i in 1:n) {
      c_assigned <- labels[i]
      inertia <- inertia + airm_distance(hpd_list[[i]], centroids[[c_assigned]])^2
    }

    if (inertia < best_inertia) {
      best_inertia <- inertia
      best_labels <- labels
      best_centroids <- centroids
      best_iter <- iter
    }
  }

  list(cluster = best_labels, centers = best_centroids, iter = best_iter, inertia = best_inertia)
}

# ==============================================================================
# Experiment 1: Phase Separation under Swelling Scaling Traps
# ==============================================================================
message("Running Experiment 1...")

set.seed(12345)
n_subjects_per_group <- 20
total_subjects <- n_subjects_per_group * 2
n_timepoints <- 300
n_channels <- 3
freq_interest <- 0.1
spans <- c(5, 5)

dataset <- list()

for (subj in 1:total_subjects) {
  if (subj <= n_subjects_per_group) {
    group_label <- 1
    phase_shift <- c(0, pi / 2, pi)
  } else {
    group_label <- 2
    phase_shift <- c(0, -pi / 2, pi)
  }

  base_ampl <- 1.5
  ar_noise_range <- c(0.4, 0.6)

  x_multi <- matrix(0, nrow = n_timepoints, ncol = n_channels)
  for (i in 1:n_channels) {
    alpha_comp <- base_ampl * cos(2 * pi * freq_interest * (1:n_timepoints) - phase_shift[i])
    local_noise <- arima.sim(model = list(ar = runif(1, ar_noise_range[1], ar_noise_range[2])), n = n_timepoints)
    x_multi[, i] <- alpha_comp + local_noise
  }

  # Apply swelling scale distortion
  subject_scale <- runif(1, 0.5, 2.0)
  for (i in 1:n_channels) {
    x_multi[, i] <- x_multi[, i] * subject_scale
  }

  spect_multi <- mvspec(x_multi, spans = spans, taper = 0.1, log = "no", plot = FALSE)
  idx_freq <- which.min(abs(spect_multi$freq - freq_interest))
  H <- matrix(complex(real = 0, imaginary = 0), nrow = n_channels, ncol = n_channels)

  for (i in 1:n_channels) H[i, i] <- spect_multi$spec[idx_freq, i]
  for (i in 1:(n_channels - 1)) {
    for (j in (i + 1):n_channels) {
      v <- spect_multi$fxx[i, j, idx_freq]
      H[i, j] <- v
      H[j, i] <- Conj(v)
    }
  }

  H <- H + diag(0.5, n_channels)
  H <- (H + t(Conj(H))) / 2
  dataset[[subj]] <- list(matrix = H, group = group_label)

  if (subj == 1) {
    x_sample_g1 <- x_multi
    spec_sample_g1 <- spect_multi
  }
  if (subj == 21) {
    x_sample_g2 <- x_multi
    spec_sample_g2 <- spect_multi
  }
}

matrices <- lapply(dataset, function(x) x$matrix)
true_labels <- sapply(dataset, function(x) x$group)
N <- length(dataset)

# 1. Save Signal Plot
png(file.path(plots_dir, "exp1_signals.png"), width = 12, height = 8, units = "in", res = 150)
par(mfrow = c(2, 2))
ts.plot(x_sample_g1, col = 1:3, main = "Time Series: Patient 1 (Group 1)", ylab = "Amplitude")
legend("topright", legend = c("Ch1", "Ch2", "Ch3"), col = 1:3, lty = 1, cex = 0.8)
ts.plot(x_sample_g2, col = 1:3, main = "Time Series: Patient 21 (Group 2)", ylab = "Amplitude")
legend("topright", legend = c("Ch1", "Ch2", "Ch3"), col = 1:3, lty = 1, cex = 0.8)
plot(spec_sample_g1$freq, spec_sample_g1$spec[, 1], type = "l", col = 1, log = "y",
     main = "Daniell Smoothed Spectrum (Patient 1)", ylab = "Power (log)", xlab = "Frequency")
lines(spec_sample_g1$freq, spec_sample_g1$spec[, 2], col = 2)
lines(spec_sample_g1$freq, spec_sample_g1$spec[, 3], col = 3)
abline(v = 0.1, col = "blue", lty = 2, lwd = 2)
plot(spec_sample_g2$freq, spec_sample_g2$spec[, 1], type = "l", col = 1, log = "y",
     main = "Daniell Smoothed Spectrum (Patient 21)", ylab = "Power (log)", xlab = "Frequency")
lines(spec_sample_g2$freq, spec_sample_g2$spec[, 2], col = 2)
lines(spec_sample_g2$freq, spec_sample_g2$spec[, 3], col = 3)
abline(v = 0.1, col = "blue", lty = 2, lwd = 2)
dev.off()

# 2. Distance Computation
dist_euclid <- matrix(0, N, N)
dist_lerm <- matrix(0, N, N)
dist_airm <- matrix(0, N, N)

for (i in 1:(N - 1)) {
  for (j in (i + 1):N) {
    M1 <- matrices[[i]]
    M2 <- matrices[[j]]
    d_eu <- euclidean_distance(M1, M2)
    d_le <- lerm_distance(M1, M2)
    d_ai <- airm_distance(M1, M2)
    dist_euclid[i, j] <- dist_euclid[j, i] <- d_eu
    dist_lerm[i, j] <- dist_lerm[j, i] <- d_le
    dist_airm[i, j] <- dist_airm[j, i] <- d_ai
  }
}

# Save Heatmaps
png(file.path(plots_dir, "exp1_heatmaps.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
image(1:N, 1:N, dist_euclid[, N:1], col = hcl.colors(50, "Cividis", rev = TRUE), axes = FALSE,
      main = "Euclidean Distance Matrix N x N\n(No clear block structure)")
box()
abline(v = 20.5, h = 20.5, col = "white", lwd = 2, lty = 2)
image(1:N, 1:N, dist_airm[, N:1], col = hcl.colors(50, "Cividis", rev = TRUE), axes = FALSE,
      main = "Riemannian Distance Matrix N x N\n(Intra-group structure appears clearly)")
box()
abline(v = 20.5, h = 20.5, col = "white", lwd = 2, lty = 2)
dev.off()

# 3. Clustering Execution
flat_data <- t(sapply(matrices, function(M) c(Re(as.vector(M)), Im(as.vector(M)))))
set.seed(123)
res_euclid <- kmeans(flat_data, centers = 2, nstart = 10)

flat_lerm_data <- t(sapply(matrices, function(H) {
  M <- iso_real(H)
  e <- eigen(M, symmetric = TRUE)
  log_M <- e$vectors %*% diag(log(pmax(e$values, 1e-12))) %*% t(e$vectors)
  as.vector(log_M)
}))
set.seed(123)
res_lerm <- kmeans(flat_lerm_data, centers = 2, nstart = 10)

set.seed(42)
res_airm <- kmeans_riemann_hpd(matrices, k = 2)

acc_euclid <- adjustedRandIndex(true_labels, res_euclid$cluster)
acc_lerm <- adjustedRandIndex(true_labels, res_lerm$cluster)
acc_airm <- adjustedRandIndex(true_labels, res_airm$cluster)

# Robust Confusion Matrix Plot Function
plot_conf_matrix <- function(truth, pred, title, ari) {
  tbl <- table(True = factor(truth, levels = c(1, 2)), Predicted = factor(pred, levels = c(1, 2)))
  if (tbl[1, 1] + tbl[2, 2] < tbl[1, 2] + tbl[2, 1]) {
    tbl <- tbl[, c(2, 1)]
  }
  plot(c(-1, 2.5), c(-1, 2.5), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = sprintf("%s (ARI = %.3f)", title, ari))
  text(x = c(0, 1), y = rep(1.5, 2), labels = c("Pred 1", "Pred 2"), font = 2)
  text(x = rep(-0.5, 2), y = c(1, 0), labels = c("True 1", "True 2"), font = 2)
  for (i in 1:2) {
    for (j in 1:2) {
      rect(j - 1.5, i - 1.5, j - 0.5, i - 0.5, col = "white", border = "gray")
      text(j - 1, i - 1, as.character(tbl[3 - i, j]), cex = 2, col = ifelse(tbl[3 - i, j] > 0, "black", "gray80"))
    }
  }
}

png(file.path(plots_dir, "exp1_confusion_matrices.png"), width = 10, height = 4, units = "in", res = 150)
par(mfrow = c(1, 3))
plot_conf_matrix(true_labels, res_euclid$cluster, "Euclidean", acc_euclid)
plot_conf_matrix(true_labels, res_lerm$cluster, "LERM", acc_lerm)
plot_conf_matrix(true_labels, res_airm$cluster, "AIRM", acc_airm)
dev.off()

# Save MDS Projections
mds_euclid <- cmdscale(as.dist(dist_euclid), k = 2)
mds_lerm <- cmdscale(as.dist(dist_lerm), k = 2)
mds_airm <- cmdscale(as.dist(dist_airm), k = 2)

png(file.path(plots_dir, "exp1_mds.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2))
plot(mds_euclid, col = col_groupes[true_labels], pch = 19, cex = 1.5,
     main = "MDS Projection - Euclidean Space", xlab = "Dimension 1", ylab = "Dimension 2")
plot(mds_lerm, col = col_groupes[true_labels], pch = 19, cex = 1.5,
     main = "MDS Projection - LERM Manifold", xlab = "Dimension 1", ylab = "Dimension 2")
plot(mds_airm, col = col_groupes[true_labels], pch = 19, cex = 1.5,
     main = "MDS Projection - AIRM Manifold", xlab = "Dimension 1", ylab = "Dimension 2")
dev.off()

# Save ARI Performance Barplot
png(file.path(plots_dir, "exp1_ari_barplot.png"), width = 6, height = 5, units = "in", res = 150)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
b <- barplot(c(Euclidean = acc_euclid, LERM = acc_lerm, AIRM = acc_airm),
             main = "K-Means Score (Experiment 1)", ylab = "Adjusted Rand Index (ARI)",
             col = c("#e74c3c", "#f39c12", "#2ecc71"), ylim = c(0, 1.2))
text(b, c(acc_euclid, acc_lerm, acc_airm) + 0.05,
     labels = sprintf("%.3f", c(acc_euclid, acc_lerm, acc_airm)), font = 2, cex = 1.2)
dev.off()

# Print Table to Terminal
comp_table <- data.frame(
  Metric = c("Euclidean", "Log-Euclidean (LERM)", "Affine-Invariant (AIRM)"),
  ARI = c(acc_euclid, acc_lerm, acc_airm)
)
print(knitr::kable(comp_table, digits = 6, caption = "Comparison of Clustering Performance (ARI) - Experiment 1"))

# ==============================================================================
# Monte Carlo Helper Functions
# ==============================================================================

align_labels <- function(truth, pred) {
  if (sum(truth == pred) < sum(truth == (3 - pred))) {
    return(3 - pred)
  } else {
    return(pred)
  }
}

log_spd <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(log(pmax(e$values, 1e-12))) %*% t(e$vectors)
}

# Run one iteration of the Monte Carlo simulation
run_one_simulation <- function(T_val = 300, d_val = 3, n_subj_per_group = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  total_subj <- 2 * n_subj_per_group
  mats <- list()
  truth <- rep(1:2, each = n_subj_per_group)

  for (subj in 1:total_subj) {
    if (truth[subj] == 1) {
      phase_shift <- seq(0, pi, length.out = d_val)
    } else {
      phase_shift <- seq(0, -pi, length.out = d_val)
    }

    x_multi <- matrix(0, nrow = T_val, ncol = d_val)
    for (i in 1:d_val) {
      alpha_comp <- 1.5 * cos(2 * pi * 0.1 * (1:T_val) - phase_shift[i])
      local_noise <- arima.sim(model = list(ar = runif(1, 0.4, 0.6)), n = T_val)
      x_multi[, i] <- alpha_comp + local_noise
    }

    # Swelling Trap
    x_multi <- x_multi * runif(1, 0.5, 2.0)

    current_span <- max(3, floor(T_val / 20))
    if (current_span %% 2 == 0) current_span <- current_span + 1

    spect_multi <- mvspec(x_multi, spans = c(current_span, current_span), taper = 0.1, plot = FALSE, log = "no")
    idx_freq <- which.min(abs(spect_multi$freq - 0.1))
    
    H <- matrix(complex(real = 0, imaginary = 0), nrow = d_val, ncol = d_val)
    for (i in 1:d_val) H[i, i] <- spect_multi$spec[idx_freq, i]
    for (i in 1:(d_val - 1)) {
      for (j in (i + 1):d_val) {
        v <- spect_multi$fxx[i, j, idx_freq]
        H[i, j] <- v
        H[j, i] <- Conj(v)
      }
    }

    H <- H + diag(0.5, d_val)
    H <- (H + t(Conj(H))) / 2
    mats[[subj]] <- H
  }

  t_eu <- system.time({
    flat_eu <- t(sapply(mats, function(M) c(Re(as.vector(M)), Im(as.vector(M)))))
    res_eu <- kmeans(flat_eu, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_le <- system.time({
    flat_le <- t(sapply(mats, function(H) {
      M <- iso_real(H)
      as.vector(log_spd(M))
    }))
    res_le <- kmeans(flat_le, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_ai <- system.time({
    res_ai <- kmeans_riemann_hpd(mats, k = 2, nstart = 5)$cluster
  })["elapsed"]

  res_eu_aligned <- align_labels(truth, res_eu)
  res_le_aligned <- align_labels(truth, res_le)
  res_ai_aligned <- align_labels(truth, res_ai)

  list(
    ARI = c(Euclidean = adjustedRandIndex(truth, res_eu),
            LERM = adjustedRandIndex(truth, res_le),
            AIRM = adjustedRandIndex(truth, res_ai)),
    ACC = c(Euclidean = mean(truth == res_eu_aligned),
            LERM = mean(truth == res_le_aligned),
            AIRM = mean(truth == res_ai_aligned)),
    TIME = c(Euclidean = t_eu, LERM = t_le, AIRM = t_ai)
  )
}

# Monte Carlo Loop Wrapper
run_monte_carlo <- function(values, varying = c("T", "d", "N"), n_mc = 30, T_fixed = 300, d_fixed = 3, n_subj_fixed = 20) {
  varying <- match.arg(varying)
  metrics <- c("Euclidean", "LERM", "AIRM")

  mean_ari <- matrix(NA, nrow = length(values), ncol = 3, dimnames = list(NULL, metrics))
  sd_ari   <- matrix(NA, nrow = length(values), ncol = 3, dimnames = list(NULL, metrics))
  mean_acc <- matrix(NA, nrow = length(values), ncol = 3, dimnames = list(NULL, metrics))
  sd_acc   <- matrix(NA, nrow = length(values), ncol = 3, dimnames = list(NULL, metrics))
  mean_time <- matrix(NA, nrow = length(values), ncol = 3, dimnames = list(NULL, metrics))

  for (i in seq_along(values)) {
    ari_mc <- matrix(NA, nrow = n_mc, ncol = 3, dimnames = list(NULL, metrics))
    acc_mc <- matrix(NA, nrow = n_mc, ncol = 3, dimnames = list(NULL, metrics))
    time_mc <- matrix(NA, nrow = n_mc, ncol = 3, dimnames = list(NULL, metrics))

    for (mc in 1:n_mc) {
      if (varying == "T") {
        res <- run_one_simulation(T_val = values[i], d_val = d_fixed, n_subj_per_group = n_subj_fixed, seed = 1000 + i * 100 + mc)
      } else if (varying == "d") {
        res <- run_one_simulation(T_val = T_fixed, d_val = values[i], n_subj_per_group = n_subj_fixed, seed = 2000 + i * 100 + mc)
      } else if (varying == "N") {
        res <- run_one_simulation(T_val = T_fixed, d_val = d_fixed, n_subj_per_group = values[i], seed = 3000 + i * 100 + mc)
      }
      ari_mc[mc, ]  <- res$ARI
      acc_mc[mc, ]  <- res$ACC
      time_mc[mc, ] <- res$TIME
    }

    mean_ari[i, ]  <- colMeans(ari_mc)
    sd_ari[i, ]    <- apply(ari_mc, 2, sd)
    mean_acc[i, ]  <- colMeans(acc_mc)
    sd_acc[i, ]    <- apply(acc_mc, 2, sd)
    mean_time[i, ] <- colMeans(time_mc)

    message("Finished MC evaluation for parameter: ", varying, " = ", values[i])
  }

  list(values = values, mean_ari = mean_ari, sd_ari = sd_ari, mean_acc = mean_acc, sd_acc = sd_acc, mean_time = mean_time)
}

# Plot MC Results
plot_mc_results <- function(mc_res, xlab_title, main_ari, main_acc) {
  values <- mc_res$values
  cols <- c(AIRM = "#2ecc71", LERM = "#f39c12", Euclidean = "#e74c3c")

  par(mfrow = c(1, 2))
  plot(values, mc_res$mean_ari[, "AIRM"], type = "b", pch = 19, col = cols["AIRM"], lwd = 3, ylim = c(-0.1, 1.1),
       xlab = xlab_title, ylab = "Mean Adjusted Rand Index (ARI)", main = main_ari)
  lines(values, mc_res$mean_ari[, "LERM"], type = "b", pch = 19, col = cols["LERM"], lwd = 3)
  lines(values, mc_res$mean_ari[, "Euclidean"], type = "b", pch = 19, col = cols["Euclidean"], lwd = 3)
  grid(col = "gray80", lty = "dotted")
  legend("bottomright", legend = c("AIRM", "LERM", "Euclidean"), col = cols, lwd = 3, pch = 19)

  plot(values, mc_res$mean_acc[, "AIRM"], type = "b", pch = 19, col = cols["AIRM"], lwd = 3, ylim = c(0.4, 1.05),
       xlab = xlab_title, ylab = "Mean Accuracy", main = main_acc)
  lines(values, mc_res$mean_acc[, "LERM"], type = "b", pch = 19, col = cols["LERM"], lwd = 3)
  lines(values, mc_res$mean_acc[, "Euclidean"], type = "b", pch = 19, col = cols["Euclidean"], lwd = 3)
  grid(col = "gray80", lty = "dotted")
  legend("bottomright", legend = c("AIRM", "LERM", "Euclidean"), col = cols, lwd = 3, pch = 19)
}

# ==============================================================================
# Experiment 2: Effect of Sample Size (T)
# ==============================================================================
message("\nRunning Experiment 2: Effect of Signal Length (T)...")
T_values <- c(50, 100, 200, 300, 500, 1000)

mc_T <- run_monte_carlo(values = T_values, varying = "T", n_mc = N_MC_RUNS, d_fixed = 3, n_subj_fixed = 20)

png(file.path(plots_dir, "exp2_sample_size.png"), width = 12, height = 5, units = "in", res = 150)
plot_mc_results(mc_T, xlab_title = "Sample Size (T timepoints)", main_ari = "Mean ARI vs. Signal Length", main_acc = "Mean Accuracy vs. Signal Length")
dev.off()

time_table_T <- data.frame(
  `Signal Length (T)` = mc_T$values,
  Euclidean = mc_T$mean_time[, "Euclidean"],
  LERM = mc_T$mean_time[, "LERM"],
  AIRM = mc_T$mean_time[, "AIRM"],
  check.names = FALSE
)
print(knitr::kable(time_table_T, digits = 6, caption = "Average Clustering Computation Time (s) vs Signal Length (T)"))

# ==============================================================================
# Experiment 3: Effect of High Dimensionality (d)
# ==============================================================================
message("\nRunning Experiment 3: Effect of Channel Dimensionality (d)...")
d_values <- 3:14

mc_d <- run_monte_carlo(values = d_values, varying = "d", n_mc = N_MC_RUNS, T_fixed = 300, n_subj_fixed = 20)

png(file.path(plots_dir, "exp3_dimension.png"), width = 12, height = 5, units = "in", res = 150)
plot_mc_results(mc_d, xlab_title = "Number of Channels (d)", main_ari = "Mean ARI in High Dimensions", main_acc = "Mean Accuracy in High Dimensions")
dev.off()

time_table_d <- data.frame(
  `Dimension (d)` = mc_d$values,
  Euclidean = mc_d$mean_time[, "Euclidean"],
  LERM = mc_d$mean_time[, "LERM"],
  AIRM = mc_d$mean_time[, "AIRM"],
  check.names = FALSE
)
print(knitr::kable(time_table_d, digits = 6, caption = "Average Clustering Computation Time (s) vs Number of Channels (d)"))

# ==============================================================================
# Experiment 4: Effect of Cohort Size (N)
# ==============================================================================
message("\nRunning Experiment 4: Effect of Cohort Size (N)...")
subj_values <- c(10, 20, 30, 40, 50, 60, 80, 100)

mc_N <- run_monte_carlo(values = subj_values, varying = "N", n_mc = N_MC_RUNS, T_fixed = 300, d_fixed = 14)
mc_N$total_N <- subj_values * 2

png(file.path(plots_dir, "exp4_cohort.png"), width = 12, height = 5, units = "in", res = 150)
plot_mc_results(mc_N, xlab_title = "Subjects per Group", main_ari = "Mean ARI vs. Cohort Size (d = 14)", main_acc = "Mean Accuracy vs. Cohort Size (d = 14)")
dev.off()

time_table_N <- data.frame(
  `Subjects per Group` = mc_N$values,
  `Total Cohort Size (2N)` = mc_N$total_N,
  Euclidean = mc_N$mean_time[, "Euclidean"],
  LERM = mc_N$mean_time[, "LERM"],
  AIRM = mc_N$mean_time[, "AIRM"],
  check.names = FALSE
)
print(knitr::kable(time_table_N, digits = 6, caption = "Average Clustering Computation Time (s) vs Cohort Size (N)"))

# ==============================================================================
# Model A: Pure Autoregressive AR(2) Process Simulation
# ==============================================================================
message("\nRunning Robustness Study: Model A (Pure AR(2) process)...")

run_one_simulation_ar2 <- function(T_val = 300, d_val = 5, n_subj_per_group = 20, seed = NULL, save_mats = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  total_subj <- 2 * n_subj_per_group
  mats  <- list()
  truth <- rep(1:2, each = n_subj_per_group)

  r <- 0.90
  ar_params <- list(
    group1 = c(phi1 = 2 * r * cos(2 * pi * 0.1), phi2 = -r^2),
    group2 = c(phi1 = 2 * r * cos(2 * pi * 0.2), phi2 = -r^2)
  )

  for (subj in 1:total_subj) {
    g    <- truth[subj]
    phi1 <- ar_params[[g]]["phi1"]
    phi2 <- ar_params[[g]]["phi2"]

    x_multi <- matrix(0, nrow = T_val, ncol = d_val)
    for (c in 1:d_val) {
      x_multi[, c] <- arima.sim(model = list(ar = c(phi1, phi2)), n = T_val)
    }

    # Swelling scale distortion
    x_multi <- x_multi * runif(1, 0.5, 2.0)

    current_span <- max(3, floor(T_val / 20))
    if (current_span %% 2 == 0) current_span <- current_span + 1
    spect_multi <- mvspec(x_multi, spans = c(current_span, current_span), taper = 0.1, plot = FALSE, log = "no")

    idx_freq <- which.min(abs(spect_multi$freq - 0.1))
    H <- matrix(complex(real = 0, imaginary = 0), nrow = d_val, ncol = d_val)
    for (i in 1:d_val) H[i, i] <- spect_multi$spec[idx_freq, i]
    for (i in 1:(d_val - 1)) {
      for (j in (i + 1):d_val) {
        v       <- spect_multi$fxx[i, j, idx_freq]
        H[i, j] <- v
        H[j, i] <- Conj(v)
      }
    }
    H <- H + diag(0.5, d_val)
    H <- (H + t(Conj(H))) / 2
    mats[[subj]] <- H
  }

  t_eu <- system.time({
    flat_eu <- t(sapply(mats, function(M) c(Re(as.vector(M)), Im(as.vector(M)))))
    res_eu  <- kmeans(flat_eu, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_le <- system.time({
    flat_le <- t(sapply(mats, function(H) {
      M <- iso_real(H)
      as.vector(log_spd(M))
    }))
    res_le <- kmeans(flat_le, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_ai <- system.time({
    res_ai <- kmeans_riemann_hpd(mats, k = 2, nstart = 5)$cluster
  })["elapsed"]

  list(
    ARI  = c(Euclidean = adjustedRandIndex(truth, res_eu),
             LERM      = adjustedRandIndex(truth, res_le),
             AIRM      = adjustedRandIndex(truth, res_ai)),
    ACC  = c(Euclidean = mean(truth == align_labels(truth, res_eu)),
             LERM      = mean(truth == align_labels(truth, res_le)),
             AIRM      = mean(truth == align_labels(truth, res_ai))),
    TIME = c(Euclidean = t_eu, LERM = t_le, AIRM = t_ai),
    PRED = list(eu = res_eu, le = res_le, ai = res_ai),
    TRUTH = truth,
    MATS  = if (save_mats) mats else NULL
  )
}

# Monte Carlo evaluation
ari_mc_a  <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
acc_mc_a  <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
time_mc_a <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
last_run_a <- NULL

for (mc in 1:N_MC_RUNS) {
  res <- run_one_simulation_ar2(T_val = 300, d_val = 5, n_subj_per_group = 20, seed = 4000 + mc, save_mats = (mc == N_MC_RUNS))
  ari_mc_a[mc, ]  <- res$ARI
  acc_mc_a[mc, ]  <- res$ACC
  time_mc_a[mc, ] <- res$TIME
  if (mc == N_MC_RUNS) last_run_a <- res
}

summary_a <- data.frame(
  Metric          = c("Euclidean", "Log-Euclidean (LERM)", "Affine-Invariant (AIRM)"),
  `Mean ARI`      = colMeans(ari_mc_a),
  `SD ARI`        = apply(ari_mc_a, 2, sd),
  `Mean Accuracy` = colMeans(acc_mc_a),
  `SD Accuracy`   = apply(acc_mc_a, 2, sd),
  `Mean Time (s)` = colMeans(time_mc_a),
  check.names = FALSE
)
print(knitr::kable(summary_a, digits = 3, caption = "Model A — Pure AR(2): Mean +/- SD over Monte Carlo runs"))

# 1. Save Boxplot and Barplot
png(file.path(plots_dir, "model_a_ari_plots.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
boxplot(ari_mc_a, col = c("#e74c3c", "#f39c12", "#2ecc71"), main = "Model A — ARI Distribution", ylab = "Adjusted Rand Index (ARI)", ylim = c(-0.1, 1.15), names = c("Euclidean", "LERM", "AIRM"))
abline(h = 1.0, lty = 2, col = "gray40", lwd = 1.5)
abline(h = 0.0, lty = 3, col = "gray70")

means_a <- colMeans(ari_mc_a)
sds_a   <- apply(ari_mc_a, 2, sd)
bp <- barplot(means_a, col = c("#e74c3c", "#f39c12", "#2ecc71"), main = "Model A — Mean ARI with +/- 1 SD", ylab = "Mean ARI", ylim = c(0, 1.3), names.arg = c("Euclidean", "LERM", "AIRM"))
arrows(bp, means_a - sds_a, bp, means_a + sds_a, angle = 90, code = 3, length = 0.08, lwd = 2)
text(bp, means_a + sds_a + 0.07, labels = sprintf("%.3f\n+/-%.3f", means_a, sds_a), font = 2, cex = 0.95)
dev.off()

# 2. Save Confusion Matrices
png(file.path(plots_dir, "model_a_confusion_matrices.png"), width = 10, height = 4, units = "in", res = 150)
par(mfrow = c(1, 3))
plot_conf_matrix(last_run_a$TRUTH, last_run_a$PRED$eu, "Euclidean", last_run_a$ARI["Euclidean"])
plot_conf_matrix(last_run_a$TRUTH, last_run_a$PRED$le, "LERM",      last_run_a$ARI["LERM"])
plot_conf_matrix(last_run_a$TRUTH, last_run_a$PRED$ai, "AIRM",      last_run_a$ARI["AIRM"])
dev.off()

# 3. Save MDS Projections
mats_last_a <- last_run_a$MATS
truth_mds_a <- last_run_a$TRUTH
N_a <- length(mats_last_a)
dist_eu_a <- dist_le_a <- dist_ai_a <- matrix(0, N_a, N_a)
for (i in 1:(N_a - 1)) {
  for (j in (i + 1):N_a) {
    dist_eu_a[i, j] <- dist_eu_a[j, i] <- euclidean_distance(mats_last_a[[i]], mats_last_a[[j]])
    dist_le_a[i, j] <- dist_le_a[j, i] <- lerm_distance(mats_last_a[[i]], mats_last_a[[j]])
    dist_ai_a[i, j] <- dist_ai_a[j, i] <- airm_distance(mats_last_a[[i]], mats_last_a[[j]])
  }
}

png(file.path(plots_dir, "model_a_mds.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2))
plot(cmdscale(as.dist(dist_eu_a), k=2), col=col_groupes[truth_mds_a], pch=19, cex=1.5,
     main="MDS — Euclidean (Model A)", xlab="Dim 1", ylab="Dim 2")
plot(cmdscale(as.dist(dist_le_a), k=2), col=col_groupes[truth_mds_a], pch=19, cex=1.5,
     main="MDS — LERM (Model A)", xlab="Dim 1", ylab="Dim 2")
plot(cmdscale(as.dist(dist_ai_a), k=2), col=col_groupes[truth_mds_a], pch=19, cex=1.5,
     main="MDS — AIRM (Model A)", xlab="Dim 1", ylab="Dim 2")
legend("topright", legend=c("Group 1 (Alpha)","Group 2 (Beta)"), col=col_groupes, pch=19, cex=0.8)
dev.off()

# ==============================================================================
# Model B: Spatially and Temporally Correlated Noise Simulation
# ==============================================================================
message("\nRunning Robustness Study: Model B (VAR(1) Spatiotemporal Noise)...")

run_one_simulation_var <- function(T_val = 300, d_val = 5, n_subj_per_group = 20, phi_ar = 0.5, rho_spatial = 0.5, seed = NULL, save_mats = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(abs(phi_ar) < 1)

  total_subj <- 2 * n_subj_per_group
  mats  <- list()
  truth <- rep(1:2, each = n_subj_per_group)

  Sigma      <- outer(1:d_val, 1:d_val, function(i, j) rho_spatial^abs(i - j))
  Sigma_chol <- chol(Sigma)

  for (subj in 1:total_subj) {
    if (truth[subj] == 1) {
      phase_shift <- seq(0,  pi, length.out = d_val)
    } else {
      phase_shift <- seq(0, -pi, length.out = d_val)
    }

    # VAR(1) generation
    W <- matrix(0, nrow = T_val, ncol = d_val)
    for (t in 2:T_val) {
      z_t    <- rnorm(d_val)
      E_t    <- as.numeric(t(Sigma_chol) %*% z_t)
      W[t, ] <- phi_ar * W[t - 1, ] + E_t
    }

    # Signal = harmonic + VAR(1) noise
    x_multi <- matrix(0, nrow = T_val, ncol = d_val)
    for (c in 1:d_val) {
      x_multi[, c] <- 1.5 * cos(2 * pi * 0.1 * (1:T_val) - phase_shift[c]) + W[, c]
    }

    # Swelling Trap
    x_multi <- x_multi * runif(1, 0.5, 2.0)

    current_span <- max(3, floor(T_val / 20))
    if (current_span %% 2 == 0) current_span <- current_span + 1
    spect_multi <- mvspec(x_multi, spans = c(current_span, current_span), taper = 0.1, plot = FALSE, log = "no")

    idx_freq <- which.min(abs(spect_multi$freq - 0.1))
    H <- matrix(complex(real = 0, imaginary = 0), nrow = d_val, ncol = d_val)
    for (i in 1:d_val) H[i, i] <- spect_multi$spec[idx_freq, i]
    for (i in 1:(d_val - 1)) {
      for (j in (i + 1):d_val) {
        v       <- spect_multi$fxx[i, j, idx_freq]
        H[i, j] <- v
        H[j, i] <- Conj(v)
      }
    }
    H <- H + diag(0.5, d_val)
    H <- (H + t(Conj(H))) / 2
    mats[[subj]] <- H
  }

  t_eu <- system.time({
    flat_eu <- t(sapply(mats, function(M) c(Re(as.vector(M)), Im(as.vector(M)))))
    res_eu  <- kmeans(flat_eu, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_le <- system.time({
    flat_le <- t(sapply(mats, function(H) {
      M <- iso_real(H)
      as.vector(log_spd(M))
    }))
    res_le <- kmeans(flat_le, centers = 2, nstart = 20)$cluster
  })["elapsed"]

  t_ai <- system.time({
    res_ai <- kmeans_riemann_hpd(mats, k = 2, nstart = 5)$cluster
  })["elapsed"]

  list(
    ARI   = c(Euclidean = adjustedRandIndex(truth, res_eu),
              LERM      = adjustedRandIndex(truth, res_le),
              AIRM      = adjustedRandIndex(truth, res_ai)),
    ACC   = c(Euclidean = mean(truth == align_labels(truth, res_eu)),
              LERM      = mean(truth == align_labels(truth, res_le)),
              AIRM      = mean(truth == align_labels(truth, res_ai))),
    TIME  = c(Euclidean = t_eu, LERM = t_le, AIRM = t_ai),
    PRED  = list(eu = res_eu, le = res_le, ai = res_ai),
    TRUTH = truth,
    MATS  = if (save_mats) mats else NULL
  )
}

# Monte Carlo evaluation
ari_mc_b  <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
acc_mc_b  <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
time_mc_b <- matrix(0, N_MC_RUNS, 3, dimnames = list(NULL, c("Euclidean","LERM","AIRM")))
last_run_b <- NULL

for (mc in 1:N_MC_RUNS) {
  res <- run_one_simulation_var(T_val = 300, d_val = 5, n_subj_per_group = 20, phi_ar = 0.5, rho_spatial = 0.5, seed = 5000 + mc, save_mats = (mc == N_MC_RUNS))
  ari_mc_b[mc, ]  <- res$ARI
  acc_mc_b[mc, ]  <- res$ACC
  time_mc_b[mc, ] <- res$TIME
  if (mc == N_MC_RUNS) last_run_b <- res
}

summary_b <- data.frame(
  Metric          = c("Euclidean", "Log-Euclidean (LERM)", "Affine-Invariant (AIRM)"),
  `Mean ARI`      = colMeans(ari_mc_b),
  `SD ARI`        = apply(ari_mc_b, 2, sd),
  `Mean Accuracy` = colMeans(acc_mc_b),
  `SD Accuracy`   = apply(acc_mc_b, 2, sd),
  `Mean Time (s)` = colMeans(time_mc_b),
  check.names = FALSE
)
print(knitr::kable(summary_b, digits = 3, caption = "Model B — Harmonic + VAR(1) Noise: Mean +/- SD over Monte Carlo runs"))

# 1. Save Boxplot and Barplot
png(file.path(plots_dir, "model_b_ari_plots.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
boxplot(ari_mc_b, col = c("#e74c3c", "#f39c12", "#2ecc71"), main = "Model B — ARI Distribution", ylab = "Adjusted Rand Index (ARI)", ylim = c(-0.1, 1.15), names = c("Euclidean", "LERM", "AIRM"))
abline(h = 1.0, lty = 2, col = "gray40", lwd = 1.5)
abline(h = 0.0, lty = 3, col = "gray70")

means_b <- colMeans(ari_mc_b)
sds_b   <- apply(ari_mc_b, 2, sd)
bp <- barplot(means_b, col = c("#e74c3c", "#f39c12", "#2ecc71"), main = "Model B — Mean ARI with +/- 1 SD", ylab = "Mean ARI", ylim = c(0, 1.3), names.arg = c("Euclidean", "LERM", "AIRM"))
arrows(bp, means_b - sds_b, bp, means_b + sds_b, angle = 90, code = 3, length = 0.08, lwd = 2)
text(bp, means_b + sds_b + 0.07, labels = sprintf("%.3f\n+/-%.3f", means_b, sds_b), font = 2, cex = 0.95)
dev.off()

# 2. Save Confusion Matrices
png(file.path(plots_dir, "model_b_confusion_matrices.png"), width = 10, height = 4, units = "in", res = 150)
par(mfrow = c(1, 3))
plot_conf_matrix(last_run_b$TRUTH, last_run_b$PRED$eu, "Euclidean", last_run_b$ARI["Euclidean"])
plot_conf_matrix(last_run_b$TRUTH, last_run_b$PRED$le, "LERM",      last_run_b$ARI["LERM"])
plot_conf_matrix(last_run_b$TRUTH, last_run_b$PRED$ai, "AIRM",      last_run_b$ARI["AIRM"])
dev.off()

# 3. Save MDS Projections
mats_b  <- last_run_b$MATS
truth_b <- last_run_b$TRUTH
N_b     <- length(mats_b)
dist_eu_b <- dist_le_b <- dist_ai_b <- matrix(0, N_b, N_b)
for (i in 1:(N_b - 1)) {
  for (j in (i + 1):N_b) {
    dist_eu_b[i, j] <- dist_eu_b[j, i] <- euclidean_distance(mats_b[[i]], mats_b[[j]])
    dist_le_b[i, j] <- dist_le_b[j, i] <- lerm_distance(mats_b[[i]], mats_b[[j]])
    dist_ai_b[i, j] <- dist_ai_b[j, i] <- airm_distance(mats_b[[i]], mats_b[[j]])
  }
}

png(file.path(plots_dir, "model_b_mds.png"), width = 12, height = 5, units = "in", res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2))
plot(cmdscale(as.dist(dist_eu_b), k = 2), col = col_groupes[truth_b], pch = 19, cex = 1.5,
     main = "MDS — Euclidean (Model B)", xlab = "Dim 1", ylab = "Dim 2")
plot(cmdscale(as.dist(dist_le_b), k = 2), col = col_groupes[truth_b], pch = 19, cex = 1.5,
     main = "MDS — LERM (Model B)", xlab = "Dim 1", ylab = "Dim 2")
plot(cmdscale(as.dist(dist_ai_b), k = 2), col = col_groupes[truth_b], pch = 19, cex = 1.5,
     main = "MDS — AIRM (Model B)", xlab = "Dim 1", ylab = "Dim 2")
legend("topright", legend = c("Group 1 (+phase)", "Group 2 (-phase)"), col = col_groupes, pch = 19, cex = 0.8)
dev.off()

message("\nAll simulations completed. Output plots saved in: ", plots_dir)
