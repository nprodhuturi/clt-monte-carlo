simulate_means <- function(generator, n, reps = B) {
  replicate(reps, mean(generator(n)))
}

skew <- function(x) {
  m3 <- mean((x - mean(x))^3)
  m2 <- mean((x - mean(x))^2)
  m3 / m2^(3/2)
}

kurt <- function(x) {
  m4 <- mean((x - mean(x))^4)
  m2 <- mean((x - mean(x))^2)
  m4 / m2^2 - 3
}

simulate_all_n <- function(generator, n_values) {
  all_means <- matrix(NA, nrow = B, ncol = length(n_values))
  colnames(all_means) <- n_values
  for (i in seq_along(n_values)) {
    all_means[, i] <- simulate_means(generator, n_values[i])
  }
  all_means
}

summarize_means <- function(all_means, mu, sigma) {
  n_vec <- as.numeric(colnames(all_means))
  sim_sd <- apply(all_means, 2, sd)
  data.frame(
    n = n_vec,
    sim_mean = apply(all_means, 2, mean),
    theory_mean = mu,
    sim_sd = sim_sd,
    theory_se = sigma / sqrt(n_vec),
    sd_ratio = sim_sd / (sigma / sqrt(n_vec)),
    skewness = apply(all_means, 2, skew),
    excess_kurtosis = apply(all_means, 2, kurt),
    row.names = NULL
  )
}

plot_sampling_distributions <- function(generator, n_values, pop_mean,
                                        pop_sd, label) {
  old_par <- par(mfrow = c(1, length(n_values)))
  on.exit(par(old_par))
  
  for (n in n_values) {
    means <- simulate_means(generator, n)
    hist(means, probability = TRUE, breaks = 35,
         col = "gray85", border = "white",
         main = paste(label, "\nn =", n), xlab = "Sample mean")
    curve(dnorm(x, mean = pop_mean, sd = pop_sd / sqrt(n)),
          add = TRUE, col = "firebrick", lwd = 3)
  }
}

plot_histograms <- function(all_means, n_show, pop_mean, pop_sd, label, step = NULL) {
  old_par <- par(mfrow = c(1, length(n_show)))
  on.exit(par(old_par))
  
  for (n in n_show) {
    means <- all_means[, as.character(n)]
    if (is.null(step)) {
      breaks <- 35
    } else {
      w <- step / n
      k <- max(1, round(diff(range(means)) / w / 35))
      breaks <- seq(min(means) - w / 2, max(means) + k * w, by = k * w)
    }
    h <- hist(means, breaks = breaks, plot = FALSE)
    top <- max(h$density, dnorm(pop_mean, pop_mean, pop_sd / sqrt(n)), na.rm = TRUE)
    plot(h, freq = FALSE, ylim = c(0, top * 1.05),
         col = "gray85", border = "white",
         main = paste(label, "\nn =", n), xlab = "Sample mean")
    curve(dnorm(x, mean = pop_mean, sd = pop_sd / sqrt(n)),
          add = TRUE, col = "firebrick", lwd = 3)
  }
}

plot_qq <- function(all_means, n_show, label) {
  old_par <- par(mfrow = c(1, length(n_show)))
  on.exit(par(old_par))
  
  for (n in n_show) {
    means <- all_means[, as.character(n)]
    qqnorm(means, main = paste(label, "Q-Q plot\nn =", n))
    qqline(means, col = "firebrick", lwd = 2)
  }
}

flag_good <- function(results, skew_max = 0.25, kurt_max = 0.5, sd_tol = 0.05) {
  ok <- abs(results$skewness) <= skew_max &
    abs(results$excess_kurtosis) <= kurt_max &
    abs(results$sd_ratio - 1) <= sd_tol
  ok[is.na(ok)] <- FALSE
  ok
}

smallest_good_n <- function(results) {
  bad <- which(!results$good)
  if (length(bad) == 0) return(results$n[1])
  if (max(bad) == nrow(results)) return(NA)
  results$n[max(bad) + 1]
}