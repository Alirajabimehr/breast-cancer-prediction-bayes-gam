
# ======================================================================
# Breast Cancer Prediction: Bayesian Logistic Regression & GAMs (R)
# Author: Ali Rajabimehr
# Supervisors: Victor Miranda, Patricio Andres Maturana Russel
# Teammate: Craig Manning
# ======================================================================

# --------------------------- 0) Setup ---------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(corrplot)
  library(car)         # VIF
  library(rjags)       # JAGS interface
  library(coda)        # MCMC diagnostics + HPD intervals
  library(mgcv)        # GAM
  library(gridExtra)
  library(knitr)
  library(kableExtra)
})

set.seed(2024)

# ---- Paths (edit for your environment) ----
DATA_PATH    <- "data/data.csv"          # <- put your WDBC file here
RESULTS_DIR  <- "results"
FIG_DIR      <- file.path(RESULTS_DIR, "figures")
dir.create(RESULTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,      showWarnings = FALSE, recursive = TRUE)

# ----------------------- 1) Load & Prepare Data -----------------------
# Expected: a CSV with diagnosis in {"M","B"} and 30 features like radius_mean, etc.
df_raw <- read_csv(DATA_PATH, show_col_types = FALSE)

# Keep the 31 analytics columns (diagnosis + 30 features) if needed:
# If your CSV already has the 31 columns in order, skip select().
df <- df_raw %>%
  mutate(diagnosis = ifelse(diagnosis == "M", 1L, 0L))

# Quick sanity check
stopifnot(all(df$diagnosis %in% c(0L, 1L)))

# ----------------------- 2) Feature Selection -------------------------
# 2.1 Correlation-based filter vs target (abs corr > 0.5)
cor_mat <- suppressWarnings(cor(df %>% mutate(diagnosis = as.numeric(diagnosis))))
target_corr <- abs(cor_mat[, "diagnosis"])
corr_feats <- names(target_corr[target_corr > 0.5 & target_corr < 1])

# 2.2 ANOVA (keep p < 0.05)
anova_p <- sapply(df %>% select(-diagnosis), function(x) {
  summary(aov(x ~ df$diagnosis))[[1]][["Pr(>F)"]][1]
})
anova_feats <- names(anova_p[anova_p < 0.05])

# 2.3 Intersect
selected_feats <- intersect(corr_feats, anova_feats)

# 2.4 Remove high VIF iteratively (threshold = 10)
remove_high_vif <- function(data, threshold = 10) {
  # Build a linear model with diagnosis as numeric to compute VIFs
  v <- vif(lm(diagnosis ~ ., data = data))
  while (max(v) > threshold) {
    drop_name <- names(v)[which.max(v)]
    data <- data %>% select(-all_of(drop_name))
    v <- vif(lm(diagnosis ~ ., data = data))
  }
  data
}

df_model <- df %>%
  select(diagnosis, all_of(selected_feats)) %>%
  { remove_high_vif(.) }

# 2.5 Scale predictors
df_scaled <- df_model %>%
  mutate(across(-diagnosis, ~ scale(.) %>% as.vector))

message("Final predictors: ", paste(colnames(df_scaled)[-1], collapse = ", "))

# Optional: correlation heatmap (saved)
png(file.path(FIG_DIR, "correlation_heatmap.png"), width = 1200, height = 900)
corrplot(cor(df %>% select(-diagnosis)), method = "color", type = "upper", tl.cex = 0.6)
dev.off()

# ------------------ 3) Bayesian Logistic Regression -------------------
# Data for JAGS
X <- as.matrix(df_scaled %>% select(-diagnosis))
y <- df_scaled$diagnosis
N <- nrow(X); K <- ncol(X)

# JAGS model template
jags_template <- function(prec) {
  sprintf('
  model {
    for (i in 1:N) {
      y[i] ~ dbern(p[i])
      logit(p[i]) <- b0 + inprod(b[1:K], X[i,])
    }
    b0 ~ dnorm(0, %f)
    for (j in 1:K) {
      b[j] ~ dnorm(0, %f)
    }
  }', prec, prec)
}

# Prior sets
priors <- list(
  non_informative = 0.01,  # N(0, 0.01)
  informative_1   = 1,     # N(0, 1)
  informative_2   = 2,     # N(0, 2)
  informative_3   = 0.1,   # N(0, 0.1)  -> best in study
  informative_4   = 0.5    # N(0, 0.5)
)

run_jags <- function(prec, n_chains=4, n_adapt=1000, n_iter=5000) {
  model_string <- jags_template(prec)
  # write to a temp file for rjags
  tf <- tempfile(fileext = ".jags")
  writeLines(model_string, tf)
  jmod <- jags.model(
    file = tf,
    data = list(N=N, K=K, X=X, y=y),
    inits = function() list(b0 = 0, b = rep(0, K)),
    n.chains = n_chains,
    n.adapt = n_adapt
  )
  update(jmod, n.iter = n_adapt)
  coda.samples(jmod, variable.names = c("b0", paste0("b[", 1:K, "]")), n.iter = n_iter)
}

# Fit all priors and compute BIC from posterior means
log_lik_fn <- function(b0, b, X, y) {
  eta <- as.numeric(b0 + X %*% b)
  p <- 1 / (1 + exp(-eta))
  sum(dbinom(y, size = 1, prob = p, log = TRUE))
}

fits <- lapply(priors, run_jags)
names(fits) <- names(priors)

bic_values <- sapply(names(fits), function(nm) {
  samp <- fits[[nm]]
  pm   <- colMeans(as.matrix(samp))
  b0   <- pm["b0"]
  bvec <- pm[grep("^b\\[", names(pm))]
  ll   <- log_lik_fn(b0, bvec, X, y)
  -2 * ll + (K + 1) * log(N)  # BIC
})

best_name <- names(which.min(bic_values))
best_samples <- fits[[best_name]]

message(sprintf("Best Bayesian model by BIC: %s  (BIC = %.3f)", best_name, min(bic_values)))

# Convergence diagnostics (saved)
gd <- gelman.diag(best_samples, multivariate = FALSE)
eff <- effectiveSize(best_samples)
sink(file.path(RESULTS_DIR, "bayes_convergence.txt"))
cat("Gelman diagnostics:\n"); print(gd)
cat("\nEffective sample sizes:\n"); print(eff)
sink()

# Credible intervals (HPD)
hpd <- HPDinterval(best_samples, prob = 0.95)

# Extract posterior means for prediction
pm_best <- colMeans(as.matrix(best_samples))
b0_best <- pm_best["b0"]
b_best  <- pm_best[grep("^b\\[", names(pm_best))]

# Predictions & confusion matrix
p_hat <- 1 / (1 + exp(-(b0_best + as.numeric(X %*% b_best))))
pred  <- ifelse(p_hat > 0.5, 1L, 0L)

conf_bayes <- table(Predicted = pred, Actual = y)
acc_bayes  <- sum(diag(conf_bayes)) / sum(conf_bayes)
sens_bayes <- conf_bayes["1","1"] / sum(conf_bayes[,"1"])
spec_bayes <- conf_bayes["0","0"] / sum(conf_bayes[,"0"])

# Save metrics
bayes_metrics <- data.frame(
  Model = "Bayesian",
  Accuracy = acc_bayes,
  Sensitivity = sens_bayes,
  Specificity = spec_bayes,
  Prior = best_name,
  BIC = min(bic_values)
)
write.csv(bayes_metrics, file.path(RESULTS_DIR, "bayesian_metrics.csv"), row.names = FALSE)

# -------------------------- 4) GAM Model -------------------------------
gam_formula <- as.formula(
  paste("diagnosis ~", paste(sprintf("s(%s)", colnames(df_scaled)[-1]), collapse = " + "))
)
gam_fit <- gam(gam_formula, data = df_scaled, family = binomial)

gam_prob <- predict(gam_fit, type = "response")
gam_pred <- ifelse(gam_prob > 0.5, 1L, 0L)
conf_gam <- table(Predicted = gam_pred, Actual = y)
acc_gam  <- sum(diag(conf_gam)) / sum(conf_gam)
sens_gam <- conf_gam["1","1"] / sum(conf_gam[,"1"])
spec_gam <- conf_gam["0","0"] / sum(conf_gam[,"0"])

gam_metrics <- data.frame(
  Model = "GAM",
  Accuracy = acc_gam,
  Sensitivity = sens_gam,
  Specificity = spec_gam
)
write.csv(gam_metrics, file.path(RESULTS_DIR, "gam_metrics.csv"), row.names = FALSE)

# ----------------------- 5) Comparison Table ---------------------------
comp <- rbind(
  bayes_metrics[, c("Model","Accuracy","Sensitivity","Specificity")],
  gam_metrics
)
write.csv(comp, file.path(RESULTS_DIR, "model_comparison.csv"), row.names = FALSE)

# Pretty HTML table (optional)
html_tbl <- kable(comp, format = "html", caption = "Model Comparison") %>%
  kable_styling(bootstrap_options = c("striped","hover","condensed"), full_width = FALSE)
save_html <- file.path(RESULTS_DIR, "model_comparison.html")
cat(as.character(html_tbl), file = save_html)

# ----------------------- 6) Helpful Plots ------------------------------
# Example: plot smooths for first 6 terms (if many predictors)
n_show <- min(6, ncol(df_scaled) - 1)
plots <- vector("list", n_show)
for (i in seq_len(n_show)) {
  pd <- plot(gam_fit, se = TRUE, select = i, plot = FALSE)
  dfp <- data.frame(x = pd[[1]]$x, fit = pd[[1]]$fit, se = pd[[1]]$se)
  plots[[i]] <- ggplot(dfp, aes(x = x, y = fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.2) +
    labs(title = sprintf("GAM smooth: %s", colnames(df_scaled)[i+1]),
         x = colnames(df_scaled)[i+1], y = "Effect") +
    theme_minimal()
}
png(file.path(FIG_DIR, "gam_smooths.png"), width = 1600, height = 1200)
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

# ----------------------- 7) Console Summary ----------------------------
message("\n=== Performance Summary ===")
print(comp)
message("\nConfusion Matrix (Bayesian):"); print(conf_bayes)
message("\nConfusion Matrix (GAM):"); print(conf_gam)

# Done. All outputs saved under /results
