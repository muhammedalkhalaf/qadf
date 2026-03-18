#' Quantile Augmented Dickey-Fuller Unit Root Test
#'
#' Performs the Quantile ADF unit root test following Koenker & Xiao (2004).
#' Tests for unit root across the conditional distribution using quantile regression.
#'
#' @param y A numeric vector or time series.
#' @param model Character. Deterministic component: "c" for constant (default),
#'   "ct" for constant and trend, "nc" for no constant.
#' @param pmax Integer. Maximum number of lags for augmentation. Default is 4.
#' @param ic Character. Information criterion for lag selection: "AIC", "BIC", or "t-stat".
#' @param tau Numeric vector. Quantiles to test. Default is seq(0.1, 0.9, 0.1).
#'
#' @return An object of class "qadf" containing:
#' \itemize{
#'   \item \code{results} - Data frame with test results for each quantile
#'   \item \code{qks} - QKS statistic (supremum of absolute t-statistics)
#'   \item \code{qks_pvalue} - Bootstrap p-value for QKS (if computed)
#'   \item \code{model} - Model specification
#'   \item \code{lags} - Selected lag order
#'   \item \code{n} - Sample size
#'   \item \code{y} - Original series
#' }
#'
#' @details
#' The Quantile ADF test extends the standard ADF test to different quantiles
#' of the conditional distribution. This allows testing for unit root behavior
#' that may vary across the distribution (e.g., different persistence in 
#' expansions vs recessions).
#'
#' The test statistic at quantile tau is:
#' \deqn{t_n(\tau) = f_z \sqrt{y'_{-1} P_X y_{-1}} (\hat{\rho}(\tau) - 1)}
#'
#' Critical values are interpolated from Hansen (1995) response surfaces.
#'
#' The QKS (Quantile Kolmogorov-Smirnov) statistic is:
#' \deqn{QKS = \sup_{\tau} |t_n(\tau)|}
#'
#' @references
#' Koenker, R., & Xiao, Z. (2004). Unit root quantile autoregression inference.
#' \emph{Journal of the American Statistical Association}, 99(467), 775-787.
#'
#' @examples
#' # Generate a random walk
#' set.seed(123)
#' y <- cumsum(rnorm(200))
#'
#' # Basic QADF test
#' result <- qadf(y)
#' print(result)
#' summary(result)
#'
#' # Test specific quantiles
#' result <- qadf(y, tau = c(0.1, 0.5, 0.9))
#'
#' # With trend
#' result <- qadf(y, model = "ct")
#'
#' @export
#' @importFrom quantreg rq
#' @importFrom stats lm residuals coef pnorm dnorm AIC BIC
qadf <- function(y, model = c("c", "ct", "nc"), pmax = 4, 
                 ic = c("AIC", "BIC", "t-stat"), 
                 tau = seq(0.1, 0.9, 0.1)) {
  

  model <- match.arg(model)
  ic <- match.arg(ic)
  
  # Convert to numeric vector
  y <- as.numeric(y)
  n <- length(y)
  
  if (n < 20) {
    stop("Series too short. Need at least 20 observations.")
  }
  
  # First difference
  dy <- diff(y)
  y_lag1 <- y[-n]  # y_{t-1}
  
  # Select optimal lag using information criterion
  p_opt <- select_lag_qadf(dy, pmax, ic, model)
  
  # Prepare data matrices
  data_list <- prepare_qadf_data(y, dy, y_lag1, p_opt, model)
  
  # Run QADF for each quantile
  results <- lapply(tau, function(t) {
    fit_qadf_quantile(data_list, t, model, n - p_opt - 1)
  })
  
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Calculate QKS statistic
  qks <- max(abs(results_df$t_stat))
  
  # Structure output
  out <- list(
    results = results_df,
    qks = qks,
    qks_cv = qadf_qks_critical_values(),
    model = model,
    lags = p_opt,
    n = n,
    ic = ic,
    y = y,
    tau = tau,
    call = match.call()
  )
  
  class(out) <- "qadf"
  return(out)
}


#' @keywords internal
select_lag_qadf <- function(dy, pmax, ic, model) {
  n <- length(dy)
  
  if (pmax == 0) return(0)
  
  ic_values <- numeric(pmax + 1)
  
  for (p in 0:pmax) {
    if (p == 0) {
      if (model == "nc") {
        X <- NULL
        fit <- lm(dy ~ 0)
      } else if (model == "c") {
        fit <- lm(dy ~ 1)
      } else {
        fit <- lm(dy ~ 1 + seq_along(dy))
      }
    } else {
      # Create lagged differences
      dy_lags <- embed(dy, p + 1)
      y_dep <- dy_lags[, 1]
      X_lags <- dy_lags[, -1, drop = FALSE]
      colnames(X_lags) <- paste0("lag", 1:p)
      
      if (model == "nc") {
        fit <- lm(y_dep ~ X_lags - 1)
      } else if (model == "c") {
        fit <- lm(y_dep ~ X_lags)
      } else {
        trend <- seq_len(nrow(X_lags))
        fit <- lm(y_dep ~ X_lags + trend)
      }
    }
    
    if (ic == "AIC") {
      ic_values[p + 1] <- AIC(fit)
    } else if (ic == "BIC") {
      ic_values[p + 1] <- BIC(fit)
    } else {
      # t-stat criterion: use last lag significance
      if (p > 0) {
        coefs <- summary(fit)$coefficients
        last_lag_row <- grep("lag", rownames(coefs))
        if (length(last_lag_row) > 0) {
          last_t <- abs(coefs[max(last_lag_row), "t value"])
          ic_values[p + 1] <- -last_t  # Negative so min gives max t
        } else {
          ic_values[p + 1] <- Inf
        }
      } else {
        ic_values[p + 1] <- Inf
      }
    }
  }
  
  if (ic == "t-stat") {
    # Find last significant lag (|t| > 1.96)
    for (p in pmax:1) {
      if (-ic_values[p + 1] > 1.96) return(p)
    }
    return(0)
  }
  
  return(which.min(ic_values) - 1)
}


#' @keywords internal
prepare_qadf_data <- function(y, dy, y_lag1, p, model) {
  n <- length(dy)
  
  if (p > 0) {
    # Embed creates matrix with current and lagged values
    dy_embed <- embed(dy, p + 1)
    y_dep <- dy_embed[, 1]
    dy_lags <- dy_embed[, -1, drop = FALSE]
    
    # Align y_lag1
    y_lag1_aligned <- y_lag1[(p + 1):n]
    
    # Build design matrix
    if (model == "nc") {
      X <- cbind(y_lag1_aligned, dy_lags)
    } else if (model == "c") {
      X <- cbind(1, y_lag1_aligned, dy_lags)
    } else {  # ct
      trend <- seq_along(y_dep)
      X <- cbind(1, trend, y_lag1_aligned, dy_lags)
    }
  } else {
    y_dep <- dy
    y_lag1_aligned <- y_lag1
    
    if (model == "nc") {
      X <- matrix(y_lag1_aligned, ncol = 1)
    } else if (model == "c") {
      X <- cbind(1, y_lag1_aligned)
    } else {  # ct
      trend <- seq_along(y_dep)
      X <- cbind(1, trend, y_lag1_aligned)
    }
  }
  
  colnames(X) <- NULL  # Clean column names for quantreg
  
  list(
    y = y_dep,
    X = X,
    y_lag1 = y_lag1_aligned,
    p = p,
    model = model
  )
}


#' @keywords internal
fit_qadf_quantile <- function(data_list, tau, model, n) {
  
  y <- data_list$y
  X <- data_list$X
  y_lag1 <- data_list$y_lag1
  p <- data_list$p
  
  # Quantile regression
  qr_fit <- quantreg::rq(y ~ X - 1, tau = tau)
  
  # OLS regression for comparison
  ols_fit <- lm(y ~ X - 1)
  
  # Extract rho coefficient (position depends on model)
  if (model == "nc") {
    rho_idx <- 1
  } else if (model == "c") {
    rho_idx <- 2
  } else {  # ct
    rho_idx <- 3
  }
  
  rho_tau <- coef(qr_fit)[rho_idx]
  rho_ols <- coef(ols_fit)[rho_idx]
  alpha_tau <- if (model != "nc") coef(qr_fit)[1] else NA
  
  # Calculate half-life
  if (abs(rho_tau) < 1 && rho_tau > 0) {
    half_life <- log(0.5) / log(rho_tau)
    if (half_life < 0) half_life <- Inf
  } else {
    half_life <- Inf
  }
  
  # Calculate delta^2 (nuisance parameter)
  resid_qr <- residuals(qr_fit)
  ind <- as.numeric(resid_qr < 0)
  phi <- tau - ind
  w <- y  # Using dy as w
  
  cov_phi_w <- cov(phi, w)
  sd_w <- sd(w)
  delta2 <- (cov_phi_w / (sd_w * sqrt(tau * (1 - tau))))^2
  delta2 <- min(max(delta2, 0.1), 1.0)  # Bound between 0.1 and 1.0
  
  # Calculate bandwidth
  h <- qadf_bandwidth(tau, n, hs = TRUE)
  
  # Adjust bandwidth if needed
  if (tau <= 0.5 && h > tau) {
    h <- qadf_bandwidth(tau, n, hs = FALSE)
    if (h > tau) h <- tau / 1.5
  }
  if (tau > 0.5 && h > 1 - tau) {
    h <- qadf_bandwidth(tau, n, hs = FALSE)
    if (h > 1 - tau) h <- (1 - tau) / 1.5
  }
  
  # Fit quantile regressions at tau +/- h for density estimation
  tau_plus <- min(tau + h, 0.99)
  tau_minus <- max(tau - h, 0.01)
  
  qr_plus <- quantreg::rq(y ~ X - 1, tau = tau_plus)
  qr_minus <- quantreg::rq(y ~ X - 1, tau = tau_minus)
  
  # Calculate sparsity (inverse of density)
  mz <- colMeans(X)
  q1 <- sum(mz * coef(qr_plus))
  q2 <- sum(mz * coef(qr_minus))
  fz <- 2 * h / max(q1 - q2, 0.001)
  if (fz < 0) fz <- 0.01
  
  # Build projection matrix for t-statistic
  # X without y_lag1 for projection
  if (p > 0) {
    if (model == "nc") {
      X_proj <- X[, -1, drop = FALSE]
    } else if (model == "c") {
      X_proj <- X[, c(-rho_idx), drop = FALSE]
    } else {
      X_proj <- X[, c(-rho_idx), drop = FALSE]
    }
  } else {
    if (model == "nc") {
      X_proj <- matrix(1, nrow = length(y), ncol = 1)
    } else {
      X_proj <- X[, -rho_idx, drop = FALSE]
    }
  }
  
  # Projection matrix P_X = I - X(X'X)^{-1}X'
  n_obs <- length(y)
  
  # Handle case where X_proj might be rank deficient
  tryCatch({
    if (ncol(X_proj) > 0 && nrow(X_proj) > ncol(X_proj)) {
      XtX_inv <- solve(t(X_proj) %*% X_proj)
      P_X <- diag(n_obs) - X_proj %*% XtX_inv %*% t(X_proj)
    } else {
      P_X <- diag(n_obs)
    }
  }, error = function(e) {
    P_X <- diag(n_obs)
  })
  
  # QADF t-statistic
  fz_crit <- fz / sqrt(tau * (1 - tau))
  y1_vec <- matrix(y_lag1, ncol = 1)
  eq_PX <- sqrt(as.numeric(t(y1_vec) %*% P_X %*% y1_vec))
  
  t_stat <- fz_crit * eq_PX * (rho_tau - 1)
  
  # Get critical values
  cv <- qadf_critical_values(delta2, model)
  
  # Return results
  list(
    tau = tau,
    lags = p,
    alpha = round(alpha_tau, 4),
    rho_tau = round(rho_tau, 4),
    rho_ols = round(rho_ols, 4),
    delta2 = round(delta2, 4),
    half_life = round(half_life, 2),
    t_stat = round(t_stat, 4),
    cv_10 = round(cv["10%"], 4),
    cv_5 = round(cv["5%"], 4),
    cv_1 = round(cv["1%"], 4)
  )
}


#' Calculate optimal bandwidth for QADF
#' @keywords internal
qadf_bandwidth <- function(tau, n, hs = TRUE, alpha = 0.05) {
  x0 <- qnorm(tau)
  f0 <- dnorm(x0)
  
  if (hs) {
    # Hall-Sheather bandwidth
    a <- n^(-1/3)
    b <- qnorm(1 - alpha/2)^(2/3)
    c <- ((1.5 * f0^2) / (2 * x0^2 + 1))^(1/3)
    h <- a * b * c
  } else {
    # Bofinger bandwidth
    h <- n^(-0.2) * ((4.5 * f0^4) / (2 * x0^2 + 1)^2)^0.2
  }
  
  return(h)
}


#' Critical values for QADF test (Hansen 1995)
#' @keywords internal
qadf_critical_values <- function(delta2, model) {
  
  # Critical value tables from Hansen (1995)
  # Rows: delta2 from 0.1 to 1.0
  # Columns: 1%, 5%, 10%
  
  cv_nc <- matrix(c(
    -2.4611512, -1.783209, -1.4189957,
    -2.494341, -1.8184897, -1.4589747,
    -2.5152783, -1.8516957, -1.5071775,
    -2.5509773, -1.895772, -1.5323511,
    -2.5520784, -1.8949965, -1.541883,
    -2.5490848, -1.8981677, -1.5625462,
    -2.5547456, -1.934318, -1.5889045,
    -2.5761273, -1.9387996, -1.602021,
    -2.5511921, -1.9328373, -1.612821,
    -2.5658, -1.9393, -1.6156
  ), nrow = 10, ncol = 3, byrow = TRUE)
  
  cv_c <- matrix(c(
    -2.7844267, -2.115829, -1.7525193,
    -2.9138762, -2.2790427, -1.9172046,
    -3.0628184, -2.3994711, -2.057307,
    -3.1376157, -2.5070473, -2.168052,
    -3.191466, -2.5841611, -2.2520173,
    -3.2437157, -2.639956, -2.316327,
    -3.2951006, -2.7180169, -2.408564,
    -3.3627161, -2.7536756, -2.4577709,
    -3.3896556, -2.8074982, -2.5037759,
    -3.4336, -2.8621, -2.5671
  ), nrow = 10, ncol = 3, byrow = TRUE)
  
  cv_ct <- matrix(c(
    -2.9657928, -2.3081543, -1.9519926,
    -3.1929596, -2.5482619, -2.1991651,
    -3.3727717, -2.7283918, -2.3806008,
    -3.4904849, -2.8669056, -2.5315918,
    -3.6003166, -2.9853079, -2.6672416,
    -3.6819803, -3.095476, -2.7815263,
    -3.7551759, -3.178355, -2.8728146,
    -3.8348596, -3.2674954, -2.973555,
    -3.8800989, -3.3316415, -3.0364171,
    -3.9638, -3.4126, -3.1279
  ), nrow = 10, ncol = 3, byrow = TRUE)
  
  colnames(cv_nc) <- colnames(cv_c) <- colnames(cv_ct) <- c("1%", "5%", "10%")
  
  # Select appropriate table
  cv_table <- switch(model,
                     "nc" = cv_nc,
                     "c" = cv_c,
                     "ct" = cv_ct)
  
  # Interpolate based on delta2
  delta2 <- max(0.1, min(1.0, delta2))
  
  if (delta2 <= 0.1) {
    return(cv_table[1, ])
  } else if (delta2 >= 1.0) {
    return(cv_table[10, ])
  } else {
    # Linear interpolation
    idx <- delta2 * 10
    lower <- floor(idx)
    upper <- ceiling(idx)
    weight <- idx - lower
    
    if (lower < 1) lower <- 1
    if (upper > 10) upper <- 10
    if (lower == upper) {
      return(cv_table[lower, ])
    }
    
    return((1 - weight) * cv_table[lower, ] + weight * cv_table[upper, ])
  }
}


#' QKS critical values
#' @keywords internal
qadf_qks_critical_values <- function() {
  c("10%" = 2.76, "5%" = 3.13, "1%" = 3.74)
}


#' Bootstrap QADF test
#'
#' Computes bootstrap p-value for the QKS statistic.
#'
#' @param object A "qadf" object.
#' @param nboot Number of bootstrap replications. Default 199.
#' @param seed Random seed for reproducibility.
#'
#' @return Updated "qadf" object with bootstrap p-value.
#' @export
qadf_bootstrap <- function(object, nboot = 199, seed = NULL) {
  
  if (!inherits(object, "qadf")) {
    stop("Input must be a 'qadf' object")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  y <- object$y
  p <- object$lags
  tau <- object$tau
  model <- object$model
  n <- length(y)
  
  # Fit AR(p) model to first differences
  dy <- diff(y)
  
  if (p > 0) {
    dy_embed <- embed(dy, p + 1)
    ar_fit <- lm(dy_embed[, 1] ~ dy_embed[, -1] - 1)
    ar_coef <- coef(ar_fit)
    ar_resid <- residuals(ar_fit)
  } else {
    ar_coef <- NULL
    ar_resid <- dy - mean(dy)
  }
  
  # Center residuals
  ar_resid <- ar_resid - mean(ar_resid)
  
  # Bootstrap
  qks_boot <- numeric(nboot)
  
  for (b in 1:nboot) {
    # Resample residuals
    resid_star <- sample(ar_resid, length(ar_resid), replace = TRUE)
    
    # Generate bootstrap dy*
    if (p > 0) {
      dy_star <- dy[1:p]
      for (i in 1:length(resid_star)) {
        dy_new <- sum(ar_coef * dy_star[(length(dy_star) - p + 1):length(dy_star)]) + 
          resid_star[i]
        dy_star <- c(dy_star, dy_new)
      }
    } else {
      dy_star <- resid_star
    }
    
    # Cumulate to get y*
    y_star <- cumsum(c(y[1], dy_star))
    
    # Run QADF on bootstrap sample
    tryCatch({
      qadf_star <- qadf(y_star, model = model, pmax = p, tau = tau)
      qks_boot[b] <- qadf_star$qks
    }, error = function(e) {
      qks_boot[b] <- NA
    })
  }
  
  # Remove NAs
  qks_boot <- qks_boot[!is.na(qks_boot)]
  
  # Calculate p-value
  pvalue <- mean(qks_boot >= object$qks)
  
  object$qks_pvalue <- pvalue
  object$qks_boot <- qks_boot
  object$nboot <- length(qks_boot)
  
  return(object)
}


#' @export
print.qadf <- function(x, ...) {
  cat("\n")
  cat("Quantile Augmented Dickey-Fuller Test\n")
  cat("=====================================\n\n")
  
  cat("Model:", switch(x$model, "nc" = "No constant", "c" = "Constant", 
                       "ct" = "Constant + Trend"), "\n")
  cat("Lags:", x$lags, "(", x$ic, ")\n")
  cat("Sample size:", x$n, "\n\n")
  
  cat("QKS Statistic:", round(x$qks, 4), "\n")
  cat("QKS Critical Values: 10%:", x$qks_cv["10%"], 
      " 5%:", x$qks_cv["5%"], 
      " 1%:", x$qks_cv["1%"], "\n")
  
  if (!is.null(x$qks_pvalue)) {
    cat("Bootstrap p-value:", round(x$qks_pvalue, 4), 
        "(", x$nboot, "replications )\n")
  }
  
  cat("\n")
  invisible(x)
}


#' @export
summary.qadf <- function(object, ...) {
  cat("\n")
  cat("Quantile Augmented Dickey-Fuller Test - Full Results\n")
  cat("====================================================\n\n")
  
  cat("Model:", switch(object$model, "nc" = "No constant", "c" = "Constant", 
                       "ct" = "Constant + Trend"), "\n")
  cat("Lags:", object$lags, "\n")
  cat("Sample size:", object$n, "\n\n")
  
  # Results table
  res <- object$results
  
  # Add significance markers
  res$sig <- ""
  res$sig[res$t_stat < res$cv_10] <- "*"
  res$sig[res$t_stat < res$cv_5] <- "**"
  res$sig[res$t_stat < res$cv_1] <- "***"
  
  print(res, row.names = FALSE)
  
  cat("\nSignificance: *** 1%, ** 5%, * 10%\n")
  cat("\n")
  cat("QKS Statistic:", round(object$qks, 4))
  
  if (object$qks > object$qks_cv["1%"]) {
    cat(" ***\n")
  } else if (object$qks > object$qks_cv["5%"]) {
    cat(" **\n")
  } else if (object$qks > object$qks_cv["10%"]) {
    cat(" *\n")
  } else {
    cat("\n")
  }
  
  cat("QKS Critical Values: 10%:", object$qks_cv["10%"], 
      " 5%:", object$qks_cv["5%"], 
      " 1%:", object$qks_cv["1%"], "\n")
  
  if (!is.null(object$qks_pvalue)) {
    cat("Bootstrap p-value:", round(object$qks_pvalue, 4), "\n")
  }
  
  cat("\n")
  invisible(object)
}


#' @export
plot.qadf <- function(x, ...) {
  
  res <- x$results
  
  # Create plot
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # 1. t-statistics across quantiles
  plot(res$tau, res$t_stat, type = "b", pch = 19, col = "blue",
       xlab = expression(tau), ylab = expression(t[n](tau)),
       main = "QADF t-statistics")
  abline(h = res$cv_5[1], lty = 2, col = "red")
  abline(h = 0, lty = 3, col = "gray")
  legend("topright", "5% CV", lty = 2, col = "red", bty = "n", cex = 0.8)
  
  # 2. Rho estimates
  plot(res$tau, res$rho_tau, type = "b", pch = 19, col = "darkgreen",
       xlab = expression(tau), ylab = expression(hat(rho)(tau)),
       main = "Persistence Parameter")
  abline(h = 1, lty = 2, col = "red")
  abline(h = res$rho_ols[1], lty = 3, col = "blue")
  legend("topright", c("Unit root", "OLS"), lty = c(2, 3), 
         col = c("red", "blue"), bty = "n", cex = 0.8)
  
  # 3. Half-lives
  hl <- res$half_life
  hl[is.infinite(hl)] <- max(hl[is.finite(hl)]) * 1.5
  plot(res$tau, hl, type = "b", pch = 19, col = "purple",
       xlab = expression(tau), ylab = "Half-life",
       main = "Half-lives")
  
  # 4. Original series
  plot(x$y, type = "l", col = "black",
       xlab = "Time", ylab = "Value",
       main = "Original Series")
  
  invisible(x)
}
