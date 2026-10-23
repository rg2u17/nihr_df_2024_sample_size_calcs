library(pmsampsize)
library(pmvalsampsize)
library(survival)

# Primary outcome (time to first recurrence)
pmsampsize(type = "s", csrsquared = 0.127, parameters = 51, rate = 0.084,
           timepoint = 5, meanfup = 7.97, shrinkage = 0.9, seed = 1234)   # 3,353
pmsampsize(type = "s", csrsquared = 0.085, parameters = 51, rate = 0.084,
           timepoint = 5, meanfup = 7.97, shrinkage = 0.9, seed = 1234)   # 5,141


# Secondary outcomes (binary, 5-year recurrence)
pmsampsize(type = "b", cstatistic = 0.70, parameters = 51, prevalence = 0.47,
           shrinkage = 0.9, seed = 1234)                                  # 3,571
pmsampsize(type = "b", cstatistic = 0.60, parameters = 51, prevalence = 0.47,
           shrinkage = 0.9, seed = 1234)                                  # 14,744


# External validation - binary
pmvalsampsize(type = "b", prevalence = 0.31, cstatistic = 0.70, lpnormal = c(-0.91, 0.79),
              oeciwidth = 0.2, csciwidth = 0.2, cstatciwidth = 0.1, seed = 1234)  # 4,096
pmvalsampsize(type = "b", prevalence = 0.31, cstatistic = 0.60, lpnormal = c(-0.82, 0.36),
              oeciwidth = 0.2, csciwidth = 0.2, cstatciwidth = 0.1, seed = 1234)  # 14,826

# External validation - TTE - no R package that does this therefore derived from Riley et al. 2022

set.seed(1234)

HORIZON      <- 5
P_FU_GE_H    <- 0.541      # share with >= 5 years of potential follow-up
P_EVENT_BY_H <- 0.185      # share with an observed event within 5 years
P_EVENT_ANY  <- 0.244      # share with an observed event at any time (check only)
C_TARGETS    <- c(0.58, 0.65, 0.70)   # anticipated C-index (UKB primary model ~0.58)
W_SLOPE <- 0.2; W_OE <- 0.2; W_C <- 0.1
N_TUNE  <- 2e5             # simulated cohort used to tune sigma / lambda0
N_BIG   <- 1e6             # simulated cohort used to estimate the SEs

## censoring: Uniform(0, Cmax) with P(C >= 5) = P_FU_GE_H
C_MAX <- HORIZON / (1 - P_FU_GE_H)

simulate <- function(n, sigma, lambda0, z = rnorm(n), u = runif(n), cc = runif(n, 0, C_MAX)) {
  lp <- sigma * z
  tt <- -log(u) / (lambda0 * exp(lp))
  data.frame(lp = lp, time = pmin(tt, cc), status = as.integer(tt <= cc), t_true = tt, cens = cc)
}

## lambda0 so that P(observed event within 5 y) matches, for a given sigma
tune_lambda <- function(sigma, z, u, cc) {
  f <- function(loglam) { d <- simulate(length(z), sigma, exp(loglam), z, u, cc)
                          mean(d$status == 1 & d$time <= HORIZON) - P_EVENT_BY_H }
  exp(uniroot(f, c(-8, 2))$root)
}
## Harrell's C truncated at the horizon
c_index <- function(d) as.numeric(survival::concordance(Surv(time, status) ~ lp, data = d,
                                                         reverse = TRUE, ymax = HORIZON)$concordance)

z0 <- rnorm(N_TUNE); u0 <- runif(N_TUNE); c0 <- runif(N_TUNE, 0, C_MAX)
tune <- function(c_target) {
  g <- function(sigma) { lam <- tune_lambda(sigma, z0, u0, c0)
                         c_index(simulate(N_TUNE, sigma, lam, z0, u0, c0)) - c_target }
  sigma <- uniroot(g, c(0.02, 3), tol = 1e-4)$root
  list(sigma = sigma, lambda0 = tune_lambda(sigma, z0, u0, c0))
}

## SE-based required n for one scenario
required_n <- function(sigma, lambda0) {
  d <- simulate(N_BIG, sigma, lambda0)
  # calibration slope: Cox model of the outcome on the LP (truncated at the horizon)
  dh <- transform(d, status = as.integer(status == 1 & time <= HORIZON), time = pmin(time, HORIZON))
  fit <- coxph(Surv(time, status) ~ lp, data = dh)
  se_slope <- sqrt(vcov(fit)[1, 1])
  # O/E at the horizon: observed = 1 - KM(5); expected = mean predicted risk
  km  <- survfit(Surv(time, status) ~ 1, data = d)
  s5  <- summary(km, times = HORIZON)
  obs <- 1 - s5$surv; se_obs <- s5$std.err                        # Greenwood
  expct <- mean(1 - exp(-lambda0 * exp(d$lp) * HORIZON))
  se_lnoe <- se_obs / obs                                          # delta method, E treated as fixed
  # C-index (truncated): SE from the concordance variance
  cc <- survival::concordance(Surv(time, status) ~ lp, data = d, reverse = TRUE, ymax = HORIZON)
  se_c <- sqrt(cc$var)
  # n for each target width: SE(n) = SE_big * sqrt(N_BIG / n)
  n_for <- function(se_big, half_width) ceiling(N_BIG * (1.96 * se_big / half_width)^2)
  oe <- obs / expct
  # ratio-scale CI width at the observed O/E: oe * (exp(1.96 se) - exp(-1.96 se)) = W_OE
  se_target_oe <- asinh(W_OE / (2 * oe)) / 1.96
  data.frame(sigma = sigma, lambda0 = lambda0,
             c_index = as.numeric(cc$concordance), slope = unname(coef(fit)),
             risk_5y_observed = obs, risk_5y_expected = expct,
             event_any_time = mean(d$status), event_by_5y = mean(d$status == 1 & d$time <= HORIZON),
             n_slope = n_for(se_slope, W_SLOPE / 2),
             n_oe    = ceiling(N_BIG * (se_lnoe / se_target_oe)^2),
             n_cindex = n_for(se_c, W_C / 2))
}

res <- do.call(rbind, lapply(C_TARGETS, function(ct) {
  p <- tune(ct); r <- required_n(p$sigma, p$lambda0); r$c_target <- ct; r
}))
res$events_by_5y_needed <- ceiling(pmax(res$n_slope, res$n_oe, res$n_cindex) * res$event_by_5y)
res$n_required <- pmax(res$n_slope, res$n_oe, res$n_cindex)
cat(sprintf("censoring: Uniform(0, %.2f) years\n\n", C_MAX))
print(within(res[, c("c_target", "c_index", "sigma", "risk_5y_observed", "event_by_5y", "event_any_time",
                     "n_slope", "n_oe", "n_cindex", "n_required", "events_by_5y_needed")],
             { c_index <- round(c_index, 3); sigma <- round(sigma, 3); risk_5y_observed <- round(risk_5y_observed, 3)
               event_by_5y <- round(event_by_5y, 3); event_any_time <- round(event_any_time, 3) }), row.names = FALSE)
cat(sprintf("\ncheck: simulated share with an event at any time vs observed %.3f\n", P_EVENT_ANY))
