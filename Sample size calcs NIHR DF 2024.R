library(pmsampsize)
library(pmvalsampsize)

# Binary Outcome model Sample size calculation
pmsampsize(
  type = "b",
  nagrsquared = NA,
  csrsquared = NA,
  rsquared = NA,
  parameters = 30,
  shrinkage = 0.9,
  prevalence = 0.20,
  cstatistic = 0.7,
  seed = 1234,
  rate = NA,
  timepoint = NA,
  meanfup = NA,
  intercept = NA,
  sd = NA,
  mmoe = 1.1
)

# Sample size calculation for External Validation
pmvalsampsize(
  type = "b",
  prevalence = 0.20,
  cstatistic = 0.8,
  lpbeta = c(1.8,1.1),
  oeciwidth = 0.2, 
  seed = 1234, 
  simobs = 10000,
  graph = TRUE,
  cstatciwidth = 0.1
)