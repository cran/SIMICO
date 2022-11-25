## ----k2, eval = T, echo=TRUE, warning = FALSE, message = FALSE, error = FALSE----
library(SIMICO)

# Set two outcomes
k = 2

# Set number of observations
n = 5000

# Set number of covariates
p = 2

# Set number of SNPs
q = 50

# Set number of quadrature nodes
d = 100

# Variance of subject-specific random effect
tauSq = 1

# Pairwise correlation
rho = 0.1

# Minor Allele Frequencey Cutoff
Causal.MAF.Cutoff = .05

# Total number of causal variants
num = 5

# Define the effect sizes
effectSizes <- c(.03, .15)

# the baseline cumulative hazard function
bhFunInv <- function(x) {x}

set.seed(1)
# Fixed effects
xMat <- cbind(rnorm(n, mean = 0, sd = 2), rbinom(n=n, size=1, prob=0.5))

# Genetic effects
gMat <- sim_gmat(n, q, rho)

# Get indices to specific select causal variants
idx <- Get_CausalSNPs_bynum(gMat, num, Causal.MAF.Cutoff)

# Subset the gMat
gMatCausal <- gMat[,idx]

# True model has nothing
fixedMat <- matrix(data=0, nrow=n, ncol=k)

# Generate the multiple outcomes
exampleDat <- simico_gen_dat(bhFunInv = bhFunInv, obsTimes = 1:3,
                             windowHalf = 0.1, n = n, p = p, k = k, 
                             tauSq = tauSq, gMatCausal = gMatCausal,
                             xMat = xMat, effectSizes = effectSizes)

# Set the initial estimate values
init_beta <-c (rep(c(0, 0, 0, 1, 0), k), 1)

# Run the newton-raphson
nullFit <- simico_fit_null(init_beta = init_beta, epsilon = 10^-5, 
                           xDats = exampleDat$fullDat$xDats, 
                           lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, 
                           k = k, d = d)

# Get the test statistics p-values
out <- simico_out(nullFit = nullFit$beta_fit, xDats = exampleDat$fullDat$xDats, 
                  lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, 
                  Itt = nullFit$jmat, a1 = 1, a2 = 25, 
                  G = gMat, k  = k, d = d)

#  results
# Score statistic
(out$multQ)

# P-values
(out$multP)

