set.seed(1995)
n <- 1e4 # Number of correlated z-scores for each simulated sample
m <- 1e3 # Number of simulation samples

## The correlation matrix has a low-rank structure
## Sigma = BB' + I
## B is n * d, d << n
d <- 5
B <- matrix(rnorm(n * d), n, d)
# Sigma <- cov2cor(B %*% t(B) + diag(n))

## Simulate Z ~ N(0, Sigma)
## An equivalent fast way to simulate Z
## Z = BX + E, and then normalize
## X is d indep N(0, 1)
## E is n indep N(0, 1)
## This way we get n correlated null z-scores
normalize.constant <- sqrt(rowSums(B^2) + 1)
Z <- (colSums(t(B) * rnorm(d)) + rnorm(n)) / normalize.constant

## We keep doing it m = 1K times
## Each time we get n = 10K correlated null z-scores
## Now we have m iid samples from N(0, Sigma)
## Z is n * m, each column is a sample from N(0, Sigma)
Z <- replicate(m,
  (colSums(t(B) * rnorm(d)) + rnorm(n)) / normalize.constant
)

## Categorize m samples into three groups
## Inflated Noise: sd > 1.05
## Deflated Noise: sd < 0.95
## In-between: the rest
sd.Z <- apply(Z, 2, sd)
group <- cut(sd.Z, breaks = c(0, 0.95, 1.05, Inf), labels = c("Deflated Noise", "In-between", "Inflated Noise"))
table(group)

## Add signals from 0.9delta0 + 0.1N(0, 4^2)
## Add the same signals to all samples
## X = theta + Z
## Computer two-sided p-values
theta <- c(rep(0, 0.9 * n), rnorm(0.1 * n, 0, 4))
X <- theta + Z
p <- 2 * pnorm(-abs(X))

## Set nominal FDR = 0.1
## Apply BH
q <- 0.1
p.BH <- apply(p, 2, p.adjust, method = "BH")
BH.rej <- apply(p.BH, 2, function(x){x <= q})

## Calculate FDP for each sample
FDP.fun <- function (rej.id, theta) {
  sum(theta[rej.id] == 0) / max(1, length(theta[rej.id]))
}
FDP <- apply(BH.rej, 2, FDP.fun, theta)

## Overall FDP for 3 groups combined
mean(FDP)
t.test(FDP)$conf.int

## Boxplot of FDP for each group
boxplot(FDP ~ group, ylab = "FDP")
abline(h = q, lty = 2, col = "red")
points(seq(nlevels(group)), tapply(FDP, group, mean), col = "blue", pch = 13)
## Mean of FDP for each group
tapply(FDP, group, mean)
## sd of FDP for each group
tapply(FDP, group, sd)
## 95% CI of the mean of FDP for each group
tapply(FDP, group, function(x){t.test(x)$conf.int})
## 95% CI of the mean of FDP overall
t.test(FDP)$conf.int
