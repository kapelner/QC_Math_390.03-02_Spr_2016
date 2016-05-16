library(MCMCpack)

#hyperparams
m = 1
n_0 = 1
sigsq_0 = 1


#data
n = 3
true_beta0 = 2
true_beta1 = -2
true_sigsq = 1
xs = seq(0, 1, length.out = n)
ys = true_beta0 + true_beta1 * xs + rnorm(n, 0, sqrt(true_sigsq))
plot(xs, ys)
mod = lm(ys ~ xs)
coef(mod)
abline(a = coef(mod)[1], b = coef(mod)[2])

#statistics of data (for speed)
sumx = sum(xs)
sumy = sum(ys)
sumxy = sum(xs * ys)
sumxsq = sum(xs^2)
n_plus_m = n + m
sumxsq_plus_m = sumxsq + m

#chains
T = 10000
beta0s = array(NA, T)
beta1s = array(NA, T)
sigsqs = array(NA, T)
#start vals
beta0s[1] = 0
beta1s[1] = 0
sigsqs[1] = 1

for (t in 2 : T){
	sigsq = sigsqs[t - 1]
	beta1 = beta1s[t - 1]
	
	#sample beta_0 first
	beta0 = rnorm(1, (sumy - sumx * beta1) / n_plus_m, sqrt(sigsq / n_plus_m))
	
	#now sample beta_1
	beta1 = rnorm(1, (sumxy - sumx * beta0) / sumxsq_plus_m, sqrt(sigsq / sumxsq_plus_m))
	
	#now sample sigsq
	sigsq = rinvgamma(1, (n + n_0) / 2 + 1, (sum((ys - beta0 - beta1 * xs)^2) + m * (beta0^2 + beta1^2) + n_0 * sigsq_0) / 2)
	
	#record
	beta0s[t] = beta0
	beta1s[t] = beta1
	sigsqs[t] = sigsq
}


###assess convergence
B = 10
par(mfrow = c(3, 1))

plot(1 : T, beta0s)
abline(h = mean(beta0s[B : T]), col = "blue")
abline(h = true_beta0, col = "red")
abline(v = B, col = "grey")

plot(1 : T, beta1s)
abline(h = mean(beta1s[B : T]), col = "blue")
abline(h = true_beta1, col = "red")
abline(v = B, col = "grey")

plot(1 : T, sqrt(sigsqs))
abline(h = mean(sqrt(sigsqs[B : T])), col = "blue")
abline(h = sqrt(true_sigsq), col = "red")
abline(v = B, col = "grey")

##assess autocorrelation

par(mfrow = c(3, 1))

acf(beta0s[B : T], xlim = c(0, 13))
acf(beta1s[B : T], xlim = c(0, 13))
acf(sqrt(sigsqs[B : T]), xlim = c(0, 13))


#burn and thin
beta0s = beta0s[B : T]
beta0s = beta0s[seq(1, T - B, by = 6)]
beta1s = beta1s[B : T]
beta1s = beta1s[seq(1, T - B, by = 6)]
sigsqs = sigsqs[B : T]
sigsqs = sigsqs[seq(1, T - B, by = 6)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(3, 1))


hist(beta0s, br = 500)
abline(v = mean(beta0s), col = "blue", lwd = 3)
abline(v = quantile(beta0s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(beta0s, 0.975), col = "grey", lwd = 3)
abline(v = true_beta0, col = "red", lwd = 3)

hist(beta1s, br = 500)
abline(v = mean(beta1s), col = "blue", lwd = 3)
abline(v = quantile(beta1s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(beta1s, 0.975), col = "grey", lwd = 3)
abline(v = true_beta1, col = "red", lwd = 3)

hist(sqrt(sigsqs), br = 500)
abline(v = sqrt(true_sigsq), col = "red", lwd = 3)
abline(v = mean(sqrt(sigsqs)), col = "blue", lwd = 3)
abline(v = quantile(sqrt(sigsqs), 0.025), col = "grey", lwd = 3)
abline(v = quantile(sqrt(sigsqs), 0.975), col = "grey", lwd = 3)

