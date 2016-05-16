#hyperparams
mu_0 = 5
tau_sq = 100000
n_0 = 0.001
sigsq_0 = 100000

T = 10000

#data
n = 5000
true_mu = 3.45
true_sigsq = 2
xs = rnorm(n, true_mu, sqrt(true_sigsq))
xbar = mean(xs)



thetas = array(NA, T)
sigsqs = array(NA, T)
#start
thetas[1] = 0
sigsqs[1] = .01

for (t in 2 : T){
	#sample mu first
	sigsq = sigsqs[t - 1]
	sigsq_p = 1 / (n / sigsq + 1 / tau_sq)
	theta_p = (n * xbar / sigsq + mu_0 / tau_sq) * sigsq_p
	theta = rnorm(1, theta_p, sqrt(sigsq_p))
	
	#now sample sigsq
	sumsqs = sum((xs - theta)^2)
	sigsq = rinvgamma(1, (n + n_0) / 2, (n_0 * sigsq_0 + sumsqs) / 2)
	
	#record
	thetas[t] = theta
	sigsqs[t] = sigsq
}

B = 500
par(mfrow = c(2, 2))
plot(1 : T, thetas)
abline(h = mean(thetas[B : T]), col = "blue")
abline(h = true_mu, col = "red")
abline(v = B, col = "grey")
plot(1 : T, sqrt(sigsqs))
abline(h = mean(sqrt(sigsqs[B : T])), col = "blue")
abline(h = sqrt(true_sigsq), col = "red")
abline(v = B, col = "grey")
#plot

hist(thetas[B : T], br = 500)
abline(v = mean(thetas[B : T]), col = "blue", lwd = 3)
abline(v = true_mu, col = "red", lwd = 3)
#hist(sqrt(sigsqs[B : T]), br = 500, xlim = c(0, max(sqrt(sigsqs[B : T]))))
hist(sqrt(sigsqs[B : T]), br = 500)
abline(v = sqrt(true_sigsq), col = "red", lwd = 3)
abline(v = mean(sqrt(sigsqs[B : T])), col = "blue", lwd = 3)
#plot
