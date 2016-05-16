n = 30

true_lambda_1 = 2
true_lambda_2 = 4
true_m = 10


x = rpois(true_m, true_lambda_1)
x = c(x, rpois(n - true_m, true_lambda_2))

plot(1 : n, x, type = "o")
#plot

## hyperparams
alpha = 1
beta = 1

#chains
T = 10000
lambda1s = array(NA, T)
lambda2s = array(NA, T)
ms = array(NA, T)
#start positions
lambda1s[1] = 1
lambda2s[1] = 1
ms[1] = 1

for (t in 2 : T){
	m = ms[t - 1]
	lambda1 = rgamma(1, alpha + sum(x[1 : m]), m + beta)
	lambda2 = rgamma(1, alpha + sum(x[(m + 1) : n]), n - m + beta)
	
	#now we need to calculate all the m dist
	ln_p_m = function(m){
		if (m == 0){			
			sum(x[1 : n]) * log(lambda2) - sum(lgamma(x + 1))
		} else if (m == 30){
			(lambda2 - lambda1) * log(m) + sum(x[1 : m]) * log(lambda1) - sum(lgamma(x + 1))
		} else {
			(lambda2 - lambda1) * log(m) + sum(x[1 : m]) * log(lambda1) + sum(x[(m + 1) : n]) * log(lambda2)  - sum(lgamma(x + 1))
		}			
	}
	ln_m_dist = array(NA, n - 1)
	for (m in 1 : (n - 1)){
		ln_m_dist[m] = ln_p_m(m)
	}
	ln_m_dist	
	ln_m_dist = ln_m_dist - max(ln_m_dist)
	m_dist = exp(ln_m_dist) / sum(exp(ln_m_dist))
	
	
	lambda1s[t] = lambda1
	lambda2s[t] = lambda2
	ms[t] = sample(1 : (n - 1), 1, prob = m_dist)
}



###assess convergence
B = 10
par(mfrow = c(3, 1))

plot(1 : T, lambda1s)
abline(h = mean(lambda1s[B : T]), col = "blue")
abline(h = true_lambda_1, col = "red")
abline(v = B, col = "grey")

plot(1 : T, lambda2s)
abline(h = mean(lambda2s[B : T]), col = "blue")
abline(h = true_lambda_2, col = "red")
abline(v = B, col = "grey")

plot(1 : T, ms)
abline(h = mean(ms[B : T]), col = "blue")
abline(h = sqrt(true_m), col = "red")
abline(v = B, col = "grey")
#plot

##assess autocorrelation

par(mfrow = c(3, 1))

acf(lambda1s[B : T], xlim = c(0, 16))
acf(lambda2s[B : T], xlim = c(0, 16))
acf(ms[B : T], xlim = c(0, 16))
#plot

#burn and thin
lambda1s = lambda1s[B : T]
lambda1s = lambda1s[seq(1, T - B, by = 15)]
lambda2s = lambda2s[B : T]
lambda2s = lambda2s[seq(1, T - B, by = 15)]
ms = ms[B : T]
ms = ms[seq(1, T - B, by = 15)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(3, 1))


hist(lambda1s, br = 500)
abline(v = mean(lambda1s), col = "blue", lwd = 3)
abline(v = quantile(lambda1s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(lambda1s, 0.975), col = "grey", lwd = 3)
abline(v = true_lambda_1, col = "red", lwd = 3)

hist(lambda2s, br = 500)
abline(v = mean(lambda2s), col = "blue", lwd = 3)
abline(v = quantile(lambda2s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(lambda2s, 0.975), col = "grey", lwd = 3)
abline(v = true_lambda_2, col = "red", lwd = 3)

hist(ms, br = 500)
abline(v = mean(ms), col = "blue", lwd = 3)
abline(v = quantile(ms, 0.025), col = "grey", lwd = 3)
abline(v = quantile(ms, 0.975), col = "grey", lwd = 3)
abline(v = true_m, col = "red", lwd = 3)
#plot
