#sim hyperparams
sigma_mh = 1 #MH sampling variance

#data hyperparams
m = 1
mu0 = 0
n0 = 1
beta = 0.5


#data
n = 50
true_theta = 15
true_sigsq = 20
xs = rnorm(n, true_theta, sqrt(true_sigsq))
xbar = mean(xs)


#chains
T = 10000
thetas = array(NA, T)
sigsqs = array(NA, T)
#start vals
thetas[1] = 0
sigsqs[1] = 1

accsigsqs = array(TRUE, T)

ln_k_sigsq_given_theta_and_xs = function(sigsqstar, theta){	
	-(n + n0 - 1) / 2 * log(sigsqstar) - sum((xs - theta)^2) / (2 * sigsqstar) - beta * sigsqstar
}

for (t in 2 : T){
	sigsq = sigsqs[t - 1]
	
	#sample theta first
	theta = rnorm(1, (m * mu0 + n * xbar) / (n + m), sqrt(sigsq / (n + m)))
	
	#sample sigsq using Metropolis
	sigsqstar = rnorm(1, sigsq, sigma_mh)
	ln_r = ln_k_sigsq_given_theta_and_xs(sigsqstar, theta) - ln_k_sigsq_given_theta_and_xs(sigsq, theta)
	if (is.nan(ln_r) || runif(1) > exp(ln_r)){ #reject
		sigsqstar = sigsq
		accsigsqs[t] = FALSE
	} #o/t accept
	
	#record
	thetas[t] = theta
	sigsqs[t] = sigsqstar
}

#assess the efficiency of the MH-steps
sum(accsigsqs) / T


###assess convergence

par(mfrow = c(2, 1))
par(mar = c(2,4,2,0.2))
plot(1 : T, thetas)
#abline(h = mean(thetas[B : T]), col = "blue")
#abline(h = true_theta, col = "red")
#abline(v = B, col = "grey")

plot(1 : T, sigsqs)
#abline(h = mean(sigsqs[B : T]), col = "blue")
#abline(h = true_beta1, col = "red")
#abline(v = B, col = "grey")
#plot

##assess autocorrelation

par(mfrow = c(2, 1))
par(mar = c(1.9,4,1.5,0.2))

acf(thetas[B : T], xlim = c(0, 55), lag.max = 150)
acf(sigsqs[B : T], xlim = c(0, 55), lag.max = 150)
#plot

#burn and thin
B = 400
T = 40
thetas = thetas[B : T]
thetas = thetas[seq(1, T - B, by = T)]
sigsqs = sigsqs[B : T]
sigsqs = sigsqs[seq(1, T - B, by = T)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(2, 1))
par(mar = c(1.9,4,1.5,0.2))


hist(thetas, br = 100)
#abline(v = mean(thetas), col = "blue", lwd = 3)
#abline(v = quantile(beta0s, 0.025), col = "grey", lwd = 3)
#abline(v = quantile(beta0s, 0.975), col = "grey", lwd = 3)
#abline(v = true_beta0, col = "red", lwd = 3)

hist(sigsqs, br = 100)
#abline(v = mean(beta1s), col = "blue", lwd = 3)
#abline(v = quantile(beta1s, 0.025), col = "grey", lwd = 3)
#abline(v = quantile(beta1s, 0.975), col = "grey", lwd = 3)
#abline(v = true_beta1, col = "red", lwd = 3)
#plot


sum(accbeta1s) / T
#print
