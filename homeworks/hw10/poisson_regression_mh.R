#hyperparams
sigma_mh = 0.5 #MH sampling variance


#data
n = 50
true_beta0 = 1
true_beta1 = 3
ts = seq(0, 2, length.out = n)
ys = array(NA, n)
for (i in 1 : n){
	ys[i] = rpois(1, true_beta0 + true_beta1 * ts[i])
}
plot(ts, ys)
#plot


mod = lm(ys ~ ts)
coef(mod)
abline(a = coef(mod)[1], b = coef(mod)[2])

#chains
T = 100000
beta0s = array(NA, T)
beta1s = array(NA, T)
#start vals
beta0s[1] = 2
beta1s[1] = 2
accbeta0s = array(TRUE, T)
accbeta1s = array(TRUE, T)

ln_p_beta_0_beta_1_given_y_t = function(b0, b1){
	sum(log(dpois(ys, b0 + b1 * ts)))
}

for (t in 2 : T){
	beta0 = beta0s[t - 1]
	beta1 = beta1s[t - 1]
	
	#sample beta_0 first
	beta0star = rnorm(1, beta0, sigma_mh)
	#calc r
#	r = p_beta_0_beta_1_given_y_t(beta0star, beta1) / dnorm(beta0, beta0star, sigma_mh) / 
#			(p_beta_0_beta_1_given_y_t(beta0, beta1) / dnorm(beta0star, beta0, sigma_mh))
	ln_r = ln_p_beta_0_beta_1_given_y_t(beta0star, beta1) - ln_p_beta_0_beta_1_given_y_t(beta0, beta1)
	if (is.nan(ln_r) || runif(1) > exp(ln_r)){
		#reject
		beta0star = beta0
		accbeta0s[t] = FALSE
	} #o/t accept
	
	#sample beta1 next
	beta1star = rnorm(1, beta1, sigma_mh)
	#calc r
#	r = p_beta_0_beta_1_given_y_t(beta0star, beta1star) / dnorm(beta1, beta1star, sigma_mh) / 
#			(p_beta_0_beta_1_given_y_t(beta0star, beta1) / dnorm(beta1star, beta1, sigma_mh))
	ln_r = ln_p_beta_0_beta_1_given_y_t(beta0star, beta1star) - ln_p_beta_0_beta_1_given_y_t(beta0star, beta1)
	if (is.nan(ln_r) || runif(1) > exp(ln_r)){
		#reject
		beta1star = beta1
		accbeta1s[t] = FALSE
	} #o/t accept
	
	#record
	beta0s[t] = beta0star
	beta1s[t] = beta1star
}


###assess convergence
B = 1000
par(mfrow = c(2, 1))

plot(1 : T, beta0s)
abline(h = mean(beta0s[B : T]), col = "blue")
abline(h = true_beta0, col = "red")
abline(v = B, col = "grey")

plot(1 : T, beta1s)
abline(h = mean(beta1s[B : T]), col = "blue")
abline(h = true_beta1, col = "red")
abline(v = B, col = "grey")
#plot

##assess autocorrelation

par(mfrow = c(2, 1))

acf(beta0s[B : T], xlim = c(0, 100))
acf(beta1s[B : T], xlim = c(0, 100))
#plot

#burn and thin
beta0s = beta0s[B : T]
beta0s = beta0s[seq(1, T - B, by = 52)]
beta1s = beta1s[B : T]
beta1s = beta1s[seq(1, T - B, by = 52)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(2, 1))


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
#plot

#assess the efficiency of the MH-steps
sum(accbeta0s) / T
sum(accbeta1s) / T
#print

