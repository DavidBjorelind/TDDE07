## 1

# a)
# Bernoulli trails, n=20, s=5
n=20
s=5
likelihood = function(theta, s, n){
  like = theta^s*(1-theta)^(n-s)
  return(like)
}

# Assuming a Beta(2,2) prior
prior = function(theta, alpha, beta){
  pri = theta^(alpha-1)*(1-theta)^(beta-1)
  return(pri)
}

# Defining random values between 0 and 1
theta = seq(0, 1, 0.01)

vec = c()
for (i in 1:500) {
  theta = seq(0, 1, 1/i)
  mean = mean(prior(theta, 2, 2) * likelihood(theta, s, n))
  vec = c(vec, mean)
}

plot(vec)

posterior = prior(theta, 2, 2) * likelihood(theta, s, n)

plot(x=theta, y=posterior, type = 'l')

# Expected value
ev = (2+5)/(2+5+2+15)
var = (2+5)*(2+15)/(2+5+2+15)^2/(2+5+2+15+1)

# b)

nDraws = 10000
draws = rbeta(nDraws, 2+5, 2+15)
pr = mean(draws>0.3)

# c)

logodds = log(draws/(1-draws))
hist(logodds) # Amount for each theta
plot(density(logodds)) # Shows us the the density function

#---------------

#2
#a)
library(invgamma)

y = c(44,25,45,52,30,63,19,50,34,67)
mean = 3.7

# Simulating from the posterior
nDraws = 10000
draws = rchisq(nDraws, df=length(y))
sample_var = sum((log(y)-mean)^2)/length(y)
sigmasq = length(y)*sample_var/draws
#theta_draw = rnorm(nDraws, mean, sigmasq)

hist(sigmasq/sum(sigmasq))

#correct = 
  
# b)
test = density(draws)
     
