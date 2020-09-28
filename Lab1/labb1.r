#setwd("~/Studier+Studentliv/labb Bay")
#1a
#likelihood with bernoulli distribution
likelihood = function(theta, s, n){
  like=theta^s*(1-theta)^(n-s)
  return(like)
}

#prior is Beta(2,2) distributed
prior= function(theta, alpha, beta){
  prio = theta^(alpha-1)*(1-theta)^(beta-1)
  return(prio)
}

theta = seq(0,1,0.01)
posterior = prior(theta,2,2) * likelihood(theta,5,20) 
plot(x=theta, y=posterior, type="l")

alpha=2+5
beta=2+15
expected=alpha/(alpha+beta) #=0.29166667
variance=alpha*beta/((alpha+beta)^2*(alpha+beta+1)) #=0.006103845

vec=c()
for (n in 1:20000){
  mean=mean(rbeta(n,2+5,2+20,ncp=0))
  vec=c(vec, mean)
}

plot(x=1:20000, y=vec, type='b')
#We can see that 

#From the graph we can see that the calculated expected value and variance seems accurate. 

#1b
#calculate the probability that alpha > 0.3
#have dist for beta 

theta=rbeta(1000,2+5,2+20,ncp=0)
larger=theta>0.3
mean(larger) #=0.245

#1c 
new_log = log(theta/(1-theta))
hist(new_log) #gives us the amount for each theta
plot(density(new_log)) #show us the density function

#-----------------------------------------------------
#install.packages("invgamma")
library(invgamma)

#2a
y=c(44,25,45,52,30,63,19,50,34,67)
length(y)
mean=3.7
y_mean =mean(y)
nDraws=10000

draws= rchisq(nDraws, df=length(y), ncp = 0)
plot(draws)
sample_var=sum((log(y)-mean)^2)/length(y)
sigma_sq=(length(y))*sample_var/draws

hist(sigma_sq, xlim=c(0,1), breaks=32)

#HOW TO COMPARE?
correct = ((sigma_sq))^(-length(y)/2-1)*exp(-1/(2*(sigma_sq)))
hist(correct, xlim=c(0,10000), breaks=20)

#2b Gini coefficient 
normal=sqrt(sigma_sq) #/sum(sqrt(sigma_sq)) # 
gini=2*pnorm(normal/sqrt(2),0,1)-1 
hist(gini)
mean_gini=mean(gini)
plot(density(gini))

#Very equal income? mean = 5.641896e-05

#2c
sort_gini=sort(gini)
#plot(density(sort_gini[501:9500]))
lower=sort_gini[500] #=3.839527e-05
upper=sort_gini[9500] #=8.335225e-05

#View whole distribution with marks of lower and upper tail
hist(sort_gini, breaks=40)
abline(v=lower, lwd=1, lty=2, col = "red")
abline(v=upper, lwd=1, lty=2, col = "red")

#View 90% of distribution (without tails)
hist(sort_gini[500:9500], breaks=30)

#Do a kernel density estimation to compute a 90% HDP for gini 
sort_density_gini=density(sort_gini)
norm_gini=cumsum(sort_density_gini$y)/sum(sort_density_gini$y) # Normalizing the gini distribution
low_index = min(which(norm_gini>0.05))
low_val=sort_density_gini$x[low_index] #=3.81087e-05
large_index=min(which(norm_gini>0.95))
high_val=sort_density_gini$x[large_index] #=8.3485e-05

hist(gini, breaks=40)
abline(v=low_val, lwd=1, lty=2, col = "red")
abline(v=high_val, lwd=1, lty=2, col = "red")
abline(v=lower, lwd=1, lty=2, col = "blue")
abline(v=upper, lwd=1, lty=2, col = "blue")

#3a
data_radian = c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
mu=2.39
#dist
von_mises = function(kappa, y, mu){
  tall=exp(kappa*cos(y-mu))
  num=2*pi*besselI(x=kappa, nu=0)
  return(tall/num)
}
kappa_vals = seq(0,7, 0.01) #Kappa deside how wuch wind outline (?) - try out values that grasp all essentails in graph below :)
prior = dexp(kappa_vals, 1)

likelihood = function(kappa) {
  return (prod(von_mises(kappa, data_radian, mu)))
}

like_values = sapply(kappa_vals, likelihood)
posterior = like_values * prior
# Normalizing
posterior = posterior/(sum(posterior))#*(kappa_vals[2]-kappa_vals[1]))
plot(x=kappa_vals, y= posterior, type="l", xlab="values for kappa", ylab= "posterior")

#3b 
mode= (which.max(posterior))
#Kappa MODE = 2.12
kappa_vals[mode] #=2.12

