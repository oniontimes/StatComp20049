## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)

## -----------------------------------------------------------------------------
plot(lm.D9)

## -----------------------------------------------------------------------------
data <- head(mtcars)

## -----------------------------------------------------------------------------
library(knitr)
library(kableExtra)
data %>% 
  kable() %>%
  kable_styling("striped", full_width = F) %>% 
  column_spec(2:4, bold = T) 

## -----------------------------------------------------------------------------
n <- 5000
u <- runif(n)
x <- 2/(sqrt(1-u))

## -----------------------------------------------------------------------------
x[1]

## -----------------------------------------------------------------------------
hist(x, breaks=200,xlim=c(0,100),  prob = TRUE, main = expression(f(x)==8/x^3))
y <- seq(2, 100, 0.1)
lines(y, 8*y^-3,col="red")

## -----------------------------------------------------------------------------
n <- 500
u <- matrix(runif(n*3,min=-1,max=1),ncol=3)
x <- matrix(0,nrow=n)
for( i in 1:n )
  if(abs(u[i,3])>=abs(u[i,1])&& abs(u[i,3])>=abs(u[i,2])){
    x[i]=u[i,2]
  }else{
    x[i]=u[i,3]
  } 

## -----------------------------------------------------------------------------
x[1]

## -----------------------------------------------------------------------------
hist(x,  prob = TRUE, main = expression(f(x)==0.75*(1-x^2)))

## -----------------------------------------------------------------------------
n <- 500
u <- runif(n)
x <- 2*((1-u)^{-1/4}-1)
hist(x,breaks=100,xlim=c(0,8), prob = TRUE, main = expression(f(x)==64(2+x)^{-5}))
y=seq(0,8,0.01)
lines(y,64*(2+y)^{-5},col="red")

## -----------------------------------------------------------------------------
set.seed(1)
m <- 1e4
x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi / 3

## -----------------------------------------------------------------------------
c(theta.hat, 1/2)

## -----------------------------------------------------------------------------
set.seed(2)
m <- 1e4
x1 <- runif(m, min=0, max=1)
theta_hat_1 <- exp(x1)

## -----------------------------------------------------------------------------
set.seed(3)
m <- 5e4
x2 <- runif(m, min=0, max=1);
theta_hat_2 <- (exp(x2) + exp(1-x2)) / 2

## -----------------------------------------------------------------------------
c(1-var(theta_hat_2) / var(theta_hat_1))

## -----------------------------------------------------------------------------
x<-seq(1.005,5,0.01)
gx<-function(x){
  ((x^2)*exp(-0.5*x^2))/((2*pi)^0.5)*(x>1)
}
plot(x,gx(x),type='l')

## -----------------------------------------------------------------------------
f1<-rnorm(100,0,1)
f2<-rchisq(100,2)
gf1<-gx(f1)/dnorm(f1)
gf2<-gx(f2)/dchisq(f2,2)
cat("the estimation of f1 is",mean(gf1)," and se is",sd(gf1))
cat(" the estimation of f2 is",mean(gf2)," and se is",sd(gf2))

## -----------------------------------------------------------------------------
M <- 100  #number of replicates
g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

## -----------------------------------------------------------------------------
k <- 5 #number of strata
N <- 50 #number of times to repeat the estimation
estimates <- matrix(0, N, 2)
thetahat <- numeric(k)
f <- function(x,j) {
  exp(-x)/(exp(-(j-1)/5)-exp(-j/5))
}
gf <- function(x,j) {
  g(x)/f(x,j)
}
for (i in 1:N) {
  for (j in 1:k) {
    xx <- runif(M/k,(j-1)/5,j/5)
    xxx <- - log(exp(-(j-1)/5) - xx * (exp(-(j-1)/5)-exp(-j/5)))
    thetahat[j] <- mean(gf(xxx,j))
  }
  estimates[i,1] <- sum(thetahat)
  estimates[i,2] <- sd(thetahat)
}
cat("thetahat is", mean(estimates[,1])," estimated standard error is", sd(estimates[,2]))

## -----------------------------------------------------------------------------
alpha = 0.05;
n = 20;
m = 100;

UCL = numeric(m)
LCL = numeric(m)

for(i in 1:m)
{
    x = rlnorm(n,0, 1)
    LCL[i] = mean(log(x)) - qt(alpha / 2, df=n-1, lower.tail = FALSE) * sd(log(x)) / sqrt(n)
    UCL[i] = mean(log(x)) + qt(alpha / 2, df=n-1, lower.tail = FALSE) * sd(log(x)) / sqrt(n)
}

cat("t-interval is", "[",mean(LCL),",",mean(UCL),"]")

## -----------------------------------------------------------------------------
alpha = 0.05;
n = 20;
m = 100;

UCL = numeric(m)
LCL = numeric(m)

for(i in 1:m)
{
    x = rchisq(n, 2)
    LCL[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE) * sd(x) / sqrt(n)
    UCL[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE) * sd(x) / sqrt(n)
}

cat("t-interval is", "[",mean(LCL),",",mean(UCL),"]"," and coverage probability is", mean(LCL < 2 & UCL > 2))

## -----------------------------------------------------------------------------
n <- 20
m <- 100
alpha <- .05
UCL = numeric(m)
LCL = numeric(m)
for(i in 1:m)
{
  x <- rnorm(n, mean=0, sd=2)
  UCL[i] <- mean(x)+ sd(x) *qt(alpha / 2, df=n-1, lower.tail = FALSE)/sqrt(n)
  LCL[i] <- mean(x)- sd(x) *qt(alpha / 2, df=n-1, lower.tail = FALSE)/sqrt(n)
}
cat("the result in example 6.4 is", "[",mean(LCL),",",mean(UCL),"]"," and coverage probability is", mean(LCL < 0 & UCL > 0))

## ----beta---------------------------------------------------------------------

set.seed(12345)

sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

# beta(a,a)
pwr_beta = function(a){
 alpha = 0.1
 n = 20
 m = 1e3
 N = length(a)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rbeta(n, a[j], a[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr * (1-pwr) / m) 
 return(list(pwr = pwr,se = se))
}

 a = c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
 pwr = pwr_beta(a)$pwr
 # plot the power
 se = pwr_beta(a)$se
 plot(a, pwr, type = "b", xlab = "a", ylab = "pwr", pch=16)
 abline(h = 0.1, lty = 2)
 lines(a, pwr+se, lty = 4)
 lines(a, pwr-se, lty = 4)

## ----t------------------------------------------------------------------------

# t(v)
pwr_t = function(v){
 
 alpha = 0.1
 n = 20
 m = 1e3
 N = length(v)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rt(n,v[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr*(1-pwr) / m) 
  return(list(pwr = pwr,se = se))
}

v = seq(1,20)
pwr = pwr_t(v)$pwr
se = pwr_t(v)$se
# plot the power
plot(v, pwr, type = "b", xlab = "v", ylab = "pwr", ylim = c(0,1),pch=16)
abline(h = 0.1, lty = 2)
lines(v, pwr+se, lty = 4)
lines(v, pwr-se, lty = 4)


## -----------------------------------------------------------------------------
count5test <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        return(as.integer(max(c(outx, outy)) > 5))
}
set.seed(12345)
alpha.hat <- 0.055
n <- c(10, 20, 50, 100, 500, 1000)
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
m <- 1e3
result <- matrix(0, length(n), 2)
for (i in 1:length(n)){
  ni <- n[i]
  tests <- replicate(m, expr={
    x <- rnorm(ni, mu1, sigma1)
    y <- rnorm(ni, mu2, sigma2)
    Fp <- var.test(x, y)$p.value
    Ftest <- as.integer(Fp <= alpha.hat)
    c(count5test(x, y), Ftest)
    })
  result[i, ] <- rowMeans(tests)
}
data.frame(n=n, C5=result[, 1], Fp=result[, 2])


## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=100
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library("bootstrap")
n <- nrow(scor)
theta_hat <- cor(law$LSAT, law$GPA)
theta_jack <- numeric(n)
for(i in 1:n){
  x <- scor[-i,]
  theta_jack[i] <- cor(law$LSAT[-i],law$GPA[-i])
}
bias_jack <- (n-1)*(mean(theta_jack)-theta_hat)
se_jack <- sqrt((n-1)*mean((theta_jack-mean(theta_jack))^2))
print(round(c(bias_jack=bias_jack,se_jack=se_jack),4))

## -----------------------------------------------------------------------------
set.seed(12345)
library(boot)
aircondit <- as.matrix(aircondit)
boot.mean <- function(x,i) mean(x[i])
boot.obj <- boot(aircondit, statistic=boot.mean, R=2000)
print(boot.ci(boot.obj, type = c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------

library(bootstrap)
set.seed(12345)

n = nrow(scor)
lambda_hat = eigen(cov(scor))$values
theta_hat = lambda_hat[1] / sum(lambda_hat)
theta_j = rep(0,n)

for (i in 1:n) {

x = scor [-i,]
lambda = eigen(cov(x))$values
theta_j[i] = lambda[1]/sum(lambda)

}
#estimated bias of theta_hat
bias_jack = (n-1)*(mean(theta_j)-theta_hat)
#estimated se of theta_hat
se_jack = (n-1)*sqrt(var(theta_j)/n)

print(round(c(bias_jack=bias_jack,se_jack=se_jack),4))


## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- numeric(n*(n-1)/2)
e2 <- numeric(n*(n-1)/2)
e3 <- numeric(n*(n-1)/2)
e4 <- numeric(n*(n-1)/2)
count <- 0
for (i in 1:(n-1))
  for (j in (i+1):n) {
    count <- count+1
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    P1 <- lm(y~x)
    y1_1 <- chemical[i]*P1$coef[2] + P1$coef[1]
    y1_2 <- chemical[j]*P1$coef[2] + P1$coef[1]
    e1[count] <- (magnetic[i]-y1_1)^2+(magnetic[j]-y1_2)^2
    
    P2 <- lm(y~x+I(x^2))
    y2_1 <- P2$coef[1] + P2$coef[2] * chemical[i] + P2$coef[3] * chemical[i]^2
    y2_2 <- P2$coef[1] + P2$coef[2] * chemical[j] + P2$coef[3] * chemical[j]^2
    e2[count] <- (magnetic[i]-y2_1)^2+(magnetic[j]-y2_2)^2
    
    P3 <- lm(log(y)~x)
    y3_1 <- exp(P3$coef[1] + P3$coef[2] * chemical[i])
    y3_2 <- exp(P3$coef[1] + P3$coef[2] * chemical[j])
    e3[count] <- (magnetic[i]-y3_1)^2+(magnetic[j]-y3_2)^2
    
    P4 <- lm(log(y)~log(x))
    y4_1 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[i]))
    y4_2 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[j]))
    e4[count] <- (magnetic[i]-y4_1)^2+(magnetic[j]-y4_2)^2
  }

e = c(mean(e1)/2,mean(e2)/2,mean(e3)/2,mean(e4)/2)
matrix(e, nrow=1,
       dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))
detach(ironslag)

## -----------------------------------------------------------------------------

set.seed(12345)

# Count Five test
count5test = function(x, y) {
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y))
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
# Count Five test permutation
count5test_permutation = function(z) {

n = length(z)
x = z[1:(n/2)]
y = z[-(1:(n/2))]
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y)) 
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > 5))
}
permutation = function(z,R) {
  n = length(z)
  out = numeric(R)
  for (r in 1: R){
      p = sample(1:n ,n ,replace = FALSE)
      out[r] = count5test_permutation(z[p])
  }
  sum(out)/R
}              


n1 = 20
n2 = 50
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1e3

alphahat1 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
count5test(x, y)
}))
alphahat2 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
permutation(z,1000) 
})<0.05)

round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)


## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0,0,0)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
  mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
rbimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1), mu1, sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2), mu2, sd2)
  return(x)
}
for(i in 1:m){
  mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
  mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=10
n2=100
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------

set.seed(3000)

lap_f = function(x) exp(-abs(x))

rw.Metropolis = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}

N = 200
sigma = c(.05, .5, 2, 16)
x0 = 25
rw1 = rw.Metropolis(sigma[1],x0,N)
rw2 = rw.Metropolis(sigma[2],x0,N)
rw3 = rw.Metropolis(sigma[3],x0,N)
rw4 = rw.Metropolis(sigma[4],x0,N)
#number of candidate points rejected
Rej = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
Acc = round((N-Rej)/N,4)
rownames(Acc) = "Accept rates"
colnames(Acc) = paste("sigma",sigma)
knitr::kable(Acc)
#plot
    rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
    }
    



## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

k <- 4    # four chains
x0 <- c(-10,-5,5,10)    # overdispersed initial values
N <- 10000    # length of chains
b <- 200    # burn-in length

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(0.5,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (1000+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(1000+1):N], type="l", xlab="sigma=0.5", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(1,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (500+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x2 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(500+1):N], type="l", xlab="sigma=1", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(4,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x3 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=4", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rw.Metropolis(16,x0[i],N)$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x4 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=16", ylab="R_hat")
abline(h=1.2, lty=2)

c(x2,x3,x4)

## -----------------------------------------------------------------------------
k = c(4:25,100,500,1000)
S = function(a,k){
 ck = sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

f = function(a,k){S(a,k)-S(a,k-1)}
#curve(f(x),xlim = c(0,sqrt(k)))
a <- seq(0, 4, by=0.01)
plot(a, f(a, k[23]), lty=1, col=1, type="l", xlim=c(0, 4), xlab="a", ylab="f(a|k)", main="f(a) with different k")
lines(a, f(a, k[24]), xlim = c(0, 4), lty=2, col=2)
lines(a, f(a, k[25]), xlim = c(0, 4), lty=3, col=3)
legend("topright", legend=c("k=100", "k=500", "k=1000"), col=1:3,lty=1:3)
# So the lower and upper bound in function uniroot should be 1 and 2 respectively

solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

root = matrix(0,2,length(k))

for (i in 1:length(k)){
  root[2,i]=round(solve(k[i]),4)
}

root[1,] = k
rownames(root) = c('k','A(k)')
root


## -----------------------------------------------------------------------------

library(nloptr)
# Mle 
eval_f0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  return(sum(x)-0.999999)
}

opts = list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle = NULL
r = matrix(0,1,2)
r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j = 2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res = nloptr( x0=c(0.2,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=444,n.B=132,nOO=361,nAB=63 )
j = j+1
r = rbind(r,res$solution)
mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
#the result of EM algorithm
r 
#the max likelihood values
plot(-mle,type = 'l')


## -----------------------------------------------------------------------------

attach(mtcars)

formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
f3 = vector("list", length(formulas))
for (i in seq_along(formulas)){
  f3[[i]] = lm(formulas[[i]], data = mtcars)
}
f3
#2 lapply
la3 = lapply(formulas, function(x) lm(formula = x, data = mtcars))
la3


## -----------------------------------------------------------------------------

set.seed(123)
trials = replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])


## -----------------------------------------------------------------------------
datalist <- list(mtcars, faithful)
lapply(datalist, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE) return(simplify2array(out))
  unlist(out)
}
mylapply(datalist, mean, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  
#  #include <cmath>
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  double f(double x) {
#    return exp(-abs(x));
#  }
#  
#  //[[Rcpp::export]]
#  NumericVector rwMetropolis (double sigma, double x0, int N) {
#    NumericVector x(N);
#    x[0] = x0;
#    NumericVector u = runif(N);
#    for (int i = 1; i < N;i++ ) {
#      NumericVector y = rnorm(1, x[i-1], sigma);
#      if (u[i] <= (f(y[0]) / f(x[i-1]))){
#        x[i] = y[0];
#      }
#      else {
#        x[i] = x[i-1];
#      }
#    }
#    return(x);
#  }
#  

## ----eval=T-------------------------------------------------------------------
    library(Rcpp)
    library(StatComp20049)
    library(microbenchmark)
    # R
    lap_f = function(x) exp(-abs(x))

    rw.Metropolis = function(sigma, x0, N){
    x = numeric(N)
    x[1] = x0
    u = runif(N)
    k = 0
    for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
    else {
    x[i] = x[i-1]
    k = k+1
     }
    }
     return(list(x = x, k = k))
    }
    
    x0 = 25
    N = 2000
    sigma = 2
    (time = microbenchmark(rwR=rw.Metropolis(sigma,x0,N),rwC=rwMetropolis(sigma,x0,N)))

## ----eval=F-------------------------------------------------------------------
#  
#  set.seed(12345)
#  rwR = rw.Metropolis(sigma,x0,N)$x[-(1:500)]
#  rwC = rwMetropolis(sigma,x0,N)[-(1:500)]
#  qqplot(rwR,rwC)
#  abline(a=0,b=1,col='black')

