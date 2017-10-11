###Simulating correlated variables
mu <- rep(0,4)
Sigma <- matrix(.7, nrow=4,ncol=4)+diag(4)*.3
library(MASS)
rawvars <- mvrnorm(n=10000,mu = mu, Sigma = Sigma)
head(rawvars)
cor(rawvars); cov(rawvars)
#Transform some variables
pvars <- pnorm(rawvars)

plot(head(pvars,100))

### Simulating with Choletsky Matrix
R = matrix(cbind(1,.80,.2,  .80,1,.7,  .2,.7,1),nrow=3)
U = t(chol(R))
nvars = dim(U)[1]
numobs = 100000
set.seed(1)
random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
orig.raw = as.data.frame(t(random.normal))
names(raw) = c("response","predictor1","predictor2")
cor(raw)
plot(head(raw, 100))
plot(head(orig.raw,100))