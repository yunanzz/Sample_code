##########Some Examples#################

########################################
####FPCA: example from FPCA package ####
########################################   

set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
plot(res) # The design plot covers [0, 1] * [0, 1] well.
CreateCovPlot(res, 'Fitted')


################################################################################
###################Example from package "refund" ###############################
###Original site of the code: https://rdrr.io/cran/refund/man/fpca.lfda.html####
################################################################################   

########################################
### Illustration with simulated data ###
########################################   


set.seed(1)
n <- 100 # number of subjects
ss <- seq(0,1,length.out=100) 
TT <- seq(0, 1, length.out=25)
mi <- runif(n, min=6, max=15)
ij <- sapply(mi, function(a) sort(sample(1:25, size=a, replace=FALSE)))

# error variances
sigma <- 0.5

lambdaTrue <- c(1,0.5)  # True eigenvalues
eta1True <- c(0.5, 0.5^2, 0.5^3) # True eigenvalues
eta2True <- c(0.5^2, 0.5^3) # True eigenvalues

phi <- sqrt(2)*cbind(sin(2*pi*ss),cos(2*pi*ss))
psi1 <- cbind(rep(1,length(TT)), sqrt(3)*(2*TT-1), sqrt(5)*(6*TT^2-6*TT+1))
psi2 <- sqrt(2)*cbind(sin(2*pi*TT),cos(2*pi*TT))

zeta1 <- sapply(eta1True, function(a) rnorm(n = n, mean = 0, sd = a))
zeta2 <- sapply(eta2True, function(a) rnorm(n = n, mean = 0, sd = a))

xi1 <- unlist(lapply(1:n, function(a) (zeta1 %*% t(psi1))[a,ij[[a]]] ))
xi2 <- unlist(lapply(1:n, function(a) (zeta2 %*% t(psi2))[a,ij[[a]]] ))
xi <- cbind(xi1, xi2)

Tij <- unlist(lapply(1:n, function(i) TT[ij[[i]]] ))
i <- unlist(lapply(1:n, function(i) rep(i, length(ij[[i]]))))
j <- unlist(lapply(1:n, function(i) 1:length(ij[[i]])))

X <- xi %*% t(phi)
meanFn <- function(s,t){ 0.5*t + 1.5*s + 1.3*s*t}
mu <- matrix(meanFn(s = rep(ss, each=length(Tij)), t=rep(Tij, length(ss)) ) , nrow=nrow(X))

Y <- mu +  X + 
  matrix(rnorm(nrow(X)*ncol(phi), 0, sigma), nrow=nrow(X)) %*% t(phi)  #correlated error


matplot(ss, t(Y[which(i==19),]), type='l', ylab="", xlab="functional argument", 
        main="observations from subject i = 2")

##########################################
# Illustration I : when covariance of scores from a mFPCA step is estimated using fpca.sc
###########################################################################################
est <- fpca.lfda(Y = Y, 
                 subject.index = i, visit.index = j, obsT = Tij, 
                 funcArg = ss, numTEvalPoints = length(TT), 
                 newdata = data.frame(i = c(1:3), Ltime = c(Tij[1], 0.2, 0.5)), 
                 fbps.knots = 35, fbps.p = 3, fbps.m = 2,
                 LongiModel.method='fpca.sc',
                 mFPCA.pve = 0.95, mFPCA.knots = 35, mFPCA.p = 3, mFPCA.m = 2, 
                 sFPCA.pve = 0.95, sFPCA.nbasis = 10, sFPCA.npc = NULL,
                 gam.method = 'REML', gam.kT = 10)

