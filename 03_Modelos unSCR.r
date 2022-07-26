#==============================================================================#
#                                                                              #
#              SCR sin marcas reconocibles (unSCR) o Spatial Counts            #
#                         José Jiménez (CSIC-IREC)                             #
#                           15/02/2020 12:51:32                                #
#                                                                              #
#==============================================================================#

# Con este codigo se trata la informacion con el script del libro de Royle et
# al. (2013)
setwd('C:/Users/Jose/OneDrive/Nimble/00 SCR-actualizados')
library(lattice)
library(coda)
source("Funciones_SCR.R")

library(nimble)
nimble:::setNimbleOption('useSafeDeparse', FALSE)
nimbleOptions('useSafeDeparse')

source("sSampler.R")
 
tr<-seq(1.5,8.5, length=10)
X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. trampas


plot(X, xlim=c(0,10), ylim=c(0,10), pch=3, cex=0.75, xlab="X", ylab="Y", asp=TRUE)

set.seed(2224)
xlim <- c(0,10); ylim <- c(0,10)
N <- 20

S <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(S, pch=16, col=2)

sigma <- 0.5
lambda0 <- 0.4
J <- nrow(X)
K <- 15
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((X[j,1]-S[,1])^2 + (X[j,2]-S[,2])^2)
  lambda <- lambda0*exp(-dist^2/(2*sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}

n <- apply(yy, c(2,3), sum)
n1<-apply(n, 1, sum)
sum(n)
M<-150

# Lo anadimos al ploteado:
tot<-apply(n, 1,sum)
symbols(X, circles=tot/10, inches=F,bg="#00000022", fg=NULL, add=T)
points(X, pch=3, cex=0.75); points(S, pch=16, col=2)
points(S, pch=16, col="red", cex=1.25)

# SC sin informacion adicional es muy poco preciso. Vamos a usar un prior
# informativo para sigma, generandolo con una distribucion gamma. Si queremos
# usar un sigma:
mode = 0.5
sd = 0.1

# Obtenemos los parámetros de ratio y forma:
ra = ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra

show(sh)
show(ra)

# Ploteamos:
X11()
x = seq(0,mode+5*sd,len=1001)
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" ,
      main=paste("dgamma, mode=",mode,", sd=",sd,sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )


library(nimble)
## define the model
code <- nimbleCode({

  sigma ~ dgamma(sh,ra)
  sigma2<-2*sigma^2
  lam0 ~ dunif(0,5)
  psi ~ dunif(0,1)

  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    dist2[i,1:J] <- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    lam[i,1:J] <- lam0*exp(-dist2[i,1:J]/sigma2)*z[i]*K
  }

  for(j in 1:J){
    bigLambda[j] <- sum(lam[1:M,j])
    n[j] ~ dpois(bigLambda[j])
  }
  N <- sum(z[1:M])

})


str(constants <- list(M=M, 
                      K=K, 
                      J=J, 
                      sh=sh, 
                      ra=ra))

str(data    <-   list(n=n1,
                      X=X, 
                      xlim=xlim, 
                      ylim=ylim))

s.start<-cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
z<-rep(1,M)
str(inits   <-   list(lam0=runif(1,0,1), 
                      sigma=runif(1,0,1),  
                      psi=runif(1,0,1),  
                      s=s.start, 
                      z=z))

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)

conf<-configureMCMC(Rmodel, monitors=c('N', 'lam0', 'psi'), thin=1)
conf$removeSampler(paste("s[1:",M,", 1:2]", sep="")) 
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',
                  control=list(i=i, 
                               xlim=Rmodel$xlim, 
                               ylim=Rmodel$ylim, 
                               scale=0.25), #scale parameter here is just the starting scale. It will be tuned.
                  silent = TRUE)
}

MCMC <- buildMCMC(conf)

Comp_MCMC <- compileNimble(MCMC, project = Rmodel)

# Ejecutamos el modelo
nb=1000      # Iteraciones a desechar
ni=5000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(Comp_MCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecucion

# Resultados
summary(outNim)
xyplot(outNim)

gelman.diag(outNim[,c('N','lam0','psi')], multivariate = FALSE)

cat("Población que simulamos = ", N, "individuos", "\n")
cat("Fotografías (todos no identificados)", sum(n), "\n")
