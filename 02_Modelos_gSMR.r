#==============================================================================#
#                                                                              #
#             MARCAJE-REAVISTAMIENTO ESPACIAL GENERALIZADO (gSMR)              #
#                       José Jiménez (CSIC-IREC)                               #
#                         19/07/2022 22:38:14                                  #
#                                                                              #
#==============================================================================#

setwd('C:/Users/Administrator/OneDrive - Universidad de Castilla-La Mancha/00 CURSOS/28 Curso Mexico')

library(nimble)
nimble:::setNimbleOption('useSafeDeparse', FALSE)
nimbleOptions('useSafeDeparse')
library(mcmcOutput)
library(scrbook)
source("SCR_functions.R")

## Ubicación de los 9 dispositivos de captura (coordenadas)
xmn <- 3
xm0 <- seq(0.45, 0.55, length=xmn)
xm <- cbind(rep(xm0, each=xmn), rep(xm0, times=xmn))
str(xm)


## Ubicación de los 36 dispositivos de re-avistamiento (coordenadas)
xrn <- 6
xr0 <- seq(0.2, 0.80, length=xrn)
xr <- cbind(rep(xr0, each=xrn), rep(xr0, times=xrn))
str(xr)

## Espacio de estados
xlim <- ylim <- 0:1

## Simulamos los centros de actividad
##set.seed(34097)
set.seed(1960)
N <- 50
s <- cbind(runif(N, xlim[1], xlim[2]),
           runif(N, ylim[1], ylim[2]))


## Datos de captura
set.seed(77824)
KTm<-10
p0m <- 0.1
sigma<- 0.06
Jm <- nrow(xm)
ym.all <- matrix(NA, N, Jm)
for(i in 1:N) {
  dist.sxm <- sqrt((s[i,1]-xm[,1])^2 + (s[i,2]-xm[,2])^2)
  pm <- p0m*exp(-dist.sxm^2 / (2*sigma^2))
  ym.all[i,] <- rbinom(Jm, KTm, pm)
}

marked <- rowSums(ym.all)>0
(nMarked <- sum(marked))
ym.obs <- ym.all[marked,]
sum(ym.obs)


## Datos de reavistamiento (marcados) y conteos (no marcados)
##set.seed(77824)
set.seed(34980)
KTr<-20
lam0r <- 0.4
Jr <- nrow(xr)
yr.all <- matrix(NA, N, Jr)
for(i in 1:N) {
  dist.sxr <- sqrt((s[i,1]-xr[,1])^2 + (s[i,2]-xr[,2])^2)
  lamr <- lam0r*exp(-dist.sxr^2 / (2*sigma^2))
  yr.all[i,] <- rpois(Jr, lamr*KTr)
}
yr.obs <- yr.all[marked,]
str(yr.obs)
# Observaciones totales de marcados
sum(yr.obs)

(n <- colSums(yr.all*(1-marked)))  # Conteos de animales no marcados
str(n)
# Observaciones totales de no marcados
sum(n)

# Ploteamos las capturas
dev.new(width=7, height=7.5)
plot(xm, pch="+", xlim=xlim, ylim=ylim, xlab="X", ylab="Y", asp=1)
points(xr, pch=3)  # Dispositivos sw reavistamiento (cámaras-trampa)
rect(xlim[1], ylim[1], xlim[2], ylim[2], lty=3)
points(s, pch=1, cex=1.5, col="red")  # Ubicación de los centros de actividad (CA)
points(s[marked,], pch=16, cex=1.5, col="red") # CA de marcados
points(xr, cex=n, pch=16, col=rgb(0,0,1,0.2)) # Algunos marcados y otros no


## Aumentado de datos
M <- 200

## Preparación de datos de marcaje
ym.aug <- array(0L, c(M, Jm))
ym.aug[1:nrow(ym.obs),] <- ym.obs

## Preparación de datos de reavistamiento
yr.aug <- array(0L, c(M, Jr))
yr.aug[1:nrow(yr.obs),] <- yr.obs


## define the model
code <- nimbleCode({

  p0m ~ dunif(0,1)
  sig ~ dunif(0,2)
  sig2 <- 2*sig^2
  lam0r ~ dunif(0,5)
  psi ~ dbeta(1,1)

  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    z[i] ~ dbern(psi)
    # Detección en trampas de marcaje
    pm[i,1:Jm] <- GetDetectionRate(s = s[i,1:2], X = xm[1:Jm,1:2], J=Jm, sigma=sig, lam0=p0m, z=z[i])
    # Detección en cámaras-trampa
    lamr[i,1:Jr] <- GetDetectionRate(s = s[i,1:2], X = xr[1:Jr,1:2], J=Jr, sigma=sig, lam0=lam0r, z=z[i])

    ## Datos en las ocasiones de marcaje (para N latente)
    for(j in 1:Jm) {  
      ym[i,j] ~ dbinom(pm[i,j],KTm)
    }
  }

  ## Datos en las ocasiones de reavistamiento de marcados (1:nMarked)
  for(i in 1:nMarked) {
    for(j in 1:Jr) {
      yr[i,j] ~ dpois(lamr[i,j]*KTr) # Aquí z se aplica a 1:nMarked
    }
  }

  ## Conteo de animales no marcados (desde n:Marked hasta N latente)
  for(j in 1:Jr) {
    # Si podemos suponer que no hay variación temporal, sumamos:
    Lam[j] <- sum(lamr[((nMarked+1):M),j]) # lamr desde nMarked:M
    n[j] ~ dpois(Lam[j]*KTr)
  }

  N <- sum(z[1:M])
})

# Función para calcular el ratio de detección, but omitiéndolo cuando z=0
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

source("sSampler.R")

# Preparamos los inicios para s y z
s.start <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
zd <- rep(1,M)

str(constants<-list(Jm=Jm, 
                    Jr=Jr,
                    KTm=KTm,
                    KTr=KTr, 
                    M=M,
                    nMarked=nMarked,
                    A=diff(xlim)*diff(ylim)))

str(data   <-  list(xm=xm, 
                    xr=xr,
                    ym=ym.aug, 
                    yr=yr.aug,
                    xlim=xlim, 
                    ylim=ylim,
                    n=n))
             
str(inits  <-  list(z=zd, 
                    s=s.start,
                    psi=runif(1,0,1), 
                    p0m=runif(1,0,1), 
                    lam0r=runif(1,0,1),
                    sig=runif(1,0,2)))

params <- c("psi", "p0m", "lam0r", "sig", "N")

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)
conf<-configureMCMC(Rmodel, monitors=params)

conf$removeSampler(paste("s[1:",M,", 1:2]", sep="")) 
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',
                  control=list(i=i, 
                               xlim=Rmodel$xlim, 
                               ylim=Rmodel$ylim, 
                               scale=0.025), #scale parameter here is just the starting scale. It will be tuned.
                  silent = TRUE)
}

MCMC <- buildMCMC(conf)

CompMCMC <- compileNimble(MCMC, project = Rmodel)

nb = 1000
ni = 5000 + nb
nc = 3

start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains = nc,inits=inits,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecucion

summary(outNim)
diagPlot(outNim)

