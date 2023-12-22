#==============================================================================#
#                                                                              #
#                    CAPTURA-RECAPTURA ESPACIAL (SCR)                          #
#                       José Jiménez (CSIC-IREC)                               #
#                         19/07/2022 22:38:14                                  #
#                                                                              #
#==============================================================================#
setwd('C:/Users/Administrator/OneDrive - Universidad de Castilla-La Mancha/00 CURSOS/28 Curso Mexico')
source("SCR_functions.R")
library(scrbook)
library(lattice)
library(coda)
library(mcmcOutput)
library(nimble)

## Simulación de datos
N <- 50  # Tamaño de población
K <- 15  # Ocasiones de muestreo
J <- 100 # número de trampas

data <- SimSCR0(N=N, K=K, array3d = TRUE, discard0=TRUE, rnd=2013)

p0<-data$p0
sigma<-data$sigma
N<-data$N

S<-data$S
y3d<-data$Y
(detections <- sum(y3d))

(nind <- dim(y3d)[1])
X <- data$traplocs
K <- data$K
J <- nrow(X) 
M <- 150

# Espacio de estados
xlim <- data$xlim
ylim <- data$ylim
area <- diff(xlim)*diff(ylim)

y <- apply(y3d,c(1,2),sum); sum(y)
yaug<-array(0,c(M,J))
yaug[1:nind,]<-y

# Ploteamos los capturados vs no capturados (puntos sólidos vs círculos)
plot(X, xlim=xlim, ylim=ylim, pch="+", asp=TRUE)
points(S, pch=1, col="red", cex=1.5)
captured<-apply(y3d,1,sum)
captured[captured>1]<-1
captured<-(1:nind)*captured

points(S[captured,], pch=16, col="red", cex=1.5)
datn<-apply(y3d, c(2,3), sum)
tot<-apply(datn, 1,sum)
symbols(X, circles=tot/5, inches=F, bg="#00000022", fg=NULL, add=T)
points(X, pch="+", cex=0.5)


# Spiderplot
spiderplotJJ4(y3d, X, buffer=2, lwd=2)


## definimos el modelo
code <- nimbleCode({
  
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,10)
  psi ~ dbeta(1,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi) 
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J] <- pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2)
    p[i,1:J] <- p0*exp(-d2[i,1:J]/(2*sigma^2))*z[i]

    for(j in 1:J){
      y[i,j] ~ dbinom(p[i,j],K)
      #y[i,j] ~ dpois(p[i,j]*K) # Alternativamente podemos usar una
                                # distribución de Poisson para trabajar
                                # con conteos
    }
  }
  N <- sum(z[1:M])
  D <- N/area
})

str(constants <- list(M=M,        # aumentado de datos 
                      K=K,        # ocasiones de muestreo
                      J=J,        # número de trampas
                      area=area)) # área del espacio de estados

str( data   <-   list(y=yaug,     # matriz de capturas reducida sobre K
                      X=X,        # matriz de coordenadas de las trampas
                      xlim=xlim,  # extremos x del espacio de estados
                      ylim=ylim)) # extremos y del espacio de estados

# Inicios para las ubicaciones latentes
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X[yaug[i,]>0,1] )
  sst[i,2] <- mean( X[yaug[i,]>0,2] )
}
zd<-c(rep(1, nind), rbinom((M-nind),1,0.2))

str( inits   <-  list(p0=runif(1,0,1),    # probabilidad basal de detección
                      sigma=runif(1,0,2), # parametrización de sigma
                      psi=runif(1,0,1),   # parámetro de aumentado de datos
                      s=sst,              # ubicaciones de inicio
                      z=zd))              # está o no en la población

params <- c('N', 'D', 'sigma','psi', 'p0','s','z')

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()  # Conveniente para ver los inicios
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)
mcmc<-configureMCMC(Rmodel, monitors=params)

# Para muestrear conjuntamente las coordenadas x, y de s
mcmc$removeSamplers("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  mcmc$addSampler(target = node,
                  type = "AF_slice",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

MCMC <- buildMCMC(mcmc)

CompMCMC <- compileNimble(MCMC, project = Rmodel)

nb = 1000
ni = 5000 + nb
nc = 3

start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains =  nc,inits=inits,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

# Resultados
summary(outNim[,c('N','D','p0','psi', 'sigma')])

# Inspeccionamos la convergencia de las cadenas de Markov
#xyplot(outNim[,c('N','p0','psi', 'sigma')])
diagPlot(outNim[,c('N','p0','psi', 'sigma')])

gelman.diag(outNim[,c('N','p0','psi', 'sigma')], multivariate = FALSE)

cat("Población que simulamos = ", N, "individuos", "\n")
cat("lambda0 simulada  = ", p0, "\n")
cat("sigma simulada  = ", sigma, "\n")
cat("Datos recogidos = ", sum(y), "foto-capturas", "\n")


samplesn<-data.matrix(outNim)

# Coeficiente variación
sd(samplesn[,2])/mean(samplesn[,2])

# Histograma
hist(samplesn[,2], col="grey90")

################################################################################

# VISUALIZACIÓN ESPACIAL
#========================
mco <- mcmcOutput(outNim)
s1<-mco$s[,,1] ; Sx<-as.matrix(s1)
s2<-mco$s[,,2] ; Sy<-as.matrix(s2)
z<-mco$z       ; z<- as.matrix(z)

delta<-1.5
Xl<-min(X[,1])-delta
Xu<-max(X[,1])+delta
Yl<-min(X[,2])-delta
Yu<-max(X[,2])+delta
obj<-list(Sx=Sx,Sy=Sy,z=z)

par(mfrow=c(1,1))
dev.new(width=8,height=8)
Spat<-SCRdensity(obj, nx=100, ny=100, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu)
points(X, pch="+", col="white")
points(S, col="red", pch=16)
