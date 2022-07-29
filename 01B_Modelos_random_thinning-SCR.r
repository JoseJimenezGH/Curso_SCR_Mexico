#==============================================================================#
#                                                                              #
#                CAPTURA-RECAPTURA ESPACIAL CON ACLAREO ALEATORIO              #
#                               José Jiménez                                   #
#                            29/07/2022 17:44:07                               #
#                                                                              #
#==============================================================================#
setwd('C:/Users/Administrator/OneDrive - Universidad de Castilla-La Mancha/00 CURSOS/28 Curso Mexico')
library(nimble)
library(coda)
library(lattice)
library(scrbook)
library(MCMCvis)
library(mcmcOutput)
library(nimble)
source("SCR_functions.R")

#Funciones a usar
e2dist <- function (x, y) {  # Función del paquete scrbook para calcular la
                             # distancias entre localizaciones de 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Simulador de datos
SimSCR0<-function (N = 100, K = 5, lam0 = 0.35, sigma = 0.5, discard0 = TRUE,
 tel=2, n.locs=50, rnd = 2013) {
  set.seed(rnd)
  traplocs <- cbind(sort(rep(1:12, 12)), rep(1:12, 12))
  Dmat <- e2dist(traplocs, traplocs)
  ntraps <- nrow(traplocs)
  buffer <- 1.25
  Xl <- min(traplocs[, 1] - buffer)
  Xu <- max(traplocs[, 1] + buffer)
  Yl <- min(traplocs[, 2] - buffer)
  Yu <- max(traplocs[, 2] + buffer)
  sx <- runif(N, Xl, Xu)
  sy <- runif(N, Yl, Yu)
  S <- cbind(sx, sy)
  D <- e2dist(S, traplocs)
  lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
  plot(traplocs, xlim=c(Xl,Xu), ylim=c(Yl,Yu), pch="+")
  points(S, col="blue", pch=16)
  Y <- array(NA, dim = c(N, ntraps, K))
  for (i in 1:nrow(Y)) {
    for (j in 1:ntraps) {
      Y[i, j, 1:K] <- rpois(K, lam[i, j])
    }
  }
  if (discard0) {
    Y2d <- apply(Y, c(1, 2), sum)
    ncaps <- apply(Y2d, 1, sum)
    Y <- Y[ncaps > 0, , ]
  }

  if (tel > 0) {
    itel <- sort(sample(1:tel, tel, replace = F))
    locs <- list()
    for (i in 1:tel) {
      lx <- rnorm(n.locs, S[itel[i], 1], sigma)
      ly <- rnorm(n.locs, S[itel[i], 2], sigma)
      locs[[i]] <- cbind(lx, ly)
    }
  }
  else {
    locs <- NULL
    itel <- NULL
  }

  list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl, Yu), N = N,
       lam0 = lam0, sigma = sigma, K = K, tel = tel, locs=locs, n.locs=n.locs, S=S)
}

## Simulación de los datos
data <- SimSCR0(N=20, lam0=0.5, K=10, discard0=TRUE, sigma=0.5, tel=2, n.locs=50, rnd=1)

# Centros de actividad
S<-data$S
y<- data$Y  # todos los historiales de captura

# Trampas
X<- data$traplocs
X<-matrix(X, ncol=2)
J<-nrow(X) # número de trampa
K<-data$K  # número de ocasiones de muestreo

(nind <- nrow(apply(y, c(1,2), sum))) # Número de individuos detectados

set.seed(1960)
## Simulación de identificación de los individuos 
id.prob <- 0.2  # ID rate
y.id1<-array(0,c(nind, J, K))
for(i in 1:nind){
  for(j in 1:J){
    for(k in 1:K){
      y.id1[i,j,k]<-rbinom(prob=id.prob, 1, y[i,j,k])
    }
  }
}

sum(y.id1) # eventos con identificación (ID)
## [1] 10
y.id<-y.id1[apply(y.id1,1,sum)>0,,]       # Descartamos las que no tienen ID
(nind2<-nrow(apply(y.id, c(1,3), sum)))   # Número de individuos identificados
M <- 80                                   # Aumentado de datos
yaug<-array(0,c(M,J,K))
yaug[1:nind2,,]<-y.id
sum(yaug)                                 # Eventos con identificación

# Frecuencias de capturas sin identificación
nnid<-apply(y,c(2,3),sum)-apply(y.id, c(2,3), sum)
sum(nnid)                                 # Eventos sin identificación
nnidd<-apply(nnid,1,sum)

# Espacio de estados
xlims <- data$xlim
ylims <- data$ylim
A <- diff(xlims)*diff(ylims)   # Tamaño del espacio de estados

# Ploteado de capturas
plot(X, pch="+",cex=1, xlim=xlims, ylim=ylims, main="", type="n")
tot<-apply(y.id, 2,sum)
## Individuos identificados
symbols(X, circles=tot/5, inches=F, bg="#0000FF3F",  fg=NULL, add=T)
## Individuos sin identificación
nID<-as.numeric(apply(nnid,1,sum))
symbols(X, circles=nID/5, inches=F,bg="#EEAD0E66",fg=NULL, add=T)
points(X, pch="+", cex=1)
points(data$S, col="blue", pch=16)
spiderplot.Over(y.id, X, buffer=0,lwd=2)


library(nimble)
#### Modelo
NimModel <- nimbleCode({

  lam0 ~ dunif(0,5)
  sig ~ dunif(0,5)
  sig2 <- 2*sig^2
  psi ~ dbeta(1,1)
  id.prob ~ dunif(0,1)

  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    z[i] ~ dbern(psi)
    d2[i,1:J] <- (s[i,1]-x[1:J,1])^2 + (s[i,2]-x[1:J,2])^2
    lam[i,1:J] <- lam0*exp(-d2[i,1:J]/sig2)*z[i]

    for(j in 1:J) {
      # Historias completas de encuentro
      y.full[i,j] ~ dpois(lam[i,j]*K)
      # Historias de encuentros aleatorios con identificación
      y.obs[i,j] ~ dbin(id.prob, y.full[i,j])
    }
  }

  # nnid (no identificados) se va a usar el el muestreador IDSampler

  N <- sum(z[1:M])
  D <- N/A
})


## MUESTREADOR PERSONALIZADO METROPOLIS-HASTINGS
## muestreador para actualizar conjuntamente y.un[1:M,j] de manera que ponemos
## en cada paso del muestreo la condición de que sumen n[j]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Definimos los componentes a usar
    nnidd<-control$nnidd
    j<-control$j
    M<-control$M
    calcNodes <- model$getDependencies(target)
  },

run = function() {
  lam.curr <- model$lam[1:M,j] # Conteos esperados de individuos por trampa

  #Muestrea y[1:M,j] reasignando n[j] usando el condicional completo
  switch.probs <- lam.curr[1:M]/sum(lam.curr[1:M])

  #propone nuevas identificaciones para nnid[j,k]
  y.latent.curr <- model$y.full[1:M,j]- model$y.obs[1:M,j]
  y.latent.prop <- rmulti(1, nnidd, switch.probs[1:M])
  model$y.full[1:M,j] <<-  model$y.obs[1:M,j] + y.latent.prop

  # modelo inicial logProb
  model_lp_initial <- model$getLogProb(calcNodes)

  # modelo propuesto logProb
  model_lp_proposed <- model$calculate(calcNodes)

  # Relación log-Metropolis-Hastings
  log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd, switch.probs, log=TRUE)) -
                  (model_lp_initial + dmulti(y.latent.prop,  nnidd, switch.probs, log=TRUE))

  # Paso Metrópolis-Hastings
  accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

## CONSTANTS
str(constants <- list(nnid=nnid,      # eventos sin identificación
                      J=J,            # nº trampas
                      M=M,            # aumentado de datos
                      K=K,            # nº ocasiones
                      A=A))           # tamaño del espacio estados

## DATA
yred<-apply(yaug,c(1,2),sum)
str(data    <-   list(y.obs=yred,     # historias con identificación 
                      xlim=xlims,     # límite X del espacio de estados
                      ylim=ylims,     # límite Y del espacio de estados
                      x=X))           # coordenadas de las trampas


## INICIOS
# Los valores iniciales en este código en NIMBLE deben ser cuidadosamente 
# adaptados para asegurar que todos los nodos, particularmente para el latente 
# y.true, comiencen con valores iniciales razonables. Cualquier error aquí 
# podría invalidar los resultados
ys<-apply(yaug,c(1,2),sum)
s.start <- cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1], ylims[2]))
d <- e2dist(s.start[1:M,], X)
lam0s<- runif(1,0,1)
sigs <- runif(1,0,1)

lam <- lam0s * exp( -(d^2)/(2 * sigs^2))

yi <- array(0, c(M, J, K)) # matriz de eventos
for (j in 1:J) {
  for (k in 1:K) {
    if (nnid[j, k] > 0) {
      probs <- lam[ ,j]
      probs <- probs / sum(probs)
      latent.id <- sample(1:M, nnid[j,k], prob = probs, replace = FALSE)
      yi[latent.id , j, k] <- 1
    }
  } # k
}   # j

yis<-apply(yi,c(1,2),sum) + apply(yaug,c(1,2),sum)
zst<-apply(yis, 1, sum); zst[zst>0]<-1
id.prob.s<-sum(yaug)/(sum(yaug)+sum(nnid))

str(inits   <-   list(z=zst,             # inicios de z
                      s=s.start,         # inicios de s
                      lam0=lam0s,        # ratio base de detección
                      sig=sigs,          # parámetro de movimiento
                      id.prob=id.prob.s, # ratio de identificación
                      psi=runif(1,0,1),  # parámetro de aumentado de datos
                      y.full=yis))       # historias completas latentes

## Parámetros
params <- c('psi', 'lam0', 'sig', 'N', 'D', 'id.prob')

# Ejecutamos el modelo
# Compilamos y ejecutamos el modelo en NIMBLE
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel,
                      constants=constants,
                      data=data,
                      inits=inits,
                      check=FALSE)
Cmodel <- compileNimble(Rmodel)
conf <- configureMCMC(Rmodel,monitors=params, thin=10, useConjugacy = TRUE)

#### -- Hemos cargado antes el muestreador personalizado tras el modelo  -- ####
# reemplazar con un nuevo muestreador para *y* (muestra sin reemplazo con suma 
# resringida a n[j,k])
conf$removeSampler("y.full")
for(j in 1:J){
  conf$addSampler(target = paste(paste("y.full[1:",M,", ",j,"]"), sep=""),
                  type = 'IDSampler',
                  control = list(nnidd = nnidd[j], j=j, M=M),
                  silent = TRUE)
}

# Muestreamos simultáneamente los nodos s[,1] y s[,2]
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  conf$addSampler(target = node,
                  type = "AF_slice",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}


Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# start.time2<-Sys.time()  # Tarda aproximademente 4h la ejecución
# outNim <- runMCMC(Cmcmc, niter = 500000, nburnin = 50000, nchains = 3, 
#                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
#                  samplesAsCodaMCMC = TRUE)
# end.time<-Sys.time()
# end.time-start.time2 # tiempo de ejecución

# Cargo el resultado previamente calculado
load("outNim_RandomThinningSCR.RData")

summary(outNim)

xyplot(outNim)
