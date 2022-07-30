#==============================================================================#
#                                                                              #
#                   MARCAJE-REAVISTAMIENTO ESPACIAL CON                        #
#                     MARCAS PARCIALMENTE IDENTIFICADAS                        #
#                         José Jiménez (CSIC-IREC)                             #
#                           15/02/2020 12:51:32                                #
#                                                                              #
#==============================================================================#

setwd('C:/Users/Administrator/OneDrive - Universidad de Castilla-La Mancha/00 CURSOS/28 Curso Mexico')
source("SCR_functions.R")
library(adehabitatHR)
library(coda)
library(lattice)

M<-750  # Aumentado de datos

# Ocasiones de muestreo
# (06/3/2014-30/3/2014)
Oper.mark<- data.matrix(read.table("Oper_trap.txt", header=FALSE))
dim(Oper.mark)
Oper.mark.R<-apply(Oper.mark,1,sum)

## 1. PROCESO DE MARCADO
##========================
mark.ch <- secr::read.capthist("ymark.txt", "trapssecr.txt", noccasions =21)
str(mark.ch)
y.mark3D<-as.array(mark.ch)
y.mark<-apply(y.mark3D,c(1,3),sum)
dim(y.mark)
(n.marked<-dim(y.mark)[1])

y.mark3D.test<-aperm(y.mark3D, c(1,3,2))
sum(apply(y.mark3D.test,c(2,3),sum)); dim(y.mark3D.test)
sum(apply(y.mark3D.test,c(2,3),sum)* Oper.mark)
which(apply(y.mark3D.test,c(2,3),sum)* Oper.mark-apply(y.mark3D.test,c(2,3),sum)<0, arr.in=TRUE)

(J.mark<-dim(Oper.mark)[1])
(K.mark<- dim(Oper.mark)[2])

y.mark.aug<-array(0,c(M,J.mark))
y.mark.aug[1:11,]<-y.mark

# Ploteado de la operatividad
x1<-as.matrix(Oper.mark)  ## Tenemos que convertir esto en matriz
image(1:K.mark,1:J.mark,t(x1), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Live trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J.mark, by=2)))
box()

## 2. PROCESO DE REAVISTAMIENTO
##===============================
# Operatividad de las cámaras-trampa
Oper.resight<-  data.matrix(read.table("Oper_cam.txt", header=FALSE))
(K <- dim(Oper.resight)[2])
(J.resight <- dim(Oper.resight)[1])

# Ploteado de la operatividad
x2<-as.matrix(Oper.resight)  ## matriz
image(1:K,1:J.resight,t(x2), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Camera trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J.resight, by=2)))
box()

## 2.1. MARCADOS PERO NO IDENTIFICADOS (TODOS)
nnid2<- secr::read.capthist("nnid.txt", "traplocs.txt", noccasions =49)
nnid2<-aperm(nnid2,c(1,3,2))
nnid<-data.matrix(apply(nnid2,c(2,3),sum))


x.resight<-data.matrix(read.table("traplocs.txt", header=FALSE))[,2:3]/1000
Xx<-mean(x.resight[,1])
Xy<-mean(x.resight[,2])
x.resight[,1]<-(x.resight[,1]-Xx)
x.resight[,2]<-(x.resight[,2]-Xy)

livetraps<-as.matrix(secr::traps(mark.ch))
x.mark<-data.matrix(livetraps)/1000

x.mark[,1]<-(x.mark[,1]-Xx)
x.mark[,2]<-(x.mark[,2]-Xy)


# Ploteado de las capturas
plot(x.resight,
     xlim=c(min(x.resight[,1])- 1,max(x.resight[,1])+ 1),
     ylim=c(min(x.resight[,2])- 1,max(x.resight[,2])+ 1),
     pch="+", 
     xlab="X",
     ylab="Y",
     asp=1)
points(x.mark, pch="+", col="red")

dim(Oper.resight)

# Probamos si hay errores en la operatividad:
which((nnid*Oper.resight)-nnid<0, arr.in=TRUE)


## 2.2. SIN MARCAJE
n<-secr::read.capthist("unmarked.txt", "traplocs.txt", noccasions =49)
n<-data.matrix(apply(n,c(3,2),sum))

# Probamos si hay errores en la operatividad:
which(n*Oper.resight -n<0, arr.in=TRUE)

# Ploteamos las detecciones
tot<-apply(nnid, 1,sum)
symbols(x.resight, circles=tot/10, inches=F,bg="#EEC90033", fg=NULL, add=T)
tot2<-apply(n,1,sum)
symbols(x.resight, circles=tot2/10, inches=F,bg="#228B2219", fg=NULL, add=T)


## MÁSCARA DE FECHA DE MARCAJE
mark.date<-data.matrix(read.table("mark_date.txt", header=FALSE))
mark.date.aug<-array(0,c(M,49))
mark.date.aug[1:11,]<-mark.date

# ratio de capturas:
(sum(nnid)+sum(n))*100/sum(Oper.resight)


# TELEMETRIA
#===============
library(adehabitatHR)
locs<-read.table("PosZorro.txt", header=FALSE)[,1:2]/1000
locs[,1]<-locs[,1]-Xx
locs[,2]<-locs[,2]-Xy
locs1<-locs[1:31,]
locs2<-locs[32:461,]
locs3<-locs[462:527,]
locs4<-locs[528:546,]

points(locs1[,1], locs1[,2], pch=16, col="#FF00004C", cex=1)
points(locs2[,1], locs2[,2], pch=16, col="#0000FF4C", cex=1)

xy1<-cbind(locs1[,1], locs1[,2])
xy1<-data.matrix(xy1)
xysp1 <- SpatialPoints(xy1)
kudl1 <- kernelUD(xysp1, grid=100)
homerange1 <- getverticeshr(kudl1); homerange1@data$area*1E6
contour(getvolumeUD(kudl1)[1], level=c(50,95), add=TRUE, col=c('blue','blue'), lwd=2)

xy2<-cbind(locs2[,1], locs2[,2])
xy2<-data.matrix(xy2)
xysp2 <- SpatialPoints(xy2)
kudl2 <- kernelUD(xysp2, grid=100)
homerange2 <- getverticeshr(kudl2); homerange2@data$area*1E6
contour(getvolumeUD(kudl2)[1], level=c(50,95), add=TRUE, col=c('blue','blue'), lwd=2)

points(locs3[,1], locs3[,2], pch=16, col="#00FF004C", cex=1)
points(locs4[,1], locs4[,2], pch=16, col="forestgreen", cex=1)

locs<-rbind(locs1,locs2,locs3,locs4)
ind<-c(rep(1,nrow(locs1)),rep(2,nrow(locs2)),rep(3,nrow(locs3)),rep(4,nrow(locs4)))
nlocs<-nrow(locs)
n.collar<-4


buff<- 2
xl<-min(x.resight[,1])-buff
xu<-max(x.resight[,1])+buff
yl<-min(x.resight[,2])-buff
yu<-max(x.resight[,2])+buff
xlims=c(xl, xu)
ylims=c(yl, yu)
area<-(xu-xl)*(yu-yl)


library(nimble)
#### Modelo en NIMBLE
NimModel <- nimbleCode({

  lam0.resight ~ dunif(0,2) # ratio base de detección en reavistamientos
  lam0.mark  ~ dunif(0,1)   # ratio base de detección en capturas
  sigma ~ dunif(0,2)        # parámetro escala de la seminormal
  sigma2 <- sigma^2
  id.prob ~ dunif(0,1)      # ratio de identificación de marcados
  psi ~ dbeta(1,1)          # parámetro del aumentado de datos

  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    # Información auxiliar para datos de marcaje (trampas)
    d.mark2[i,1:J.mark] <- (s[i,1]-x.mark[1:J.mark,1])^2 + (s[i,2]-x.mark[1:J.mark,2])^2
    lam.mark[i,1:J.mark] <- lam0.mark*exp(-d.mark2[i,1:J.mark]/(2*sigma2))*z[i]
    # Información auxiliar para datos de reavistamiento (cámaras-trampa)
    d.resight2[i,1:J.resight] <- (s[i,1]-x.resight[1:J.resight,1])^2 + (s[i,2]-x.resight[1:J.resight,2])^2
    lam[i,1:J.resight] <- lam0.resight * exp(-d.resight2[i,1:J.resight]/(2*sigma2))*z[i]
    
    # Proceso de marcaje
    for(j in 1:J.mark) {
      y.mark[i,j] ~ dbinom(lam.mark[i,j], Oper.mark.R[j])
    }
  }
   # Telemetría para los 4 animales con GPS
  for(r in 1:n.locs){
    locs[r,1]~dnorm(s[ind[r],1], 1/(sigma2))
    locs[r,2]~dnorm(s[ind[r],2], 1/(sigma2))
  }

  # Proceso de reavistamiento
  for(i in 1:n.marked)  {
    for(j in 1:J.resight) {
      for(k in 1:K){
        # Modelo de historiales completos de los marcados. Atención a 
        # la máscara de fechas de marcaje
        y.full[i,j,k] ~ dpois(lam[i,j]*mark.date[i,k]*Oper.resight[j,k])
        # Observaciones reconocidas de marcados con probabilidad id.prob
        y.obs[i,j,k] ~ dbin(id.prob, y.full[i,j,k])
      }
    }
  }

  # Conteos de individuos no marcados
  for(j in 1:J.resight) {
    for(k in 1:K){
      # Atención a la máscara de fechas de marcaje
      Lam[j,k] <- sum(lam[1:M,j]*(1-mark.date[1:M,k]))
      n[j,k] ~ dpois(Lam[j,k]*Oper.resight[j,k])
    }
  }

  # Los datos de individuos marcados no identificados se usan en 
  # el muestreador MH-Hasting

  N <- sum(z[1:M])
  D<-N/area
 })


## Muestreador MH para actualizar y[1:n.marked,j,k] de forma que 
## sume siempre nnid[j,k]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Definimos componentes
    nnid<-control$nnid
    j<-control$j
    k<-control$k
    n.marked<-control$n.marked
    #nnid<-control$nnid
    calcNodes <- model$getDependencies(target)
  },

  run = function() {
    # conteos esperados de individuos por trampa
    lam.curr <- model$lam[1:n.marked,j]

    # Muestreo de y[1:n.marked,j] reasignando n[j] usando la
    # probabilidad condicional completa
    switch.probs <- lam.curr[1:n.marked]/sum(lam.curr[1:n.marked])

    # propone nueva identificaciones para nnid[j,k]
    y.latent.curr <- model$y.full[1:n.marked,j,k] - model$y.obs[1:n.marked,j,k]
    y.latent.prop <- rmulti(1, nnid, switch.probs[1:n.marked])
    model$y.full[1:n.marked,j,k] <<- model$y.obs[1:n.marked,j,k] + y.latent.prop

    # Modelo inicial logProb
    model_lp_initial <- model$getLogProb(calcNodes)

    # Modelo propuesto logProb
    model_lp_proposed <- model$calculate(calcNodes)

    # log-Metropolis-Hastings ratio
    log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnid, 
                                                switch.probs, log=TRUE)) -
                    (model_lp_initial + dmulti(y.latent.prop,  nnid, 
                                               switch.probs, log=TRUE))

    # Pasos Metropolis-Hastings
    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, 
           logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, 
           logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

## CONSTANTRS
constants      <- list(nnid=nnid,
                       J.resight=J.resight,
                       J.mark=J.mark,
                       M=M,
                       K=K,
                       xlim=xlims,
                       ylim=ylims,
                       n.marked=n.marked,
                       n.collar=4,
                       n.locs=nlocs,
                       ind=ind,
                       area=area)
str(constants)

## DATOS
data     <-    list   (y.obs=array(0,c(n.marked,J.resight,K)), # None ID!!!
                       y.mark=y.mark.aug,
                       n=n,
                       locs=locs,
                       x.resight=x.resight,
                       x.mark=x.mark,
                       mark.date=mark.date.aug,
                       Oper.resight=Oper.resight,
                       Oper.mark.R=Oper.mark.R)
                
str(data)

## INICIOS
# Marcados ID
s.start <- cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1], ylims[2]))
d <- e2dist(s.start[1:n.marked,], x.resight)
lam0s<- 0.1
sigs <- 0.8
lam <- lam0s * exp( -(d^2)/(2 * sigs^2))
lam[lam=="NaN"]<-0

# Marcados no ID
yi <- array(0, c(n.marked, J.resight, K)) # matriz de reavistamientos
for (j in 1:J.resight) {
  for (k in 21:K) {  # empezamos con el dia 21, con algunos marcados
    if (nnid[j, k] > 0) {
      probs <- lam[ ,j]
      probs <- probs / sum(probs)
      latent.id <- sample(1:n.marked, nnid[j,k], prob = probs, replace = FALSE)
      yi[latent.id , j, k] <- 1
    }
  } # end of k
}   # end of j
datYknown<-0
yi <- yi + datYknown
apply(yi,1,sum)  # para evitar errores en inicios, testamos si: 
                 # sum(yi)==sum(apply(yi,c(1,3),sum)*mark.date)
sum(yi)
sum(apply(yi,c(1,3),sum)*mark.date)

zst<-rep(1,M)
id.prob.s<-sum(datYknown)/(sum(datYknown)+sum(nnid))
n.collar<-4

inits   <-        list(z=zst,
                       s=s.start,
                       lam0.resight=0.1,
                       lam0.mark=0.002,
                       sigma=sigs,
                       id.prob=id.prob.s,
                       psi=0.2,
                       y.full=yi)
str(inits)

params <- c('psi','lam0.resight','lam0.mark','sigma','N','D','id.prob')


# Compilamos y ejecutamos el modelo
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel,
                      constants=constants,
                      data=data,
                      inits=inits,
                      check=FALSE)
Rmodel$initializeInfo()
Cmodel <- compileNimble(Rmodel)
conf <- configureMCMC(Rmodel,monitors=params, thin=1, useConjugacy = TRUE)
#### --- Usamos ahora el MH  --- ####
conf$removeSampler("y.full")
for(j in 1:J.resight){
  for(k in 1:K){
    conf$addSampler(target = paste(paste("y.full[1:",n.marked,",",j,",",k,"]"), sep=""),
                    type = 'IDSampler',
                    control = list(nnid = nnid[j,k], j=j,k=k,n.marked=n.marked),
                    silent = TRUE)
  }
}

Rmcmc <- buildMCMC(conf)


Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
nb = 1000
ni = 5000 + nb
nc = 3

outNim <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)


options(scipen=999)
summary(outNim)
xyplot(outNim)
