logit <- function(x) {log(x/(1-x))}

SimSCR0<-function (N = 100, K = 20, p0=0.1, sigma=0.5, discard0 = TRUE,
  array3d = FALSE, rnd = 2013){
    set.seed(rnd)
    traplocs <- cbind(sort(rep(1:10, 10)), rep(1:10, 10))
    Dmat <- e2dist(traplocs, traplocs)
    ntraps <- nrow(traplocs)
    delta <- 2.5*sigma
    Xl <- min(traplocs[, 1] - delta)
    Xu <- max(traplocs[, 1] + delta)
    Yl <- min(traplocs[, 2] - delta)
    Yu <- max(traplocs[, 2] + delta)
    sx <- runif(N, Xl, Xu)
    sy <- runif(N, Yl, Yu)
    S <- cbind(sx, sy)
    D <- e2dist(S, traplocs)
    alpha0 <- logit(p0)
    sigma <- 0.5
    alpha1 <- 1/(2 * sigma * sigma)
    probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
    Y <- matrix(NA, nrow = N, ncol = ntraps)
    for (i in 1:nrow(Y)) {
        Y[i, ] <- rbinom(ntraps, K, probcap[i, ])
    }
    if (discard0) {
        totalcaps <- apply(Y, 1, sum)
        Y <- Y[totalcaps > 0, ]
    }
    dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), sep = ""))
    if (array3d) {
        Y <- array(NA, dim = c(N, ntraps, K))
        for (i in 1:nrow(Y)) {
            for (j in 1:ntraps) {
                Y[i, j, 1:K] <- rbinom(K, 1, probcap[i, j])
            }
        }
        if (discard0) {
            Y2d <- apply(Y, c(1, 2), sum)
            ncaps <- apply(Y2d, 1, sum)
            Y <- Y[ncaps > 0, , ]
        }
    }

    list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl,
        Yu), N = N, p0=p0, alpha1 = alpha1, sigma = sigma,
        K = K, S=S)
}



spatial.plot<- function (x, y, add = FALSE, cx = 1, col = "gray")
{
    nc <- as.numeric(cut(y, 10))
    if (!add)
        plot(x, pch = " ", asp = 1)
    if (col == "gray") {
        cc <- seq(3, 17, , 10)/20
        cc <- gray(cc)
    }
    else cc <- terrain.colors(10)
    points(x, pch = 20, col = cc[nc], cex = cx)
    image.scale(y, col = cc)
}

#Hay 5 tipos de spiderplot:
#1. spiderplot: clasico
#2. spiderplotJJ: posibilidad de buffer alrededor de las trampas
#3. spiderplotJJ2: diferentes colores por animal, diferentes grosores en los segmentos y buffer.
#4. spiderplotJJ3: spiderplot a un plot existente, con diferentes colores, grosores y buffer
#5. spiderplotJJ4: spiderplot clasico a un plot existente

spiderplot<-function (y, traplocs)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy
    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ<-function (y, traplocs, buffer)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy

    Xl<-min(traplocs[,1])-buffer
    Xu<-max(traplocs[,1])+buffer
    Yl<-min(traplocs[,2])-buffer
    Yu<-max(traplocs[,2])+buffer
    xlim<-c(Xl,Xu)
    ylim<-c(Yl,Yu)

    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5, xlim=xlim, ylim=ylim)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ2<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
            las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
            cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ3<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ4<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = 1, lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = 'violet')
    }
    par(op)
}


make.grid<-function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
nx = 40, ny = NULL, buffer = 0)
  {
    if (is.null(ny))
    ny <- nx
    if (!is.na(ll)) {
    minx <- min(ll[, 1])
    maxx <- max(ll[, 1])
    miny <- min(ll[, 2])
    maxy <- max(ll[, 2])
    bx <- (maxx - minx) * buffer
    by <- (maxy - miny) * buffer
    minx <- minx - bx
    maxx <- maxx + bx
    miny <- miny - by
    maxy <- maxy + by
  }
  x <- sort(rep(seq(minx, maxx, , nx), ny))
  y <- rep(seq(maxy, miny, , ny), nx)
  cbind(x, y)
}


rot<-function (m)
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

make.scrFrame<-function (caphist, traps, indCovs = NULL, trapCovs = NULL, sigCovs = NULL,
    trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")
{
    if (any(is.null(caphist), is.null(traps)))
        stop("caphist and trap must be provided")
    if (!is.list(caphist))
        stop("caphist must be a list")
    n.sessions <- length(caphist)
    caphist.dimensions <- sapply(caphist, dim)
    if (nrow(caphist.dimensions) == 2)
        caphist.dimensions <- rbind(caphist.dimensions, 1)
    for (i in 1:n.sessions) {
        caphist[[i]] <- array(caphist[[i]], dim = caphist.dimensions[,
            i])
        all.zero <- apply(apply(caphist[[i]], c(1, 3), sum),
            1, sum)
        if (any(all.zero == 0)) {
            cat("At least one individual has an all-zero encounter history",
                fill = TRUE)
            cat("Make sure this is ok...", fill = TRUE)
        }
    }
    if (!is.null(indCovs)) {
        if (!is.list(indCovs))
            stop("indCovs must be a list")
        if (any(!sapply(indCovs, is.data.frame)))
            stop("indCovs must be a list of dataframes")
        if (length(indCovs) != length(caphist))
            stop("number of sessions in indCovs does not match caphist")
        check.dim <- sapply(indCovs, nrow)
        if (any(check.dim != caphist.dimensions[1, ]))
            stop("number of individuals in indCovs does not match caphist")
        if (!("rmv" %in% indCovs[[1]])) {
            for (i in 1:length(indCovs)) {
                indCovs[[i]]$removed <- dim(caphist[[i]])[3]
            }
        }
    }
    else {
        indCovs <- list()
        for (i in 1:length(caphist)) {
            indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],
                dim(caphist[[i]])[1]))
        }
    }
    if (!is.list(traps))
        stop("traps must be a list")
    if (length(traps) != length(caphist))
        stop("number of sessions in traps does not match caphist")
    check.dim <- sapply(traps, nrow)
    if (!all(check.dim == caphist.dimensions[2, ]))
        stop("number of traps does not match caphist")
    if (!is.null(trapCovs)) {
        if (!is.list(trapCovs))
            stop("trapCovs must be a list")
        if (any(!sapply(trapCovs, is.list)))
            stop("trapCovs must be a list of lists")
        if (any(!unlist(sapply(trapCovs, function(x) sapply(x,
            is.data.frame)))))
            stop("trapCovs must be a list of dataframes")
        if (length(trapCovs) != length(caphist))
            stop("number of sessions in trapCovs does not match caphist")
        check.dim <- lapply(trapCovs, function(x) sapply(x, nrow))
        for (i in 1:length(check.dim)) {
            if (!all(check.dim[[i]] == caphist.dimensions[2,
                i]))
                stop("number of traps does not match caphist")
        }
    }
    if (!is.null(sigCovs)) {
        if (nrow(sigCovs) != length(caphist))
            stop("number of rows in sigCovs does not match number of sessions")
        if (!"session" %in% colnames(sigCovs)) {
            sigCovs$session <- factor(1:n.sessions)
        }
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    else {
        sigCovs <- data.frame(session = factor(1:n.sessions))
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    if (!is.null(trapOperation)) {
        if (!is.list(trapOperation))
            stop("trapOperation must be a list")
        if (length(trapOperation) != length(caphist))
            stop("number of sessions in trapOperation does not match caphist")
        check.dim <- sapply(trapOperation, nrow)
        if (!all(check.dim == caphist.dimensions[2, ]))
            stop("number of traps does not match caphist")
    }
    max.dist <- NULL
    for (i in 1:length(caphist)) {
        for (j in 1:nrow(caphist[[i]])) {
            if (dim(caphist[[i]])[3] > 1) {
                where <- apply(caphist[[i]][j, , ], 1, sum) >
                  0
            }
            else {
                where <- caphist[[i]][j, , ] > 0
            }
            if (sum(where) > 1)
                max.dist <- c(max.dist, max(0, dist(traps[[i]][where,
                  c("X", "Y")]), na.rm = T))
        }
    }
    mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
    mdm <- max(max.dist, na.rm = T)
    if (!is.null(telemetry)) {
        if (!is.list(telemetry$fixfreq))
            stop("telemetry$fixfreq must be a list")
        fixfreq.dimensions <- sapply(telemetry$fixfreq, dim)
        if (nrow(fixfreq.dimensions) == 2)
            fixfreq.dimensions <- rbind(fixfreq.dimensions, 1)
        if (!is.null(telemetry$indCovs)) {
            if (!is.list(telemetry$indCovs))
                stop("telemetry$indCovs must be a list")
            if (any(!sapply(telemetry$indCovs, is.data.frame)))
                stop("telemetry$indCovs must be a list of dataframes")
            if (length(telemetry$indCovs) != length(telemetry$fixfreq))
                stop("number of sessions in telemetry$indCovs does not match telemetry$fixfreq")
            check.dim <- sapply(telemetry$indCovs, nrow)
            if (any(check.dim != fixfreq.dimensions[1, ]))
                stop("number of individuals in telemetry$indCovs does not match telemetry$fixfreq")
            if (any(!names(indCovs[[1]]) %in% c(names(telemetry$indCovs[[1]]),
                "removed")))
                stop("indCovs do not match between capture and telemetry data")
        }
        if (!is.null(telemetry$cap.tel)) {
            if (!is.list(telemetry$cap.tel))
                stop("telemetry$indCovs must be a list")
            warning("make sure captured individuals w/ collars sorted first!")
        }
    }
    if (!is.null(rsfDF)) {
        library(FNN)
        rsfCovs <- names(rsfDF[[1]][, -c(1:2), drop = F])
        if (is.null(trapCovs)) {
            trapCovs <- list()
            length(trapCovs) <- n.sessions
            for (s in 1:n.sessions) {
                trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                  c("X", "Y")], traps[[s]][, c("X",
                  "Y")], 1)$nn.index)
                trapCovs[[s]] <- list()
                length(trapCovs[[s]]) <- caphist.dimensions[3,
                  s]
                for (k in 1:caphist.dimensions[3, s]) {
                  trapCovs[[s]][[k]] <- data.frame(rsfDF[[s]][trap.grid,
                    rsfCovs])
                  names(trapCovs[[s]][[k]]) <- rsfCovs
                }
            }
        }
        else {
            for (s in 1:n.sessions) {
                if (any(!rsfCovs %in% trapCovs[[s]][[1]])) {
                  miss.rsfCovs <- rsfCovs[which(!rsfCovs %in%
                    trapCovs[[s]][[1]])]
                  trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                    c("X", "Y")], traps[[s]][, c("X",
                    "Y")], 1)$nn.index)
                  for (k in 1:caphist.dimensions[3, s]) {
                    newtrapCovs <- data.frame(rsfDF[[s]][trap.grid,
                      miss.rsfCovs])
                    names(newtrapCovs) <- miss.rsfCovs
                    trapCovs[[s]][[k]] <- data.frame(trapCovs[[s]][[k]],
                      newtrapCovs)
                  }
                }
            }
        }
    }
    scrFrame <- list(caphist = caphist, traps = traps, indCovs = indCovs,
        trapCovs = trapCovs, sigCovs = sigCovs, trapOperation = trapOperation,
        occasions = caphist.dimensions[3, ], type = type, mmdm = mmdm,
        mdm = mdm, telemetry = telemetry)
    class(scrFrame) <- "scrFrame"
    return(scrFrame)
}



e2dist <- function (x, y) {  # Function from scrbook package to calculate the distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


SCRdensity<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col = terrain.colors(ncolors))
    image.scale(Dn, col = terrain.colors(ncolors))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


# Con este comando sacamos un distribución real del espacio entre el eje de
# las X y el de las Y. Para que sea regular el pixel, lo ajusto con asratio
# plot(X, asp=1)
# asratio <- (Yu-Yl)/(Xu-Xl)
# nx<-75
# ny<-75* asratio

SCRdensityJJ1<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=gray.colors(ncolors, start=1, end=0),asp=1)
    image.scale(Dn, col=gray.colors(ncolors, start=1, end=0))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


sim.pID.Data <- function (N = N, K = K, sigma = sigma, lam0 = lam0,
    knownID = knownID, X = X, xlims = xlims, ylims = ylims,
    obsmod = c("pois", "bern"),nmarked = c("known", "unknown"), rat = 1,
    tel = 0, nlocs = 0)
{
    if (tel > knownID)
        stop("tel cannot be bigger than knownID")
    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)
    npts <- dim(X)[1]
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    D <- e2dist(S, X)
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    Y <- array(NA, c(N, npts, K))
    for (i in 1:N) {
        for (j in 1:npts) {
            if (identical(obsmod, "bern")) {
                Y[i, j, ] <- rbinom(K, 1, lam[i, j])
            }
            else if (identical(obsmod, "pois")) {
                Y[i, j, ] <- rpois(K, lam[i, j])
            }
        }
    }
    n <- apply(Y, c(2, 3), sum)
    Yknown <- Y[1:knownID, , ]
    if (identical(nmarked, "unknown")) {
        iobs <- which(apply(Yknown > 0, 1, any))
        Yobs <- Yknown[iobs, , ]
    }
    else if (identical(nmarked, "known")) {
        Yobs <- Yknown
    }
    YknownR <- Yobs
    counter <- array(0, c(dim(Yobs)[1], dim(X)[1], K))
    for (i in 1:dim(Yobs)[1]) {
        for (j in 1:dim(X)[1]) {
            for (k in 1:K) {
                if (identical(obsmod, "bern")) {
                  if (YknownR[i, j, k] == 1) {
                    IDed <- rbinom(1, 1, rat)
                    if (IDed == 0) {
                      YknownR[i, j, k] <- 0
                      counter[i, j, k] <- 1
                    }
                  }
                }
                else if (identical(obsmod, "Ypois")) {
                  if (Yobs[i, j, k] > 0) {
                    IDed <- sum(rbinom(Yobs[i, j, k], 1, rat))
                    YknownR[i, j, k] <- IDed
                    if (IDed != Yobs[i, j, k]) {
                      counter[i, j, k] <- Yobs[i, j, k] - IDed
                    }
                  }
                }
            }
        }
    }
    n <- n - apply(counter, 2:3, sum)
    if (tel > 0) {
        itel <- sort(sample(1:knownID, tel, replace = F))
        locs <- list()
        for (i in 1:tel) {
            lx <- rnorm(nlocs, S[itel[i], 1], sigma)
            ly <- rnorm(nlocs, S[itel[i], 2], sigma)
            locs[[i]] <- cbind(lx, ly)
        }
    }
    else {
        locs <- NULL
        itel <- NULL
    }
    list(n = n, Y = Y, Yknown = Yknown, Yobs = Yobs, YknownR = YknownR,
        counter = sum(counter), locs = locs, telID = itel, S=S)
}


