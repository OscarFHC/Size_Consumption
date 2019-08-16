
## Function returns the energy content, handling time, and encounter rate
## given body masses, and scaling exponents
## for the ratio model of handling times
## M: vector of species masses in ascending order
## e: scaling exponent of energy with mass
## h.a and h.b: parameters of the scaling of handling time with mass
## a, ai, aj: parameters of the scaling of attack rate with mass
## n, ni: parameters of the scaling of abundance with mass
## (i is resource, j is consumer)
## N can be a vector of abundances, in which case n and ni are not used.
Get.EHL.RatioH <- function(M,
                           e,
                           h.a, h.b,
                           a, ai, aj,
                           n, ni,
                           N){
  
  Mi = M
  Mj = M
  
  get.h <- function(Mi, Mj, h.a, h.b)
    ifelse((h.b-Mi/Mj)>0,
           h.a/(h.b-Mi/Mj),
           Inf)
  
  ## in matrix H resources are rows and consumers are columns
  if(!h.b==0)
    H <- outer(M, M, get.h, h.a, h.b)
  if(h.b==0)
    H = matrix(h.a, length(M),length(M))
  
  ## ENCOUNTER RATES: consumer - resource mass specific encounter rates
  get.a <- function(Mi, Mj,
                    a, ai, aj)
    a * Mi^ai * Mj^aj
  A <- outer(M, M, get.a,
             a=a, ai=ai, aj=aj)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A* n*Mi^ni
  
  
  if(sum(order(M)==1:length(M))!=length(M))
    stop("Body sizes not sorted")
  
  ## energy values
  E <- e*M
  
  list(E=E, H=H, L=L)
  
}




## Function returns the energy content, handling time, and encounter rate
## given body masses, and scaling exponents
## for the power model of handling times
## M: vector of species masses in ascending order
## e: scaling exponent of energy with mass
## h.a and h.b: parameters of the scaling of handling time with mass
## a, ai, aj: parameters of the scaling of attack rate with mass
## n, ni: parameters of the scaling of abundance with mass
## (i is resource, j is consumer)
## N can be a vector of abundances, in which case n and ni are not used.
Get.EHL.powerH <- function(M,
                           e,
                           h, hi, hj,
                           a, ai, aj,
                           n, ni,
                           N){
  
  Mi = M
  Mj = M
  
  
  ## HANDLING TIMES
  get.h <- function(Mi, Mj, h, hi, hj)
    h * Mi^hi * Mj^hj
  ## in matrix H resources are rows and consumers are columns
  H <- outer(M, M, get.h, h=h, hi=hi, hj=hj)
  
  ## ENCOUNTER RATES: consumer - resource mass specific encounter rates
  get.a <- function(Mi, Mj,
                    a, ai, aj)
    a * Mi^ai * Mj^aj
  A <- outer(M, M, get.a,
             a=a, ai=ai, aj=aj)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A* n*Mi^ni
  
  if(sum(order(M)==1:length(M))!=length(M))
    stop("Body sizes not sorted")
  
  ## energy values
  E <- e*M
  
  return = list(E=E, H=H, L=L)
  
}





Get.EHL.ratioH.temperature <- function(M,
                                       e, ei,
                                       h.a, h.b, h.E,
                                       a, ai, aj, a.E,
                                       n, ni,
                                       N,
                                       temperature.C) {
  
  tempK <- 273.15 + temperature.C
  kT <- (tempK-T0)/(k*tempK*T0)  ## used in the Arrhenius equation later
  
  E <- e*M^ei
  
  #N <- n*M^ni
  
  H.f <- function(Mi, Mj, kT, h.a, h.b, h.E)
    ifelse(h.a / (h.b - Mi/Mj) * exp(h.E*kT)>0,
           h.a / (h.b - Mi/Mj) * exp(h.E*kT),
           Inf)
  H <- outer(M, M, H.f, kT,
             h.a, h.b, h.E)
  
  A.f <- function(Mi, Mj, kT, a, ai, aj, a.E)
    A <- a * Mi^ai * Mj^aj * exp(a.E*kT)
  A <- outer(M, M, A.f, kT,
             a, ai, aj, a.E)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A* n*Mi^ni
  
  EHL <- list(E, H, L)
  EHL
}



Get.EHL.powerH.temperature <- function(M,
                                       e, ei,
                                       h, hi, hj, h.E,
                                       a, ai, aj, a.E,
                                       n, ni,
                                       N,
                                       temperature.C) {
  
  tempK <- 273.15 + temperature.C
  kT <- (tempK-T0)/(k*tempK*T0)  ## used in the Arrhenius equation later
  
  E <- e*M^ei
  
  #N <- n*M^ni
  
  H.f <- function(Mi, Mj, kT, h, hi, hj, h.E)
    H <- h * Mi^hi * Mj^hj * exp(h.E*kT)
  H <- outer(M, M, H.f, kT,
             h, hi,
             hj, h.E)
  
  A.f <- function(Mi, Mj, kT, a, ai, aj, a.E)
    A <- a * Mi^ai * Mj^aj * exp(a.E*kT)
  A <- outer(M, M, A.f, kT,
             a, ai, aj, a.E)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A* n*Mi^ni
  
  EHL <- list(E, H, L)
  EHL
}




## Returns the binary food web, and the energy intake rates of consumers, if requested
Get.web <- function(EHL, energy.intake=F){
  
  ##E <- EHL[[1]]
  ##H <- EHL[[2]]
  ##L <- EHL[[3]]
  S <- length(EHL[[1]])
  
  web <- matrix(0, S, S)
  overall.energy <- numeric(S)
  per.species.energy <- matrix(0, S, S)
  
  ## in matrix P, columns are cosumers and contain profit of that consumer
  ## feeding on each prey (row)
  P <- EHL[[1]]/EHL[[2]]
  
  ## split code depending on whether encounter rates are predator specific or not
  if(is.matrix(EHL[[3]])){
    for(j in 1:S){
      
      ## ordering of p required
      p <- P[,j]
      order.by.p <- order(p, decreasing=T)
      p <- p[order.by.p]
      Lj <- EHL[[3]][,j][order.by.p]
      hj <- EHL[[2]][,j][order.by.p]
      Ej <- EHL[[1]][order.by.p]
      
      cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
      ##cumulative.profit[length(E)] <- NA
      
      if(!all(p==0)){
        web[,j] <- c(1, cumulative.profit[1:(length(EHL[[1]])-1)] < p[2:length(EHL[[1]])])[order(order.by.p)]
        overall.energy[j] <- cumulative.profit[sum(web[,j])]
      }
      energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
      all.energies <- c(energies, rep(0, S-length(energies)))
      
      per.species.energy[,j] <- all.energies[order(order.by.p)]
      
    }
  }
  if(is.vector(EHL[[3]])){
    for(j in 1:S){
      
      ## ordering of p required
      p <- P[,j]
      order.by.p <- order(p, decreasing=T)
      p <- p[order.by.p]
      Lj <- EHL[[3]][order.by.p]
      hj <- EHL[[2]][,j][order.by.p]
      Ej <- EHL[[1]][order.by.p]
      
      cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
      
      dj <- max(which(cumulative.profit==max(cumulative.profit)))
      web[,j] <- c(rep(1, dj), rep(0, S-dj))[order(order.by.p)]
      
      overall.energy[j] <- cumulative.profit[sum(web[,j])]    
      
      energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
      all.energies <- c(energies, rep(0, S-length(energies)))
      
      per.species.energy[,j] <- all.energies[order(order.by.p)]
      
    }
  }
  
  ##web[,!these] <- 0
  
  if(energy.intake)
    result <- list(web=web, overall.flux=overall.energy, per.species.flux=per.species.energy)
  else
    result <- web
  result
}


## Function for easily plotting a food web matrix
Plot.matrix <- function(web, title=" ", point.cex=0.5, trait.cex=1,
                        diag.line=T, traits=F, by.consumer=T, axes.labels=F, sp.pt.ch=NA){
  
  S <- length(web[,1])
  
  ##point.cex <- 30/30
  ##trait.cex <- 30/30
  
  dimnames(web) <- list(1:S, 1:S)
  consumer <- rep(1:S, each=S)
  resource <- rep(1:S, S)
  web.list <- Matrix.to.list(web)
  par(xpd=T)
  plot(consumer, resource, pch=19, type="n", cex=0.1,
       ann=F, axes=F,
       xlim=c(1, S), ylim=c(1, S))
  if(length(traits)==1)
    points(web.list[,2], S+1-as.numeric(web.list[,1]),
           type="p", pch=19, cex=point.cex)
  if(length(traits)==length(web)){
    
    colours.to.use <- rev(heat.colors(30)[c(-1:-5, -26:-30)])
    ##colours.to.use <- rev(gray(1:30/30)[c(-1:-5, -26:-30)])
    
    if(by.consumer){
      integer.traits <- matrix(0, S, S)
      for(i in 1:S){
        traits.01 <- traits[,i]-min(traits[,i])
        traits.01 <- traits.01/max(traits.01)
        integer.traits[,i] <- round(traits.01*19)+1
        integer.traits[traits[,i]==0,i] = NaN
        
      }
    }
    
    if(!by.consumer){
      colours.to.use <- heat.colors(20)
      traits.01 <- traits-min(traits)
      traits.01 <- traits.01/max(traits.01)
      integer.traits <- round(traits.01*19)+1
    }
    
    if(point.cex>trait.cex){
      points(web.list[,2], S+1-as.numeric(web.list[,1]),
             type="p", pch=19, cex=point.cex, col="black")##colours.to.use[integer.traits])
      points(rep(1:S, each=S), rep(S:1, times=S), 
             pch=19, cex=trait.cex, col=colours.to.use[integer.traits])
    }
    if(point.cex<trait.cex){
      points(rep(1:S, each=S), rep(S:1, times=S), 
             pch=19, cex=trait.cex, col=colours.to.use[integer.traits])
      points(web.list[,2], S+1-as.numeric(web.list[,1]),
             type="p", pch=19, cex=point.cex, col="black")##colours.to.use[integer.traits])
    }
    
    if(!is.na(sp.pt.ch))
      points(web.list[,2], S+1-as.numeric(web.list[,1]),
             type="p", pch=sp.pt.ch, cex=point.cex, col="black")##colours.to.use[integer.traits])
    
    
  }
  par(xpd=F)
  mtext(side=3, text=title, font=2, line=2, cex=0.5)
  if(axes.labels){
    mtext(side=2, "Resource", line=0.5, cex=1)
    mtext(side=3, "Consumer", line=0.5, cex=1)
  }
  if(diag.line==T)
    lines(1:S, S:1, lty="dashed")
  box()
}



## takes a food web in matrix format and coverts it to list format
Matrix.to.list <- function(web.matrix, predator.first=TRUE){
  if(length(dimnames(web.matrix)[[1]])==length(web.matrix[1,]))
    species.names <- dimnames(web.matrix)[[1]]
  else
    species.names <- 1:length(web.matrix[,1])
  web.list <- matrix(0, sum(web.matrix), 2)
  counter <- 1
  for(i in 1:length(web.matrix[,1]))
    for(j in 1:length(web.matrix[,1]))
      if(web.matrix[i,j]==1){
        web.list[counter,] <- c(species.names[i],species.names[j])
        counter <- counter + 1
      }
  if(!predator.first)
    web.list <- cbind(web.list[,2], web.list[,1])
  web.list
}
