library(tidyverse)

#########################################################################################
##### Loading data and universally used functions #######################################
#########################################################################################
se_func <- function(x){
  perm_se <- sd(replicate(5000, mean(sample(x, size = length(x), replace = TRUE))))
  perm_se
}

taxa_name <- c("Calaniod", "Oithonid", "Corycaeid", "Oncaeid", "Harpacticoid", "CN", "ON", "HN", "other")
zp_list <- as.data.frame(matrix(0, length(taxa_name), 3)) %>%
  rename(zp_taxa = V1, Length = V2, Biomass = V3) %>%
  mutate(zp_taxa = taxa_name,
         Length = c(237, 165, 132, 132, 159, 237,	109, 159, 150), 
         a = c(-4.309, -5.3245, -7.458, -7.458, -6.4, -4.309, -5.3245, -6.4, -6.4),
         b = c(1.871, 1.997, 2.875, 2.875, 2.62, 1.871, 1.997, 2.62, 2.62)) %>%
  mutate(Biomass = 10^(a + log10(Length)*b)) # Biomass in "ug"

zp_all <- read.table(file = "D:/Research/Size_Consumption/Zoopl.csv", sep = ",", header = TRUE) %>%
  gather(key = zp_taxa, value = Den, -c("Cruise", "Station", "rep")) %>%
  mutate(Biom_ind = ifelse(zp_taxa == "Calaniod", zp_list$Biomass[1], # biomass in "ug"
                           ifelse(zp_taxa == "Oithonid", zp_list$Biomass[2],
                                  ifelse(zp_taxa == "Corycaeid", zp_list$Biomass[3],
                                         ifelse(zp_taxa == "Oncaeid", zp_list$Biomass[4],
                                                ifelse(zp_taxa == "Harpacticoid", zp_list$Biomass[5],
                                                       ifelse(zp_taxa == "CN", zp_list$Biomass[6],
                                                              ifelse(zp_taxa == "ON", zp_list$Biomass[7],
                                                                     ifelse(zp_taxa == "HN", zp_list$Biomass[8],
                                                                            ifelse(zp_taxa == "other", zp_list$Biomass[9], 0)))))))))) %>%
  mutate(Biom_all = Den * Biom_ind) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)

phyto_all <- read.table(file = "D:/Research/Size_Consumption/Phyto_10.csv", sep = ",", header = TRUE) %>%
  gather(key = phyto_ESD, value = Den, -c("Cruise", "Station", "state", "rep")) %>%
  spread(state, Den) %>%
  mutate(phyto_ESD = as.numeric(gsub("X", "", phyto_ESD))) %>%
  mutate(Biom_ind = ifelse(phyto_ESD < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_ESD/2)^3)) * 10^-6, # biomass in "ug"
                           ifelse(phyto_ESD > 20 & phyto_ESD < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_ESD/2)^3)) * 10^-6, 0)),
         Biom_T0 = Biom_ind * c0,
         Biom_C24 = Biom_ind * c24,
         Biom_T24 = Biom_ind * t24) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)

# to check if zp is in phyto
phyto <- phyto_all[which(phyto_all$ID %in% unique(zp_all$ID)),]
zp <- zp_all[which(zp_all$ID %in% unique(phyto_all$ID)),] %>%
  filter(zp_taxa != "Total")

Consumption <- phyto %>%
  #select(phyto_ESD, Cruise, Station, rep, ID) %>%
  mutate(Gc = log(Biom_C24/Biom_T0)) %>% # from C24
  mutate(Gt = log(Biom_T24/Biom_T0)) %>% # from T24
  mutate(Gz = log(Biom_C24/Biom_T24))
#########################################################################################
##### Loading data and universally used functions #######################################
#########################################################################################

#########################################################################################
##### setup functions and parameters for fixed optimal PPMR model########################
#########################################################################################
Interfer_func <- function(m_pd){
  c <- c_0 * (m_pd^phi) * (m_pd^phi)
  return(c)
}

Handling_func <- function(m_py, m_pd, Temp){
  h_0 * (m_pd^h_pd) * (m_py^h_py) * exp((E_h * (Temp - T_0)) / (k * Temp * T_0))
}

Attack_func <- function(m_py, m_pd, Temp, R_opt, gamma){
  pref <- ( (m_pd / (m_py * R_opt)) * exp(1 - (m_pd / (m_py * R_opt))) ) ^ gamma
  atta <- a_0 * (m_pd^a_pd) * (m_py^a_py) * pref * exp((E_a * (Temp - T_0)) / (k * Temp * T_0))
  return(atta)
}

R_opt_func <- function(R){
  mean(zp_list$Biomass) / 
    (exp(-0.583 + 0.86 * log(4/3 * pi * ((mean(c(237,165,132,132,159,237,109,159))/R)/2)^3)) * 10^-6)  
}
#########################################################################################
##### setup functions and parameters for fixed optimal PPMR model########################
#########################################################################################
range <- data.frame(
  c_0_range = seq(from = 0, to  = 1, length = 500),
  phi_0_range = seq(from = -1, to  = 1, length = 500),
  
  a_0_range = seq(from = -50, to = 50, length = 500),
  a_py_range = seq(from = 0, to  = 1, length = 500),
  a_pd_range = seq(from = 0, to  = 3, length = 500),
  E_a_range = seq(from = 0, to  = 1, length = 500),
  R_opt_range = apply(matrix(seq(from = 0.1, to  = 100, length = 500), 1, 500), 1, FUN = R_opt_func),
  gamma_range = seq(from = 0.5, to  = 100, length = 500),
  
  h_0_range = seq(from = -50, to = 50, length = 500),
  h_py_range = seq(from = 0, to  = 1, length = 500),
  h_pd_range = seq(from = -3, to  = 0, length = 500),
  E_h_range = seq(from = -1, to  = 0, length = 500)
)

param_set <- list()
for(i in 1:ncol(range)){
  param_basic <- data.frame(
    c_0 = rep(1, 500),
    phi = rep(0.25, 500), # scaling exponent
    # Attack (Capture) rate
    a_0 = rep(1, 500),#mean(c(-28.13, -27.68))
    a_py = rep(mean(c(1/3, 2/3)), 500),
    a_pd = rep(1/4 + 2/3, 500),
    E_a = rep(0.65, 500),
    R_opt = rep(4103.13, 500), #1000
    # optimum ESD ratio from Hansen 1994 is 18:1
    # mean(zp_list$Biomass) / (exp(-0.583 + 0.86 * log(4/3 * pi * ((mean(c(237,165,132,132,159,237,109,159))/18)/2)^3)) * 10^-6)
    gamma = rep(2, 500), # assymetrical hump-shaped curve
    # Handling time
    h_0 = rep(1, 500),
    h_pd = rep(-0.75, 500),
    h_py = rep(0.5, 500),
    E_h = rep(-0.65, 500),
    
    T_0 = rep(273.15, 500), #(K)
    k = rep(8.617333262145 * (10^(-5)), 500),
    q = rep(0, 500) # q = 0 as type II, q = 1 as type III))
  )
  
  param_basic[,i] <- range[,i]
  param_set[[i]] <- param_basic
}
#param_set[[10]]
#########################################################################################
##### setup functions and parameters for fixed optimal PPMR model########################
#########################################################################################
#########################################################################################
##### fixed optimal PPMR model ##########################################################
#########################################################################################
param_RSS <- as.data.frame(matrix(0, nrow(range), ncol(range)))
colnames(param_RSS) <- substr(colnames(range), 1, nchar(colnames(range))-6) 
#t1 <- Sys.time()
for (i in 1:ncol(range)){
  param <- param_set[[i]]
  RSS <- c()
  for (j in 1:nrow(param)){
    # Interference term
    c_0 <- param$c_0[j]
    phi <- param$phi[j]
    
    a_0 <- param$a_0[j]
    a_py <- param$a_py[j]
    a_pd <- param$a_pd[j]
    E_a <-  param$E_a[j]
    R_opt <- param$R_opt[j]
    gamma <- param$gamma[j]
    
    h_0 <- param$h_0[j]
    h_py <- param$h_py[j]
    h_pd <- param$h_pd[j]
    E_h <- param$E_h[j]
    
    T_0 <- param$T_0[j]
    k <- param$k[j]
    q <- param$q[j]
    
    Consump_PPMR_w <- as.data.frame(matrix(0, length(unique(phyto$ID)), length(unique(phyto$phyto_ESD))))
    
    for (k in 1:length(unique(phyto$ID))){
      phyto_temp <- phyto[which(phyto$ID == unique(phyto$ID)[k]),]
      zp_temp <- zp[which(zp$ID == unique(phyto$ID)[k]),]
      
      intf <- Interfer_func(m_pd = zp_temp$Biom_ind)
      handle <- outer(phyto_temp$Biom_ind, zp_temp$Biom_ind, Handling_func, 
                      Temp = 25 + 273.15)
      attack <- outer(phyto_temp$Biom_ind, zp_temp$Biom_ind, Attack_func, 
                      Temp = 25 + 273.15, R_opt <- R_opt, gamma <- gamma)
      
      w <- 1/length(phyto_temp$phyto_ESD) 
      zp_consume <- 
        t(zp_temp$Den * # need to transpose the matrix b/c default matrix times/divides by vector is by "row"
            (t(w * attack * (phyto_temp$c0)^(1 + q)) /
               ( 1 + intf * (zp_temp$Den - 1) + colSums(w * handle * attack * (phyto_temp$c0)^(1 + q)) ) 
            ))
      Consump_PPMR_w[k,] <- rowSums(zp_consume)
    }
    
    PPMR_Diff <- as.data.frame(Consump_PPMR_w) %>%
      rename("7.5" = V1, "15" = V2, "25" = V3, "35" = V4, "45" = V5) %>%
      mutate(ID = unique(phyto$ID)) %>%
      gather(key = phyto_ESD, value = G_PPMR, -ID) %>%
      mutate(phyto_ESD = as.numeric(phyto_ESD)) %>%
      inner_join(Consumption, by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
      filter(c0>0 & c0>0 & t24>0) %>%
      filter(Gz != Inf & Gz != -Inf) %>%
      mutate(diff = (G_PPMR - Gz)^2)

      RSS <- c(RSS, sum(PPMR_Diff$diff))
  }
  param_RSS[,i] <- RSS
}
#t2 <- Sys.time()
#t2-t1
plot(x = range$E_a_range, y = param_RSS$E_a)

write.table(param_RSS, file = "D:/Research/Size_Consumption/Sensitivity/Param_RSS.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)


c_0 <- 1
phi <- 0.25 # scaling exponent

# Attack (Capture) rate
a_0 <- 1#mean(c(-28.13, -27.68))
a_py <- mean(c(1/3, 2/3))
a_pd <- 1/4 + 2/3
E_a <-  0.65
R_opt <- 18
# optimum ESD ratio from Hansen 1994 is 18:1
# mean(zp_list$Biomass) / (exp(-0.583 + 0.86 * log(4/3 * pi * ((mean(c(237,165,132,132,159,237,109,159))/18)/2)^3)) * 10^-6)
gamma <- 2 # assymetrical hump-shaped curve

# Handling time
h_0 <- 1
h_pd <- -0.75
h_py <- 0.5
E_h <- -0.65

T_0 <- 273.15 #(K)
k <- 8.617333262145 * (10^(-5))
q <- 0 # q = 0 as type II, q = 1 as type III

n.sim <- 50000
rand <- data.frame(
  c_0_range = runif(min = c_0 * 0.5, max = c_0 * 1.5, n = n.sim),
  phi_range = runif(min = phi * 0.5, max = phi * 1.5, n = n.sim),
  
  a_0_range = runif(min = -50, max = 50, n = n.sim),
  a_py_range = runif(min = a_py * 0.5, max = a_py * 1.5, n = n.sim),
  a_pd_range = runif(min = a_pd * 0.5, max = a_pd * 1.5, n = n.sim),
  E_a_range = runif(min = E_a * 0.5, max = E_a * 1.5, n = n.sim),
  R_opt_range = runif(min = R_opt * 0.5, max = R_opt * 1.5, n = n.sim),
  gamma_range = runif(min = gamma * 0.5, max = gamma * 1.5, n = n.sim),
  
  h_0_range = runif(min = -50, max = 50, n = n.sim),
  h_py_range = runif(min = h_py * 0.5, max = h_py * 1.5, n = n.sim),
  h_pd_range = runif(min = h_pd * 1.5, max = h_pd * 0.5, n = n.sim),
  E_h_range = runif(min = E_h * 1.5, max = E_h * 0.5, n = n.sim),
  
  T_0 = rep(T_0, n.sim), #(K)
  k = rep(k, n.sim),
  q = rep(q, n.sim), # q = 0 as type II, q = 1 as type III))
  
  RSS = rep(0, n.sim)
)

t1 <- Sys.time()
for (i in 1:nrow(rand)){
  # Interference term
  c_0 <- rand$c_0[i]
  phi <- rand$phi[i]
  
  a_0 <- rand$a_0[i]
  a_py <- rand$a_py[i]
  a_pd <- rand$a_pd[i]
  E_a <-  rand$E_a[i]
  R_opt <- R_opt_func(rand$R_opt_range[i])
  gamma <- rand$gamma[i]
  
  h_0 <- rand$h_0[i]
  h_py <- rand$h_py[i]
  h_pd <- rand$h_pd[i]
  E_h <- rand$E_h[i]
  
  T_0 <- rand$T_0[i]
  k <- rand$k[i]
  q <- rand$q[i]
  
  Consump_PPMR_w <- as.data.frame(matrix(0, length(unique(phyto$ID)), length(unique(phyto$phyto_ESD))))
  
  for (k in 1:length(unique(phyto$ID))){
    phyto_temp <- phyto[which(phyto$ID == unique(phyto$ID)[k]),]
    zp_temp <- zp[which(zp$ID == unique(phyto$ID)[k]),]
    
    intf <- Interfer_func(m_pd = zp_temp$Biom_ind)
    handle <- outer(phyto_temp$Biom_ind, zp_temp$Biom_ind, Handling_func, 
                    Temp = 25 + 273.15)
    attack <- outer(phyto_temp$Biom_ind, zp_temp$Biom_ind, Attack_func, 
                    Temp = 25 + 273.15, R_opt <- R_opt, gamma <- gamma)
    
    w <- 1/length(phyto_temp$phyto_ESD) 
    zp_consume <- 
      t(zp_temp$Den * # need to transpose the matrix b/c default matrix times/divides by vector is by "row"
          (t(w * attack * (phyto_temp$c0)^(1 + q)) /
             ( 1 + intf * (zp_temp$Den - 1) + colSums(w * handle * attack * (phyto_temp$c0)^(1 + q)) ) 
          ))
    Consump_PPMR_w[k,] <- rowSums(zp_consume)
  }
  
  PPMR_Diff <- as.data.frame(Consump_PPMR_w) %>%
    rename("7.5" = V1, "15" = V2, "25" = V3, "35" = V4, "45" = V5) %>%
    mutate(ID = unique(phyto$ID)) %>%
    gather(key = phyto_ESD, value = G_PPMR, -ID) %>%
    mutate(phyto_ESD = as.numeric(phyto_ESD)) %>%
    inner_join(Consumption, by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
    filter(c0>0 & c0>0 & t24>0) %>%
    filter(Gz != Inf & Gz != -Inf) %>%
    mutate(diff = (G_PPMR - Gz)^2)
  
  rand$RSS[i] <- sum(PPMR_Diff$diff)
}
t2 <- Sys.time()
t2 - t1

rand[which(rand$RSS == min(rand$RSS)),]
plot(RSS ~ c_0_range, rand[which(rand$RSS < 300),])
plot(RSS ~ phi_range, rand[which(rand$RSS < 300),])
plot(RSS ~ a_0_range, rand[which(rand$RSS < 300),])
plot(RSS ~ a_py_range, rand[which(rand$RSS < 300),])
plot(RSS ~ a_pd_range, rand[which(rand$RSS < 300),])



#########################################################################################
##### fixed optimal PPMR model ##########################################################
#########################################################################################




#########################################################################################
##### setup functions and parameters for ADBM model######################################
#########################################################################################
Get.EHL.ratioH.temperature <- function(M,
                                       e, ei,
                                       h.a, h.b, h.E,
                                       a_0, a_py, a_pd, a.E,
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
  H <- outer(M, M, 
             H.f, kT, h.a, h.b, h.E)
  
  A.f <- function(Mi, Mj, kT, a_0, a_py, a_pd, a.E)
    A <- a_0 * Mi^a_py * Mj^a_pd * exp(a.E * kT)
  A <- outer(M, M, A.f, kT,
             a_0, a_py, a_pd, a.E)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A * n * Mi^ni
  
  EHL <- list(E, H, L)
  EHL
}

Get.EHL.powerH.temperature <- function(M,
                                       e, ei,
                                       h_0, h_py, h_pd, h.E,
                                       a_0, a_py, a_pd, a.E,
                                       n, ni,
                                       N,
                                       temperature.C) {
  
  tempK <- 273.15 + temperature.C
  kT <- (tempK-T0)/(k*tempK*T0)  ## used in the Arrhenius equation later
  
  E <- e*M^ei
  
  #N <- n*M^ni
  
  H.f <- function(Mi, Mj, kT, h_0, h_py, h_pd, h.E){
    H <- h_0 * Mi^h_py * Mj^h_pd * exp(h.E*kT)
  }
  H <- outer(M, M, 
             H.f, kT, h_0, h_py, h_pd, h.E)
  
  A.f <- function(Mi, Mj, kT, a_0, a_py, a_pd, a.E)
    A <- a_0 * Mi^a_py * Mj^a_pd * exp(a.E*kT)
  A <- outer(M, M, A.f, kT,
             a_0, a_py, a_pd, a.E)
  
  if(length(N)>1) {
    ##print(x)
    L <- A * N
  }
  if(length(N)==1)
    L <- A* n*Mi^ni
  
  EHL <- list(E, H, L)
  EHL
}

Get.web <- function(EHL, energy.intake = TRUE){
  
  ##E <- EHL[[1]]
  ##H <- EHL[[2]]
  ##L <- EHL[[3]]
  S <- length(EHL[[1]])
  
  web <- matrix(0, S, S)
  overall.energy <- numeric(S)
  per.species.energy <- matrix(0, S, S)
  per.species.consump <- matrix(0, S, S)
  
  ## in matrix P, columns are cosumers and contain profit of that consumer
  ## feeding on each prey (row)
  P <- EHL[[1]]/EHL[[2]]
  
  ## split code depending on whether encounter rates are predator specific or not
  if(is.matrix(EHL[[3]])){
    for(j in 1:S){
      if (sum(EHL[[3]][j,]) == 0) {
        web[,j] = 0
        overall.energy[j] = 0
        per.species.energy[,j] = 0
        per.species.consump[,j] = 0
      } else{
        ## ordering of p required
        p <- P[,j]
        order.by.p <- order(p, decreasing=T)
        p <- p[order.by.p]
        Lj <- EHL[[3]][,j][order.by.p]
        hj <- EHL[[2]][,j][order.by.p]
        Ej <- EHL[[1]][order.by.p]
        allprey <- ifelse(is.nan(Lj * hj), 0, Lj * hj)
        
        cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum(allprey))
        ##cumulative.profit[length(E)] <- NA
        
        if(!all(p==0)){
          web[,j] <- c(1, cumulative.profit[1:(length(EHL[[1]])-1)] < p[2:length(EHL[[1]])])[order(order.by.p)]
          overall.energy[j] <- cumulative.profit[sum(web[,j])] 
        }
        energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum(allprey)[sum(web[,j])]) 
        all.energies <- c(energies, rep(0, S-length(energies)))
        
        consump <- c(Lj)[1:sum(web[,j])] / (1 + cumsum(allprey)[sum(web[,j])]) 
        all.consump <- c(consump, rep(0, S-length(consump)))
        
        per.species.energy[,j] <- all.energies[order(order.by.p)]
        per.species.consump[,j] <- all.consump[order(order.by.p)] 
      }
    }
  }
  if(energy.intake)
    result <- list(web=web, 
                   overall.flux=overall.energy, per.species.flux=per.species.energy,
                   per.species.consump = per.species.consump)
  else
    result <- web
  result
}

e <- 1
ei <- 1
h.a <- 10^-5
h.b <- 0.4
h_0 <- 1
h_pd <- -0.75
h_py <- 0.5
h.E <- 0.65 
a_0 <- 1
a_pd <- 0.92
a_py <- 0.66
a.E <- mean(c(-0.46, -0.96))
k <- 8.617333262145 * (10^(-5))
temperature.C <- 25
T0 <- 273.15
#########################################################################################
##### setup functions and parameters for ADBM model######################################
#########################################################################################

#########################################################################################
##### ADBM model ########################################################################
#########################################################################################
zp <- zp %>% filter(zp_taxa != "Total")

wb_ratio <- list()
wb_pwr <- list()
Consump_ADBM <- as.data.frame(c(phyto$Biom_ind, zp$Biom_ind)) %>%
  mutate(ID = c(phyto$ID, zp$ID),
         TrophicLevel = c(rep("phyto", length(phyto$ID)), rep("zp", length(zp$ID))),
         G_ADBM_ratioH = 0,
         G_ADBM_pwrH = 0,
         phyto_ESD = c(phyto$phyto_ESD, zp$Biom_ind)) %>%
  rename(size = "c(phyto$Biom_ind, zp$Biom_ind)")

for (i in 1:length(unique(phyto$ID))){
  phyto_temp <- phyto[which(phyto$ID == unique(phyto$ID)[i]), ]
  zp_temp <- zp[which(zp$ID == unique(zp$ID)[i]), ]
  
  N <- c(phyto_temp$c0[order(phyto_temp$Biom_ind)], zp_temp$Den[order(zp_temp$Biom_ind)])
  M <- c(sort(phyto_temp$Biom_ind), sort(zp_temp$Biom_ind))
  
  EHL_ratio <- Get.EHL.ratioH.temperature(M = M,
                                          e = e, ei = ei,
                                          h.a = h.a, h.b = h.b, h.E = h.E,
                                          a_0 = a_0, a_py = a_py, a_pd = a_pd, a.E = a.E,
                                          temperature.C = temperature.C,
                                          N = N)
  EHL_pwr <- Get.EHL.powerH.temperature(M = M,
                                        e = e, ei = ei,
                                        h_0 = h_0, h_py = h_py, h_pd = h_pd, h.E = h.E,
                                        a_0 = a_0, a_py = a_py, a_pd = a_pd, a.E = a.E,
                                        n, ni,
                                        N,
                                        temperature.C)
  
  wb_ratio[[i]] <- Get.web(EHL = EHL_ratio)
  wb_pwr[[i]] <- Get.web(EHL = EHL_pwr)
  
  Consump_ADBM$G_ADBM_ratioH[which(Consump_ADBM$ID == unique(phyto_temp$ID))] <- 
    rowSums(wb_ratio[[i]]$per.species.consump)
  #rowSums(sweep(wb_ratio[[i]]$per.species.consump, MARGIN = 2, N, FUN = "*"))
  Consump_ADBM$G_ADBM_pwrH[which(Consump_ADBM$ID == unique(phyto_temp$ID))] <- 
    rowSums(wb_pwr[[i]]$per.species.consump)
  #rowSums(sweep(wb_pwr[[i]]$per.species.consump, MARGIN = 2, N, FUN = "*"))
}

#########################################################################################
##### ADBM model ########################################################################
#########################################################################################

#########################################################################################
##### Plotting empirical measuremens vs theoretical predictions #########################
#########################################################################################
Consumption_sum <- Consumption %>%
  inner_join(Consump_PPMR, by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
  inner_join(Consump_ADBM[which(Consump_ADBM$TrophicLevel == "phyto"),], by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
  filter(c0>0 & c0>0 & t24>0) %>%
  filter(Gz != Inf & Gz != -Inf) %>%
  group_by(Cruise, Station, phyto_ESD) %>%
  summarize(meanGc = mean(Gc),
            meanGt = mean(Gt),
            meanGz = mean(Gz),
            meanG_PPMR = mean(G_PPMR),
            meanG_ADBM_ratioH = mean(G_ADBM_ratioH),
            meanG_ADBM_pwrH = mean(G_ADBM_pwrH),
            seGc = se_func(Gc),
            seGt = se_func(Gt),
            seGz = se_func(Gz),
            seG_PPMR = se_func(G_PPMR),
            seG_ADBM_ratioH = se_func(G_ADBM_ratioH),
            seG_ADBM_pwrH = mean(G_ADBM_pwrH)
  )

Consumption_p <- Consumption_sum %>%
  select(-c(meanGc, seGc, meanGt, seGt, meanG_ADBM_pwrH, seG_ADBM_pwrH)) %>%
  gather(key = Gtype, value = G, -c(Cruise, Station, phyto_ESD, seGz, seG_PPMR, seG_ADBM_ratioH)) %>%
  gather(key = SEtype, value = SE, -c(Cruise, Station, phyto_ESD, Gtype, G)) %>%
  mutate(CrSt = paste0(as.character(Cruise), "_", as.character(Station))) %>%
  ggplot(aes(x = phyto_ESD, y = G, 
             color = factor(Gtype, level = c("meanGz", "meanG_PPMR", "meanG_ADBM_ratioH", "meanG_ADBM_pwrH", "meanGc", "meanGt")))) +
  geom_point() +
  geom_jitter(width = 0.5) +
  geom_errorbar(aes(x = phyto_ESD, ymin = G - SE, ymax = G + SE), width = 0.2, size = 0.5) +
  scale_color_manual(values=c("#000033", "#0072B2", "#D55E00", "#56B4E9", "#999999", "#666666"), name = "") +
  facet_grid(Cruise ~ Station) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
# theme(panel.border = element_blank(), 
#       panel.grid.major = element_blank(), 
#       panel.grid.minor = element_blank(), 
#       plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
#       axis.line = element_line(colour = "black"), 
#       axis.text = element_text(size = 20), 
#       axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Consumption_p
#########################################################################################
##### Plotting empirical measuremens vs theoretical predictions #########################
#########################################################################################

Consumption_Diff <- Consumption %>%
  inner_join(Consump_PPMR, by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
  inner_join(Consump_ADBM[which(Consump_ADBM$TrophicLevel == "phyto"),], by = c("ID" = "ID", "phyto_ESD" = "phyto_ESD")) %>%
  filter(c0>0 & c0>0 & t24>0) %>%
  filter(Gz != Inf & Gz != -Inf) %>%
  mutate(diff_PPMR = (G_PPMR - Gz)^2,
         diff_ADBM_ratioH = (G_ADBM_ratioH - Gz)^2)




# 
# Consumption %>%
#   ggplot(aes(x = phyto_class, y = IR_pred)) +
#   geom_point()



# Get.EHL.powerH.temperature <- function(M,
#                                        e, ei,
#                                        h, hi, hj, h.E,
#                                        a, ai, aj, a.E,
#                                        n, ni,
#                                        N,
#                                        temperature.C) {
#   
#   tempK <- 273.15 + temperature.C
#   kT <- (tempK-T0)/(k*tempK*T0)  ## used in the Arrhenius equation later
#   
#   E <- e*M^ei
#   
#   #N <- n*M^ni
#   
#   H.f <- function(Mi, Mj, kT, h, hi, hj, h.E)
#     H <- h * Mi^hi * Mj^hj * exp(h.E*kT)
#   H <- outer(M, M, H.f, kT,
#              h, hi,
#              hj, h.E)
#   
#   A.f <- function(Mi, Mj, kT, a, ai, aj, a.E)
#     A <- a * Mi^ai * Mj^aj * exp(a.E*kT)
#   A <- outer(M, M, A.f, kT,
#              a, ai, aj, a.E)
#   
#   if(length(N)>1) {
#     ##print(x)
#     L <- A * N
#   }
#   if(length(N)==1)
#     L <- A* n*Mi^ni
#   
#   EHL <- list(E, H, L)
#   EHL
# }
head(phyto_examine)
zeros <- read.table(file = "D:/Research/Size_Consumption/Phyto_10.csv", sep = ",", header = TRUE) %>%
  filter(X7.5 == 0 | X15 == 0 |	X25 == 0 | X35 == 0 |	X45 == 0) %>%
  mutate(ID = paste0(Cruise, "_", Station, "_", rep))
unique(zeros$ID)


