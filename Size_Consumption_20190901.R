library(tidyverse)

se_func <- function(x){
  perm_se <- sd(replicate(5000, mean(sample(x, size = length(x), replace = TRUE))))
  perm_se
}

# Interference term
C_0 <- 1
phi <- 0.25 # scaling exponent

# Handling time
h_0 <- 1
eta_pd <- -0.75
eta_py <- 0.5
E_h <- 0.65

# Attack (Capture) rate
a_0 <- 1
alpha_pd <- 0.92
alpha_py <- 0.66
E_a <- mean(c(-0.46, -0.96)) # range from -0.46 to -0.96

T_0 <- 293.15 #(K)
k <- 8.617333262145 * (10^(-5))
R_opt <- 4103.13 #1000
# optimum ESD ratio from Hansen 1994 is 18:1
# mean(zp_list$Biomass) / (exp(-0.583 + 0.86 * log(4/3 * pi * ((mean(c(237,165,132,132,159,237,109,159))/18)/2)^3)) * 10^-6)
gamma <- 2 # assymetrical hump-shaped curve

q <- 0 # q = 0 as type II, q = 1 as type III

Interfer_func <- function(m_pd){
  c <- C_0 * (m_pd^phi) * (m_pd^phi)
  return(c)
}

Handling_func <- function(m_pd, m_py, Temp){
  h_0 * (m_pd^eta_pd) * (m_py^eta_py) * exp((E_h * (Temp - T_0)) / (k * Temp * T_0))
}

Pref_func <- function(m_pd, m_py){
  p <- ( (m_pd / (m_py * R_opt)) * exp(1 - (m_pd / (m_py * R_opt))) ) ^ gamma
  return(p)
}

Attack_func <- function(m_pd, m_py, Temp){
  atta <- a_0 * (m_pd^alpha_pd) * (m_py^alpha_py) * pref * exp((E_h * (Temp - T_0)) / (k * Temp * T_0))
  return(atta)
}


phyto_all <- read.table(file = "D:/Research/Size_Consumption/Phyto.csv", sep = ",", header = TRUE) %>%
  filter(state == "c0") %>%
  gather(key = phyto_class, value = Den, -c("Cruise", "Station", "state", "rep")) %>%
  mutate(phyto_class = as.numeric(gsub("X", "", phyto_class))) %>%
  mutate(Biom_ind = ifelse(phyto_class < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, # biomass in "ug"
                      ifelse(phyto_class > 20 & phyto_class < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, 0)),
         Biom_all = Den * Biom_ind) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)
  
# phyto_sum <- phyto_m %>%
#   arrange(Cruise, Station) %>%
#   group_by(Cruise, Station, phyto_class) %>%
#   summarize(Den = mean(Den),
#             Biom_ind = mean(Biom_ind),
#             Biom_all = mean(Biom_all)) %>%
#   mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep)))

taxa_name <- c("Calaniod", "Oithonid", "Corycaeid", "Oncaeid", "Harpacticoid", "CN", "ON", "HN", "other")
zp_list <- as.data.frame(matrix(0, length(taxa_name), 3)) %>%
  rename(zp_taxa = V1, Length = V2, Biomass = V3) %>%
  mutate(zp_taxa = taxa_name,
         Length = c(237, 165, 132, 132, 159, 237,	109, 159, 150), 
         a = c(-4.309, -5.3245, -7.458, -7.458, -6.4, -4.309, -5.3245, -6.4, -6.4),
         b = c(1.871, 1.997, 2.875, 2.875, 2.62, 1.871, 1.997, 2.62, 2.62)) %>%
  mutate(Biomass = 10^(a + log10(Length)*b)) # Biomass in "ug"

# zp_all.test <- read.table(file = "D:/Research/Size_Consumption/Zoopl.csv", sep = ",", header = TRUE) %>%
#   gather(key = zp_taxa, value = Den, -c("Cruise", "Station", "rep")) %>%
#   mutate(biomass = NA)
# 
# for (i in 1:nrow(zp_all.test)){
#   zp_all.test[i, ]$biomass <- zp_list[zp_list$zp_taxa == zp_all.test$zp_taxa[i], ]$Biomass
# }

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
  
# to check if zp is in phyto
phyto_sum <- phyto_all[which(phyto_all$ID %in% unique(zp_all$ID)),]
zp_sum <- zp_all[which(zp_all$ID %in% unique(phyto_all$ID)),]

#plot(x = phyto_m^(1/3), y = attack)

Consump_pred_w <- as.data.frame(matrix(0, length(unique(phyto_sum$ID)), length(unique(phyto_sum$phyto_class))))

for (k in 1:length(unique(phyto_sum$ID))){
  phyto <- phyto_sum[which(phyto_sum$ID == unique(phyto_sum$ID)[k]),]
  zp <- zp_sum[which(zp_sum$ID == unique(phyto_sum$ID)[k]),]
  for (i in 1:(length(zp$zp_taxa)-1)){ # "-1" because total zp is not used
    intf <- Interfer_func(m_pd = zp$Biom_ind[i])
    handle <- Handling_func(m_pd = zp$Biom_ind[i], m_py = phyto$Biom_ind, Temp = 25 + 273.15)
    pref <- Pref_func(m_pd = zp$Biom_ind[i], m_py = phyto$Biom_ind)
    attack <- Attack_func(m_pd = zp$Biom_ind[i], m_py = phyto$Biom_ind, Temp = 25 + 273.15)
    
    single_zp_consume <- c()
    for (j in 1:length(phyto$phyto_class)){
      w <- 1/length(phyto$phyto_class) 
      single_zp_consume <- c(single_zp_consume,
                             zp$Den[i] * (
                               w * attack[j] * (phyto$Den[j])^(1 + q) / 
                               ( 1 + intf * (zp$Den[i] - 1) + sum(w * handle * attack * (phyto$Den)^(1 + q)) )
                             )
                             )
    } 
    Consump_pred_w[k,] <- Consump_pred_w[k,] + single_zp_consume
  }
}

#plot(x = phyto$phyto_class, y = Consump_pred_w[k,])

Consump_pred <- as.data.frame(Consump_pred_w) %>%
  mutate(ID = unique(phyto_sum$ID)) %>%
  gather(key = cl, value = IR_pred, -ID) %>%
  mutate(phyto_class = ifelse(cl == "V1", 7.5,
                         ifelse(cl == "V2", 12.5,
                           ifelse(cl == "V3", 17.5,
                             ifelse(cl == "V4", 22.5, 
                               ifelse(cl == "V5", 27.5,
                                 ifelse(cl == "V6", 32.5,
                                   ifelse(cl == "V7", 37.5,
                                     ifelse(cl == "V8", 42.5,
                                       ifelse(cl == "V9", 47.5, 0)))))))))) %>%
  select(-cl)


#############################################################################
##### Empirical consumption #################################################
#############################################################################
phyto_ctrl <- read.table(file = "D:/Research/Size_Consumption/Phyto.csv", sep = ",", header = TRUE) %>%
  filter(state == "c24") %>%
  gather(key = phyto_class, value = Den, -c("Cruise", "Station", "state", "rep")) %>%
  mutate(phyto_class = as.numeric(gsub("X", "", phyto_class))) %>%
  mutate(Biom_ind = ifelse(phyto_class < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, # biomass in "ug"
                           ifelse(phyto_class > 20 & phyto_class < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, 0)),
         Biom_all = Den * Biom_ind) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)

phyto_trmt <- read.table(file = "D:/Research/Size_Consumption/Phyto.csv", sep = ",", header = TRUE) %>%
  filter(state == "t24") %>%
  gather(key = phyto_class, value = Den, -c("Cruise", "Station", "state", "rep")) %>%
  mutate(phyto_class = as.numeric(gsub("X", "", phyto_class))) %>%
  mutate(Biom_ind = ifelse(phyto_class < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, # biomass in "ug"
                           ifelse(phyto_class > 20 & phyto_class < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, 0)),
         Biom_all = Den * Biom_ind) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)

Consumption <- phyto_all %>%
  select(phyto_class, Cruise, Station, ID) %>%
  mutate(IR_C24 = log(phyto_ctrl$Biom_all)) %>% # from C24
  mutate(IR_T24 = log(phyto_trmt$Biom_all)) %>% # from T24
  mutate(IR_emp = log(phyto_ctrl$Biom_all/phyto_trmt$Biom_all)) %>%
  inner_join(Consump_pred, by = c("ID" = "ID", "phyto_class" = "phyto_class"))

Consumption_sum <- Consumption %>%
  filter(IR_emp != Inf & IR_emp != -Inf) %>%
  group_by(Cruise, Station, phyto_class) %>%
  summarize(meanIR_C24 = mean(IR_C24),
            meanIR_T24 = mean(IR_T24),
            meanIR_emp = mean(IR_emp),
            meanIR_pred = mean(IR_pred),
            seIR_C24 = se_func(IR_C24),
            seIR_T24 = se_func(IR_T24),
            seIR_emp = se_func(IR_emp),
            seIR_pred = se_func(IR_pred)
            )

Consumption_p <- Consumption_sum %>%
  gather(key = IRtype, value = IR, -c(Cruise, Station, phyto_class, seIR_C24, seIR_T24, seIR_emp, seIR_pred)) %>%
  gather(key = SEtype, value = SE, -c(Cruise, Station, phyto_class, IRtype, IR)) %>%
  mutate(CrSt = paste0(as.character(Cruise), "_", as.character(Station))) %>%
  ggplot(aes(x = phyto_class, y = IR, color = factor(IRtype, level = c("meanIR_pred", "meanIR_emp", "meanIR_C24", "meanIR_T24")))) +
    geom_point() + 
    scale_color_manual(values=c("#0072B2", "#D55E00", "#56B4E9", "#999999"), name = "") + 
    facet_grid(Cruise ~ Station)
Consumption_p
  
Consumption %>%
  ggplot(aes(x = phyto_class, y = IR_pred)) +
           geom_point()
###################################################################################################
##### ADBM model ##################################################################################
###################################################################################################
phyto_sum <- phyto_all[which(phyto_all$ID %in% unique(zp_all$ID)),]
zp_sum <- zp_all[which(zp_all$ID %in% unique(phyto_all$ID)),]

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
  
  E <- e*M#^ei
  
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

Get.web <- function(EHL, energy.intake = TRUE){
  
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



for (i in 1:length(unique(phyto_sum$ID))){
  phyto <- phyto_all[which(phyto_all$ID == unique(phyto_sum$ID)[i]), ]
  zp <- zp_all[which(zp_all$ID == unique(zp_sum$ID)[i]), ]

  e <- 1
  h <- 1
  h.a <- 10^-5
  h.b <- 0.4
  h.E <- 0.65 
  a <- 1
  ai <- 0.92
  aj <- 0.66
  a.E <- mean(c(-0.46, -0.96))
  temperature.C <- 25
  T0 <- 273.15
  
  N <- c(phyto$Den[order(phyto$Biom_ind)], zp[-nrow(zp),]$Biom_ind[order(zp[-nrow(zp),]$Biom_ind)])
  M <- c(sort(phyto$Biom_ind), sort(zp[-nrow(zp),]$Biom_ind))
  
  EHL <- Get.EHL.ratioH.temperature(M = M,
                                    e = e,
                                    h.a = h.a,
                                    h.b = h.b,
                                    h.E = h.E,
                                    a = a,
                                    ai = ai,
                                    aj = aj,
                                    a.E = a.E,
                                    temperature.C = temperature.C,
                                    N = N)
  
  wb <- Get.web(EHL = EHL)
  test1 <- rowSums(sweep(wb$per.species.flux, MARGIN = 2, N, FUN = "*"))
}

M <- c(seq(from = 1, to = 10, length = 10), seq(from = 10, to = 100, length = 10))
e <- 1
h.a <- 10^-5
h.b <- 0.4
a <- 4.9568 * 10^-8
ai <- -0.01
aj <- 0.58
#n, ni,
N <- M^(-3/4)
###################################################################################################
##### ADBM model ##################################################################################
###################################################################################################




