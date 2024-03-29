library(tidyverse)

#########################################################################################
##### Loading data ######################################################################
#########################################################################################
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

phyto_all <- read.table(file = "D:/Research/Size_Consumption/Phyto.csv", sep = ",", header = TRUE) %>%
  gather(key = phyto_ESD, value = Den, -c("Cruise", "Station", "state", "rep")) %>%
  spread(state, Den) %>%
  mutate(phyto_ESD = as.numeric(gsub("X", "", phyto_ESD))) %>%
  mutate(Biom_ind = ifelse(phyto_ESD < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_ESD/2)^3)) * 10^-6, # biomass in "ug"
                           ifelse(phyto_ESD > 20 & phyto_ESD < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_ESD/2)^3)) * 10^-6, 0)),
         Biom_T0 = Biom_ind * c0,
         Biom_C24 = Biom_ind * c24,
         Biom_T24 = Biom_ind * c24) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station), "_", as.character(rep))) %>%
  arrange(Cruise, Station, rep)

# to check if zp is in phyto
phyto <- phyto_all[which(phyto_all$ID %in% unique(zp_all$ID)),]
zp <- zp_all[which(zp_all$ID %in% unique(phyto_all$ID)),]

Consumption <- phyto %>%
  #select(phyto_ESD, Cruise, Station, rep, ID) %>%
  mutate(Gc = log(Biom_C24/Biom_T0)) %>% # from C24
  mutate(Gt = log(Biom_T24/Biom_T0)) %>% # from T24
  mutate(Gz = log(Biom_C24/Biom_T24))
#########################################################################################
##### Loading data ######################################################################
#########################################################################################

#########################################################################################
##### setup fixed PPMR functions and parameters #########################################
#########################################################################################
se_func <- function(x){
  perm_se <- sd(replicate(5000, mean(sample(x, size = length(x), replace = TRUE))))
  perm_se
}

# Interference term
C_0 <- 1
phi <- 0.25 # scaling exponent

# Handling time
h_0 <- 1
h_pd <- -0.75
h_py <- 0.5
E_h <- 0.65

# Attack (Capture) rate
a_0 <- 1
a_pd <- 0.92
a_py <- 0.66
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
  h_0 * (m_pd^h_pd) * (m_py^h_py) * exp((E_h * (Temp - T_0)) / (k * Temp * T_0))
}

Pref_func <- function(m_pd, m_py){
  p <- ( (m_pd / (m_py * R_opt)) * exp(1 - (m_pd / (m_py * R_opt))) ) ^ gamma
  return(p)
}

Attack_func <- function(m_pd, m_py, Temp){
  atta <- a_0 * (m_pd^a_pd) * (m_py^a_py) * pref * exp((E_h * (Temp - T_0)) / (k * Temp * T_0))
  return(atta)
}
#plot(x = phyto_m^(1/3), y = attack)
#########################################################################################
##### setup fixed PPMR functions and parameters #########################################
#########################################################################################

#########################################################################################
##### fixed optimal PPMR model ##########################################################
#########################################################################################
Consump_PPMR_w <- as.data.frame(matrix(0, length(unique(phyto$ID)), length(unique(phyto$phyto_ESD))))

for (k in 1:length(unique(phyto$ID))){
  phyto_temp <- phyto[which(phyto$ID == unique(phyto$ID)[k]),]
  zp_temp <- zp[which(zp$ID == unique(phyto$ID)[k]),]
  for (i in 1:(length(zp_temp$zp_taxa)-1)){ # "-1" because total zp is not used
    intf <- Interfer_func(m_pd = zp_temp$Biom_ind[i])
    handle <- Handling_func(m_pd = zp_temp$Biom_ind[i], m_py = phyto_temp$Biom_ind, Temp = 25 + 273.15)
    pref <- Pref_func(m_pd = zp_temp$Biom_ind[i], m_py = phyto_temp$Biom_ind)
    attack <- Attack_func(m_pd = zp_temp$Biom_ind[i], m_py = phyto_temp$Biom_ind, Temp = 25 + 273.15)
    
    single_zp_consume <- c()
    for (j in 1:length(phyto_temp$phyto_ESD)){
      w <- 1/length(phyto_temp$phyto_ESD) 
      single_zp_consume <- c(single_zp_consume,
                             zp_temp$Den[i] * (
                               w * attack[j] * (phyto_temp$c0[j])^(1 + q) / 
                               ( 1 + intf * (zp_temp$Den[i] - 1) + sum(w * handle * attack * (phyto_temp$c0)^(1 + q)) )
                             )
                             )
    } 
    Consump_PPMR_w[k,] <- Consump_PPMR_w[k,] + single_zp_consume
  }
}

#plot(x = phyto$phyto_class, y = Consump_pred_w[k,])

Consump_PPMR <- as.data.frame(Consump_PPMR_w) %>%
  mutate(ID = unique(phyto$ID)) %>%
  gather(key = cl, value = G_PPMR, -ID) %>%
  mutate(phyto_ESD = ifelse(cl == "V1", 7.5,
                       ifelse(cl == "V2", 12.5,
                         ifelse(cl == "V3", 17.5,
                           ifelse(cl == "V4", 22.5, 
                             ifelse(cl == "V5", 27.5,
                               ifelse(cl == "V6", 32.5,
                                 ifelse(cl == "V7", 37.5,
                                   ifelse(cl == "V8", 42.5,
                                     ifelse(cl == "V9", 47.5, 0)))))))))) %>%
  select(-cl)

#########################################################################################
##### fixed optimal PPMR model ##########################################################
#########################################################################################

#########################################################################################
##### ADBM model ########################################################################
#########################################################################################
zp <- zp %>% filter(zp_taxa != "Total")

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
  
  H.f <- function(Mi, Mj, kT, h_0, h_py, h_pd, h.E)
    H <- h_0 * Mi^h_py * Mj^h_pd * exp(h.E*kT)
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
k <- 8.6173 * 10^-5
temperature.C <- 25
T0 <- 273.15

wb_ratio <- list()
wb_pwr <- list()
Consump_ADBM <- as.data.frame(c(phyto$Biom_ind, zp$Biom_ind)) %>%
  mutate(ID = c(phyto$ID, zp$ID),
         TrophicLevel = c(rep("phyto", length(phyto$ID)), rep("zp", length(zp$ID))),
         G_ADBM_ratioH = 0,
         G_ADBM_pwrH = 0) %>%
  rename(size = "c(phyto$Biom_ind, zp$Biom_ind)")

for (i in 1:length(unique(phyto_sum$ID))){
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
  Consump_ADBM$IR_ADBM[which(Consump_ADBM$ID == phyto$ID)] <- 
    rowSums(sweep(wb[[i]]$per.species.consump, MARGIN = 2, N, FUN = "*"))
}

#########################################################################################
##### ADBM model ########################################################################
#########################################################################################

Consumption_sum <- Consumption %>%
  inner_join(Consump_PPMR, by = c("ID" = "ID", "phyto_class" = "phyto_class")) %>%
  inner_join(Consump_ADBM[which(Consump_ADBM$TrophicLevel == "phyto"),], by = c("ID" = "ID", "phyto_class" = "size")) %>%
  filter(IR_emp != Inf & IR_emp != -Inf) %>%
  group_by(Cruise, Station, phyto_class) %>%
  summarize(meanIR_C24 = mean(IR_C24),
            meanIR_T24 = mean(IR_T24),
            meanIR_emp = mean(IR_emp),
            meanIR_PPMR = mean(IR_PPMR),
            meanIR_ADBM = mean(IR_ADBM),
            seIR_C24 = se_func(IR_C24),
            seIR_T24 = se_func(IR_T24),
            seIR_emp = se_func(IR_emp),
            seIR_PPMR = se_func(IR_PPMR),
            seIR_ADBM = se_func(IR_ADBM)
  )

Consumption_p <- Consumption_sum %>%
  gather(key = IRtype, value = IR, -c(Cruise, Station, phyto_class, seIR_C24, seIR_T24, seIR_emp, seIR_PPMR, seIR_ADBM)) %>%
  gather(key = SEtype, value = SE, -c(Cruise, Station, phyto_class, IRtype, IR)) %>%
  mutate(CrSt = paste0(as.character(Cruise), "_", as.character(Station))) %>%
  ggplot(aes(x = phyto_class, y = IR, 
             color = factor(IRtype, level = c("meanIR_emp", "meanIR_PPMR", "meanIR_ADBM", "meanIR_C24", "meanIR_T24")))) +
  geom_point() +
  geom_jitter(width = 0.5) +
  geom_errorbar(aes(x = phyto_class, ymin = IR - SE, ymax = IR + SE), width = 0.2, size = 0.5) +
  scale_color_manual(values=c("#0072B2", "#D55E00", "#56B4E9", "#999999", "#666666"), name = "") +
  facet_grid(Cruise ~ Station)
Consumption_p
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



