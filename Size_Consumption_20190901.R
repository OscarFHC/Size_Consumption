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
# zp_sum <- zp_m %>%
#   arrange(Cruise, Station, rep) %>%
#   group_by(Cruise, Station, zp_taxa) %>%
#   summarize(Den = mean(Den),
#             Biom_ind = mean(Biom_ind),
#             Biom_all = mean(Biom_all)) %>%
#   mutate(ID = paste0(as.character(Cruise), "_", as.character(Station)))

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




