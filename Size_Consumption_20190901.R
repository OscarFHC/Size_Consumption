library(tidyverse)

# number of prey species / size class of prey 

  
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
R_opt <- 100 # optimum size ratio from Hansen 1994
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


phyto_m <- read.table(file = "D:/Research/Size_Consumption/Phyto.csv", sep = ",", header = TRUE) %>%
  filter(state == "c0") %>%
  gather(key = phyto_class, value = Den, -c("Cruise", "Station", "state")) %>%
  mutate(phyto_class = as.numeric(gsub("X", "", phyto_class))) %>%
  mutate(Biom_ind = ifelse(phyto_class < 20, exp(-0.583 + 0.86 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, # biomass in "ug"
                      ifelse(phyto_class > 20 & phyto_class < 50, exp(-0.665 + 0.939 * log(4/3 * pi * (phyto_class/2)^3)) * 10^-6, 0)),
         Biom_all = Den * Biom_ind)
phyto_sum <- phyto_m %>%
  arrange(Cruise, Station) %>%
  group_by(Cruise, Station, phyto_class) %>%
  summarize(Den = mean(Den),
            Biom_ind = mean(Biom_ind),
            Biom_all = mean(Biom_all)) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station)))

taxa_name <- c("Calaniod", "Oithonid", "Corycaeid", "Oncaeid", "Harpacticoid", "CN", "ON", "HN", "other")
zp_list <- as.data.frame(matrix(0, length(taxa_name), 3)) %>%
  rename(zp_taxa = V1, Length = V2, Biomass = V3) %>%
  mutate(zp_taxa = taxa_name,
         Length = c(237, 165, 132, 132, 159, 237,	109, 159, 150), 
         a = c(-4.309, -5.3245, -7.458, -7.458, -6.4, -4.309, -5.3245, -6.4, -6.4),
         b = c(1.871, 1.997, 2.875, 2.875, 2.62, 1.871, 1.997, 2.62, 2.62)) %>%
  mutate(Biomass = 10^(a + log10(Length)*b)) # Biomass in "ug"

zp_m <- read.table(file = "D:/Research/Size_Consumption/Zoopl.csv", sep = ",", header = TRUE) %>%
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
  mutate(Biom_all = Den * Biom_ind)

zp_sum <- zp_m %>%
  arrange(Cruise, Station, rep) %>%
  group_by(Cruise, Station, zp_taxa) %>%
  summarize(Den = mean(Den),
            Biom_ind = mean(Biom_ind),
            Biom_all = mean(Biom_all)) %>%
  mutate(ID = paste0(as.character(Cruise), "_", as.character(Station)))

#plot(x = phyto_m^(1/3), y = attack)

Consump <- as.data.frame(matrix(0, length(unique(phyto_sum$ID)), length(unique(phyto_sum$phyto_class))))

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
      w <- length(phyto$phyto_class) 
      single_zp_consume <- c(single_zp_consume,
                             w * attack[j] * (phyto$Den[j])^(1 + q) / 
                             ( 1 + intf * (zp$Den[i] - 1) + w * handle[j] * sum(attack * (phyto$Den)^(1 + q)) )
                             )
    } 
    Consump[k,] <- Consump[k,] + single_zp_consume
  }
}




plot(x = phyto$phyto_class, y = Consump[k,])


