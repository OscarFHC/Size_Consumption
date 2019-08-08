R_opt <- 100 # optimal bodymass ratio of prey
gamma <- 2 # assymetrical hump-shaped curve

b_0 <- 50 # attack rate constant
beta_i <- rnorm(1, mean = 0.47, sd = 0.04) #scaling exponent of how attack rate scales on biomass of prey 
beta_j <- rnorm(1, mean = 0.15, sd = 0.03) #scaling exponent of how attack rate scales on biomass of predator

h_0 <- 0.4
eta_i <- rnorm(1, mean = -0.48, sd = 0.03) #scaling exponent of how handling time scales on biomass of prey 
eta_j <- rnorm(1, mean = -0.66, sd = 0.02) #scaling exponent of how handling time scales on biomass of predator

Pref_func <- function(m_i, m_j){
  p <- ( (m_i / m_j * R_opt) * exp(1 - (m_i / m_j * R_opt)) ) ^ gamma
  return(p)
}

Attack_func <- function(m_i, m_j){
  atta <- b_0 * (m_i^beta_i) * (m_j^beta_j) * pref
  return(atta)
}

handling_func <- function(m_i, m_j){
  h_0 * (m_i^eta_i) * (m_j^eta_j)
}


phyto_m <- (seq(from = 0.2, to = 50, length = 100))^3
zp_m <- 100^3


pref <- Pref_func(m_i  = phyto_m, m_j = zp_m)
#plot(x = phyto_m^(1/3), y = pref)
attack <- Attack_func(m_i  = phyto_m, m_j = zp_m)
#plot(x = phyto_m^(1/3), y = attack)
handle <- handling_func(m_i  = phyto_m, m_j = zp_m)
#plot(x = phyto_m^(1/3), y = handle)

q <- 0
c <- rnorm(1, mean = 0.8, sd = 0.2)
phyto_D <- 10^5 * phyto_m^(3/4)
zp_D <- 5
Consump <- c()

for(i in 1:length(phyto_m)){
  w <- 1 / length(phyto_m)
  Consump <- c(Consump, 
    (w * attack[i] * (phyto_D[i])^(1 + q) ) / 
    (1 + c * zp_D + w * handle[i] * (sum(phyto_D[i]) - phyto_D[i])^(1+q) )
  )
}

plot(x = phyto_m^(1/3), y = Consump)
