# number of prey species / size c;ass of prey 
w <- 100
  
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
R_opt <- (1/0.07)^3
gamma <- 2 # assymetrical hump-shaped curve

q <- 0

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


phyto_m <- (seq(from = 0.2, to = 50, length = w))^3
phyto_D <- 10^5 * phyto_m^(-3/4)
#plot(x = log(phyto_m^(1/3)), y = log(phyto_D))
zp_m <- (seq(from = 100, to = 150, length = 50))^3
zp_D <- 1


#plot(x = phyto_m^(1/3), y = attack)

Consump <- as.data.frame(matrix(0, length(zp_m), length(phyto_m)))

for (i in 1:length(zp_m)){
  
  intf <- Interfer_func(m_pd = zp_m[i])
  handle <- Handling_func(m_pd = zp_m[i], m_py = phyto_m, Temp = 25 + 273.15)
  pref <- Pref_func(m_pd = zp_m[i], m_py = phyto_m)
  attack <- Attack_func(m_pd = zp_m[i], m_py = phyto_m, Temp = 25 + 273.15)
  
  for (j in 1:length(phyto_m)){
    Consump[i, j] <- w * attack[j] * (phyto_D[j])^(1 + q) / 
                    ( 1 + intf * zp_D + w * handle[j] * (sum(phyto_D[j]) - phyto_D[j])^(1+q) )
  }  
}


plot(x = phyto_m^(1/3), y = Consump[1,])
plot(x = phyto_m^(1/3), y = colSums(Consump))

