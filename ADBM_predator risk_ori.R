# Functions & packages ----------------------------------------------------
library(igraph)      # for clustering and modularity (see http://cneurocvs.rmki.kfki.hu/igraph/doc/R/modularity.html)
library(vegan)       # for nestedness (switch from using nestedtemp to nestednodf)

# Calculate OF diet matrix
Diet.OF = function(SP. = SP, basal. = basal.loc, A. = A, S. = S, HD. = HD){
  D = matrix(0, nrow(SP.), nrow(SP.)); colnames(D) = rownames(SP.); rownames(D) = rownames(SP.) # the 0/1 matrix indicating diet (col eats row).
  All.IR = matrix(0, nrow(SP.), nrow(SP.)) # col = focusing consumer (ID), row = intake rate of eating till that ith prey on the profitability ranking.
  All.PreyOrder = c() # col = focusing consumer (ID), row = resource ID in the order of profitability ranking.
  P. = S./HD.
  IR = rep(NA, nrow(SP.)) # a vector that stores the intake rate that's given by the final output diet.
  for(j in 1:ncol(D)){ # for each consumer species
    if(j %in% basal.) {
      All.IR[, j] = 0
      PreyOrder = rep(NA, nrow(SP.))
      All.PreyOrder = cbind(All.PreyOrder, PreyOrder) # col = focusing consumer (ID), row = prey order (i.e. 1st row is the ID of the most profitable prey for each consumer)
      next} # skip calculation of intake rate for basal species
    PreyOrder = order(P.[, j], decreasing = T)
    # PreyOrder = PreyOrder[-which(PreyOrder == j)] # inactivate this line to allow cannibalism.
    All.PreyOrder = cbind(All.PreyOrder, PreyOrder)
    All.IR[, j] = unlist(lapply(PreyOrder, function(x){
      to_i = PreyOrder[1:which(PreyOrder == x)]
      sum(S.[to_i, j] * A.[to_i, j] * SP.[to_i, 2])/(1 + sum(S.[to_i, j] * A.[to_i, j] * SP.[to_i, 2] * HD.[to_i, j]))
    } # end of the function within lapply
    ))
    max = PreyOrder[which.max(All.IR[, j])]
    to_max = PreyOrder[1:which(PreyOrder == max)]
    D[to_max, j] = 1
    IR[j] = sum(S.[to_max, j] * A.[to_max, j] * SP.[to_max, 2])/(1 + sum(S.[to_max, j] * A.[to_max, j] * SP.[to_max, 2] * HD.[to_max, j]))
  } # end of j for-loop
  return(list(D = D, IR = IR, All.IR = All.IR, All.PreyOrder = All.PreyOrder))
}

# Calculate the HD.constant that give OF web with desired connectance
FindHDconst.OF = function(SP. = SP, basal. = basal.loc, A. = A, S. = S, HD. = HD, tot.link){
  TestHD1 = data.frame(H=NA, C=NA) # for roughly going through different scales of HD.const.
  TestHD2 = data.frame(H=NA, C=NA) # for sequentially getting the more precise HD.const
  tot.link.wanted = tot.link
  
  for(H in 1:50){ # H increases, C will decreases.
    HD.const = 10^H
    TestHD1[H, "H"] = H
    TestHD1[H, "C"] = sum(Diet.OF(SP., basal., A., S., HD. = HD.const * HD.)$D) # sum of links
  } # end of H for-loop
  close.HD = TestHD1[min(which(TestHD1[, "C"] - tot.link.wanted < 0)), "H"] # the smallest HD.const exponent that gives fewer links than desired.
  
  for(i in c(1, 0.1, 0.01, 0.001)){
    HD2 = seq(close.HD, close.HD-i, -0.1*i) # go trough smaller scales of exponent to find the exact HD.const
    for(H in 1:length(HD2)){
      HD.const = 10^HD2[H]
      TestHD2[H, "H"] = HD2[H]
      TestHD2[H, "C"] = sum(Diet.OF(SP., basal., A., S., HD. = HD.const * HD.)$D)
    } # end of H for-loop
    close.HD = TestHD2[max(which(TestHD2[, "C"] - tot.link.wanted < 0)), "H"]
  } # end of i for-loop
  
  HD.const = 10^TestHD2[which.min(abs(TestHD2[, "C"] - tot.link.wanted)), "H"] # the HD.const that gives the closest sum of links as desired.
  return(HD.const)
}

# Given the intake threshold (T), risk-response shape (eP), and diet breadth changing direction (di), calculate the predation risk compromised diet matrix
Diet.PR = function(SP. = SP, basal. = basal.loc, A. = A, S. = S, HD. = HD, T. = T, eP. = eP, di. = di){
  # this function is Diet.OF dependent. Start from the result of OF
  OF = Diet.OF(SP., basal., A., S., HD.)
  D = OF$D
  All.IR = OF$All.IR
  All.PreyOrder = OF$All.PreyOrder
  IR.OF = OF$IR # take the intake rate of OF as the upper bound of intake rate (no predation risk)
  max.risk = sum(SP[-basal., 2]) # the total biomass of all consumers as maximum possible predation risk
  
  IR = rep(NA, nrow(SP.)) # a vector that stores the intake rate that's given by the final output diet.
  
  repeat{ # recursively update the diet and predation risk, till the diet that is balanced.
    pred.risk.D = t(SP.[, 2] * t(D)) # diet matrix mutiplied by predators' biomasses
    risk.index = apply(pred.risk.D, 1, sum)/max.risk # total biomass of consumers of the focusing species/total biomass of all consumers
    ### below is the line where eP is involved.
    New.T = IR.OF - (IR.OF - T.)*(risk.index)^eP. # the new energy intake threshold (taken risk into account). NA for basal species.
    New.D = D
    cons.ID = (1:nrow(SP.))[-basal.] # consumers' ID
    if(di. == "broaden"){closest.loc = unlist(lapply(cons.ID, function(x){max(which(All.IR[, x] >= New.T[x]))}))} # which prey (location in profitability order, in broadening direction) gives the closest-larger intake rate to the threshold.
    if(di. == "narrow"){closest.loc = unlist(lapply(cons.ID, function(x){min(which(All.IR[, x] >= New.T[x]))}))} # which prey (location in profitability order, in narrowing direction) gives the closest-larger intake rate to the threshold.
    if(di. == "either"){closest.loc = unlist(lapply(cons.ID, function(x){
      BvsN = c(max(which(All.IR[, x] >= New.T[x])), min(which(All.IR[, x] >= New.T[x]))) # the broadest- and narrowest-larger diet 
      return(BvsN[which.min(All.IR[BvsN, x] - New.T[x])])}))} # which prey (location in profitability order, in either direction) gives the closest-larger intake rate to the threshold.
    
    if(is.infinite(prod(closest.loc)) | length(closest.loc)<length(cons.ID)) {New.D = NA; IR = NA; break} # Both means there's Inf/-Inf in closest.loc, meaning that the OF intake rate (upper bound) is already lower than threshold (i.e. threshold too high), then break the loop and return NA.
    for(I in 1:length(cons.ID)){
      to_closest.ID = All.PreyOrder[1:closest.loc[I], cons.ID[I]]
      New.D[to_closest.ID, cons.ID[I]] = 1
      IR[cons.ID[I]] = sum(S.[to_closest.ID, cons.ID[I]] * A.[to_closest.ID, cons.ID[I]] * SP.[to_closest.ID, 2])/(1 + sum(S.[to_closest.ID, cons.ID[I]] * A.[to_closest.ID, cons.ID[I]] * SP.[to_closest.ID, 2] * HD.[to_closest.ID, cons.ID[I]]))
    } # end of I for-loop
    if(all(New.D == D)) {break} else {D = New.D} # recursively do the above - adding links <-> predation risk 
  } # end of repeat loop
  return(list(D = New.D, IR = IR))
} # end of Diet.PR function

# Calculate the HD.constant and eP combination that give PR web with desired connectance.
FindHDconsteP.PR = function(SP. = SP, basal. = basal.loc, A. = A, S. = S, HD. = HD, tot.link, eP. = eP, di. = di){
  repeat{ # a repeat loop that increase eP bit by bit if wanted connectance is unreached
    
    TestHD1 = data.frame(H=NA, C=NA) # for roughly going through different scales of HD.const.
    TestHD2 = data.frame(H=NA, C=NA) # for sequentially getting the more precise HD.const
    tot.link.wanted = tot.link
    
    for(H in 1:50){ # H increases, C wil decreases. When H becomes larger than a certain degree, Diet.PR will failed as T is higher than OF CR.
      HD.const = 10^H
      TestHD1[H, "H"] = H
      TestHD1[H, "C"] = sum(Diet.PR(HD. = HD.const * HD., eP. = eP., di. = di.)$D) # sum of links
    } # end of H for-loop
    close.HD = TestHD1[max(which(TestHD1[, "C"] - tot.link.wanted > 0)), "H"] # the largest HD.const exponent that still gives more links than desired.
    
    for(i in c(1, 0.1, 0.01, 0.001)){
      HD2 = seq(close.HD, close.HD+i, 0.1*i) # go trough smaller scales of exponent to find the exact HD.const
      for(H in 1:length(HD2)){
        HD.const = 10^HD2[H]
        TestHD2[H, "H"] = HD2[H]
        TestHD2[H, "C"] = sum(Diet.PR(HD. = HD.const * HD., eP. = eP., di. = di.)$D)
      } # end of H for-loop
      close.HD = TestHD2[max(which(TestHD2[, "C"] - tot.link.wanted > 0)), "H"]
    } # end of i for-loop
    
    HD.const = 10^TestHD2[which.min(abs(TestHD2[, "C"] - tot.link.wanted)), "H"] # the HD.const that gives the closest sum of links as desired.
    
    if(min(abs(TestHD2[, "C"] - tot.link.wanted), na.rm = T) > 0.1 * tot.link) { eP. = eP. + 0.25   # 0.1 * tot.link as the connectance tolerance.
    } else {break} # if diff in connectance is small enough, break the repeat loop; if not, increase eP (decrease model sensitivity to the risk)
    if(eP. > 3) stop("eP is larger than 3") # 3 as eP upper-bound.
  } # end of eP-increasing repeat loop
  return(list(HD.const = HD.const, eP = eP.))
}

# Produce random diet matrix as null model - use the shuffling method
Diet.null = function(SP. = SP, basal. = basal, A. = A, S. = S, HD. = HD, tot.link){
  repeat{
    D = matrix(0, nrow(SP.), nrow(SP.)); colnames(D) = rownames(SP.); rownames(D) = rownames(SP.) # the 0/1 matrix indicating diet (col eats row).
    D[, basal.] = 999 # mark basal species who has no prey so columns should be all 0.
    ShuffledLocation = which(D != 999) # Shuffle all the consumers
    D[sample(ShuffledLocation, tot.link)] = 1
    D[, basal.] = 0 # change basal species columns back to 0.
    if( prod(apply(D, 2, sum)[-basal.]) > 0 & prod(apply(D, 1, sum)[basal.]) > 0 ){break} # make sure no additional basal nor isolated basal nodes.
  }
  IR = (t(D*S.*A.)/(1+apply(D*S.*A.*(SP.[, 2])*HD., 2, sum))) %*% (SP.[, 2]) # a vector that stores the intake rate that's given by the final output diet.
  return(list(D = D, IR = IR))
}

# for HD.3 (ratio handling time) as there're many Inf handling time combination, which shouldn't happen even in shuffled case.
Diet.null2 = function(SP. = SP, basal. = basal, A. = A, S. = S, HD. = HD, tot.link){
  loop = 1
  repeat{
    loop = loop + 1
    D = matrix(0, nrow(SP.), nrow(SP.)); colnames(D) = rownames(SP.); rownames(D) = rownames(SP.) # the 0/1 matrix indicating diet (col eats row).
    D[, basal.] = 999 # mark basal species who has no prey so columns should be all 0.
    D[which(HD==Inf)] = 999 # mark resource with Inf handling time so should be 0.
    ShuffledLocation = which(D != 999) # Shuffle all the consumers
    D[sample(ShuffledLocation, tot.link)] = 1
    D[which(D==999)] = 0 # change 999 back to 0.
    if( prod(apply(D, 2, sum)[-basal.]) > 0 & prod(apply(D, 1, sum)[basal.]) > 0 ){break} # make sure no additional basal nor isolated basal nodes.
    if(loop > 500) break # For HD.2 it's possible this loop become endless, so add a break here and let it redraw the community later.
  }
  HD.[which(HD.==Inf)] = 0 # to be able to calculate IR, turned those Inf in HD to 0 as they're not going to be mutiplied anyway.
  IR = (t(D*S.*A.)/(1+apply(D*S.*A.*(SP.[, 2])*HD., 2, sum))) %*% (SP.[, 2]) # a vector that stores the intake rate that's given by the final output diet.
  return(list(D = D, IR = IR))
}

# Calculate species pair's similarity (in terms of shared resource and consumer)
Simi = function(DietMatrix){
  S = matrix(0, nrow(DietMatrix), nrow(DietMatrix))
  for(i in 1:nrow(DietMatrix)){
    for(j in 1:nrow(DietMatrix)){
      SharedR = (DietMatrix[, i] & DietMatrix[, j])+0
      TotalR = (DietMatrix[, i] | DietMatrix[, j])+0
      SharedC = (DietMatrix[i, ] & DietMatrix[j, ])+0
      TotalC = (DietMatrix[i, ] | DietMatrix[j, ])+0
      S[i, j] = (sum(SharedR) + sum(SharedC))/(sum(TotalR) + sum(TotalC))
    } # For isolated basal species, NaN is produced. But this situation should be skipped so no change is made here.
  }
  return(S)
}

# Calculate unweighted TL (prey-averaged) - see Williams & Martinez 2004
# and Omnivory index & incoherence index based on this TL.
Unw.TL.Omn.q = function(DietMatrix, Spe.Richness. = Spe.Richness, basal. = basal){
  TL1 = -t(DietMatrix)/apply(DietMatrix, 2, sum)
  TL1[basal.,] = 0
  diag(TL1) = 1
  TL2 = c(rep(1, Spe.Richness.))
  TL = solve(TL1, TL2)
  
  Omn = rep(NA, Spe.Richness.)
  for(s in 1:Spe.Richness.){
    if(s %in% basal.) next
    if(apply(DietMatrix, 2, sum)[s] == 1){ Omn[s] = 0
    }else{
      Omn[s] = sd(TL[which(DietMatrix[, s]==1)])
    }
  }
  
  Distance1 = matrix(TL, Spe.Richness., Spe.Richness., byrow = T)
  Distance2 = t(Distance1)
  Distance = DietMatrix * (Distance2 - Distance1)
  q = sd(Distance[which(Distance != 0)])
  
  return(list(TL = TL, Omn = Omn, q = q))
}

# A summary function calculate SE and other stuff (found online, not original work).
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}
#####
# Nestedness is evaluated using nestednodf() of the vegan package.
# Modularity is evaluated using a combination of modularity(), cluster_walktrap(), and graph.adjacency() of the igraph package.
# Clustering coefficient is evaluated using a combination of transitivity() and graph.adjacency() of the igraph package.
# Incoherence index and Similarity are evaluated using the Unw.TL.Omn.q() and Simi() functions above.
# All other structural properties can easily be calculated using basic R commands with output diet matrices.



# For generating OF and PR counterparts from a given empirical web ------------------------
# (Below is the code for generating 1 web each).

# NOTE that info from an empirical web must be loaded first from EmpiricalFoodWebInfo.RData
# web.name, D.Empi, basal.loc, PreferredRatio, Spe.Bodymass, Spe.Richness need to extracted and assigned as inputs below.

tot.link = sum(D.Empi) # the (desired) total number of links for generating counterparts.

repeat{ # use a repeat loop here to drop the iteration where a PR web couldn't be generated with given parameter combination

  ### Assign species biomasses and foraging parameters
  # Species biomass densities
  SB1 = 10^-4 # biomass scaling constant
  SB2 = 0.25  # biomass scaling exponent
  SB1.draw = rnorm(1, SB1, 0.1*SB1)
  SB2.draw = rnorm(1, SB2, 0.1*SB2)
  # Spe.Biomass = 10^-4 * Spe.Bodymass^0.25
  Spe.Biomass = SB1.draw * Spe.Bodymass^SB2.draw
  rm(SB1, SB2, SB1.draw, SB2.draw)
  
  # Summarising species body mass & biomass as a SP matrix is critical, as most of our self-defined functions take SP as an input.
  SP = matrix(NA, nrow = Spe.Richness, ncol = 2)
  colnames(SP) = c("Body-mass", "Biomass"); rownames(SP) = 1:nrow(SP)
  SP[, 1] = Spe.Bodymass
  SP[, 2] = Spe.Biomass
  
  # Search rates
  A.2D = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(A.2D) = rownames(SP); rownames(A.2D) = rownames(SP)
  A1 = 1.069 # search rate scaling constant
  A2 = -0.58 # search rate scaling exponent for consumer mass
  A3 = 0.21 # search rate scaling exponent for resource mass
  A1.draw = rnorm(1, A1, 0.1*A1)
  A2.draw = rnorm(1, A2, 0.1*abs(A2))
  A3.draw = rnorm(1, A3, 0.1*A3)
  for(i in 1:nrow(SP)){
    A.2D[, i] = A1.draw * SP[i, 1]^A2.draw * SP[, 1]^A3.draw
    # A.2D[, i] = 1.069 * SP[i, 1]^-0.58 * SP[, 1]^0.21
  }; rm(i, A1, A2, A3, A1.draw, A2.draw, A3.draw)
  
  # Attack success probabilities
  S.1 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(S.1) = rownames(SP); rownames(S.1) = rownames(SP)
  S1 = 0.25 # boundary of attack success rate slope term
  S2 = 0.33 # curvarture of attack success rate slope term
  S3 = 0.2 # broadness of attack success rate dome-shaped term (large = steep)
  S1.draw = rnorm(1, S1, 0.1*S1)
  S2.draw = rnorm(1, S2, 0.1*S2)
  S3.draw = rnorm(1, S3, 0.1*S3)
  for(i in 1:nrow(SP)){
    S.1[, i] = (1/(1 + S1.draw*exp(-SP[i, 1]^S2.draw))) * (1/(1+ log10((1/PreferredRatio)*SP[, 1]/SP[i, 1])^2))^S3.draw
    # S.1[, i] = (1/(1 + 0.25*exp(-SP[i, 1]^0.33))) * (1/(1+ log10((1/PreferredRatio)*SP[, 1]/SP[i, 1])^2))^0.2
  }; rm(i, S1, S2, S3, S1.draw, S2.draw, S3.draw)
  
  # Handling times
  HD.1 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(HD.1) = rownames(SP); rownames(HD.1) = rownames(SP)
  HD1 = 0.25 # scaling exponent for HD
  HD2 = 1 # broadness of HD U-shaped term (large = broad)
  HD1.draw = rnorm(1, HD1, 0.1*HD1)
  HD2.draw = rnorm(1, HD2, 0.1*HD2)
  for(i in 1:nrow(SP)){
    HD.1[, i] = SP[i, 1]^HD1.draw *
      (dnorm(log10(PreferredRatio), log10(PreferredRatio), HD2.draw)/dnorm(log10(SP[, 1]/SP[i, 1]), log10(PreferredRatio), HD2.draw))
    # HD.1[, i] = SP[i, 1]^0.25 *
    #   (dnorm(log10(PreferredRatio), log10(PreferredRatio), 1)/dnorm(log10(SP[, 1]/SP[i, 1]), log10(PreferredRatio), 1))
  }; rm(i, HD1, HD2, HD1.draw, HD2.draw)
  
  # Energy requirement threshold (T)
  B0 = 4.15*10^-8 # arbitrarily give a constant at the moment...but if choose appropriate value it matches intake rate pretty well.
  T1 = -0.25 # treshold scaling exponent
  B0.draw = rnorm(1, B0, 0.1*B0)
  T1.draw = rnorm(1, T1, 0.1*abs(T1))
  T = B0.draw * SP[, 1]^T1.draw
  rm(T1, T1.draw, B0.draw)
  
  ### Confirm the inputs
  A= A.2D
  S = S.1
  HD = HD.1
  di = "broaden" # or could be "narrow", "either", determining the which of the diet-broadening/-narrowing strategies is used.
  # eP is assigned below.

  ### Derive OF and PR counterparts, with connectance as close to the desired as possible
  rm(HDeP) # remove it from the last iteration
  eP = 1 # eP starts from 1, and may be increased by the following function when approaching wanted connectance.
  tryCatch({HDeP = FindHDconsteP.PR(basal. = basal.loc, tot.link = tot.link)}, error=function(e){}) # sometime it might be NA by chance.
  if(exists("HDeP")) {} else {next} # if HD.const.PR doesn't exist, skip to the next repeat iteration.
  New.eP = HDeP$eP # the eP that produces PR web with 
  HD.const.PR = HDeP$HD.const
  PR = Diet.PR(basal. = basal.loc, HD. = HD.const.PR * HD, eP. = New.eP) # The diet that balanced between OF and Predation Risk, and have close to set connectance.
  D.PR = PR$D # This is the outpout PR web (as a diet matrix)
  
  HD.const.OF = FindHDconst.OF(basal. = basal.loc, tot.link = sum(D.PR)) # or tot.link = tot.link
  OF = Diet.OF(basal. = basal.loc, HD. = HD.const.OF * HD)
  D.OF = OF$D # This is the outpout OF web (as a diet matrix)
  break # only when HD.const.PR exist this line will be gone through and break the repeat loop.
} # end of check-point (no NA HD.const.PR) repeat-loop
#####


# For generating OF, PR, and Null webs with a given scheme ------------------------

### Food-web background parameters (some schemes)
if(Scheme == 2){Spe.Richness = 150} else {Spe.Richness = 50} # 50, 150
if(Scheme == 3){Basal.proportion = 0.4} else {Basal.proportion = 0.1} # 0.1, 0.4
if(Scheme == 4){MassMean = 10^-10} else {MassMean = 10^-5} # 10^-5, 10^-10, roughly empirical range
if(Scheme == 5){MassSD = 10^3} else {MassSD = 10} # 10^1, 10^3, roughly empirical range
if(Scheme == 6){BioExp = 0} else {BioExp = 0.25} # 0.25, 0, Damuth's law vs constant.
Set.Connectance = 0.1 # not far from empirical

### Simulate synthetic community  

repeat{ # use a repeat loop here to get rid of, if happens, isolated nodes (OF shouldn't but PR may). Check the end of this loop.
  basal = sample(1:Spe.Richness, round(Basal.proportion*Spe.Richness)) # randomly picked basal species.
  PreferredRatio = 0.1 # 0.1 not far from empirical
  Spe.Bodymass = rlnorm(Spe.Richness, log(MassMean), log(MassSD)) # mean & sd not far from empirical. Use log not log10 for suit rlnorm.
  Spe.Biomass = 10^-4 * Spe.Bodymass^BioExp # 10^-4 not far from empirical
    
    
  ### Foraging model parameters (some other schemes)
  # Summarising species body mass & biomass as a SP matrix is critical, as most of our self-defined functions take SP as an input.
  SP = matrix(NA, nrow = Spe.Richness, ncol = 2)
  colnames(SP) = c("Body-mass", "Biomass"); rownames(SP) = 1:nrow(SP)
  SP[, 1] = Spe.Bodymass
  SP[, 2] = Spe.Biomass
    
  # Search rate matrix (cols eat rows)
  A.2D = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(A.2D) = rownames(SP); rownames(A.2D) = rownames(SP)
  for(i in 1:nrow(SP)){
    A.2D[, i] = 1.069 * SP[i, 1]^-0.58 * SP[, 1]^0.21
  }; rm(i)
    
  A.3D = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(A.3D) = rownames(SP); rownames(A.3D) = rownames(SP)
  for(i in 1:nrow(SP)){
    A.3D[, i] = 2.721 * SP[i, 1]^-0.37 * SP[, 1]^0.42
  }; rm(i)
    
  # Attack success rate matrix (cols eat rows)
  S.1 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(S.1) = rownames(SP); rownames(S.1) = rownames(SP)
  for(i in 1:nrow(SP)){
    S.1[, i] = (1/(1 + 0.25*exp(-SP[i, 1]^0.33))) * (1/(1+ log10((1/PreferredRatio)*SP[, 1]/SP[i, 1])^2))^0.2
  }; rm(i)
  
  S.2 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(S.2) = rownames(SP); rownames(S.2) = rownames(SP)
  for(i in 1:nrow(SP)){
    S.2[, i] = (1/(1+ log10((1/PreferredRatio)*SP[, 1]/SP[i, 1])^2))^0.2
  }; rm(i)
    
  S.3 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(S.2) = rownames(SP); rownames(S.2) = rownames(SP)
  for(i in 1:nrow(SP)){
    S.3[, i] = (1/(1+ log10((1/PreferredRatio)*SP[, 1]/SP[i, 1])^2))^1
  }; rm(i)
    
  S.4 = matrix(1 , nrow = nrow(SP), ncol = nrow(SP)); colnames(S.3) = rownames(SP); rownames(S.3) = rownames(SP)
  
  # Mass-specific handling time among species (cols eat rows) - HD.constant depends on the aimmed connectance.
  # Here produce HD of different scheme without * HD.const (or say HD.const = 1).
  # HD.2 with a broader U-Shape, HD.3 as ADBM's ratio handling time.
   HD.1 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(HD.1) = rownames(SP); rownames(HD.1) = rownames(SP)
   for(i in 1:nrow(SP)){
     HD.1[, i] = SP[i, 1]^0.25 *
       (dnorm(log10(PreferredRatio), log10(PreferredRatio), 1)/dnorm(log10(SP[, 1]/SP[i, 1]), log10(PreferredRatio), 1))
   }; rm(i)
   
   HD.2 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(HD.2) = rownames(SP); rownames(HD.2) = rownames(SP)
   for(i in 1:nrow(SP)){
     HD.2[, i] = SP[i, 1]^0.25 *
       (dnorm(log10(PreferredRatio), log10(PreferredRatio), 1.1)/dnorm(log10(SP[, 1]/SP[i, 1]), log10(PreferredRatio), 1.1))
   }; rm(i) # 1.1 makes handling time increase much slower when deviates from PreferredRatio
    
   HD.3 = matrix(NA, nrow = nrow(SP), ncol = nrow(SP)); colnames(HD.3) = rownames(SP); rownames(HD.3) = rownames(SP)
   B = 1 # upper limit of acceptable MR/MC ratio
   for(i in 1:nrow(SP)){
     HD.3[, i] = 1 / (B - SP[, 1]/SP[i, 1]) * (SP[i, 1]/SP[, 1])
   }; rm(i)
   HD.3[which(HD.3 < 0)] = Inf
   ### possible isolated nodes: the smallest one is a consumer then it may have no suitable prey.
    
   # Energy requirement threshold (T)
   B0 = 4.15*10^-8 # related to respiration
   T = B0 * SP[, 1]^-0.25
   
   ### Parameter schemes
   A = A.2D
   if(Scheme == 7) A = A.3D
   S = S.1
   if(Scheme == 8) S = S.2
   if(Scheme == 9) S = S.3
   if(Scheme == 10) S = S.4
   HD = HD.1
   if(Scheme == 11) HD = HD.2
   if(Scheme == 12) HD = HD.3
   eP = 1
   if(Scheme == 13) eP = 0.5
   if(Scheme == 14) eP = 2

   ### Wiring food webs
  rm(HD.const.PR) # remove it from the last iteration
  tryCatch({HD.const.PR = FindHDconst.PR(tot.link = Set.Connectance * Spe.Richness^2)}, error=function(e){})
  if(exists("HD.const.PR")) {} else {next} # if HD.const.PR doesn't exist, skip to the next repeat iteration.
  PR = Diet.PR(HD. = HD.const.PR * HD) # The diet that balanced between OF and Predation Risk, and have close to set connectance.
  D.PR = PR$D # This is the outpout PR web (as a diet matrix)
  
  HD.const.OF = FindHDconst.OF(tot.link = sum(D.PR))
  OF = Diet.OF(HD. = HD.const.OF * HD)
  D.OF = OF$D # This is the outpout OF web (as a diet matrix)
    
  null = Diet.null(tot.link = sum(D.PR))
  if(Scheme == 12){null = Diet.null2(tot.link = sum(D.PR))} # for ratio handling time scheme 12, use special function to deal with Inf.
  D.null = null$D # This is the outpout Null web (as a diet matrix)
  
  if(min(apply(D.OF[basal, ], 1, sum)) > 0 & min(apply(D.PR[basal, ], 1, sum)) > 0 & # no isolated basal nodes
     min(apply(D.OF[, -basal], 2, sum)) > 0 & min(apply(D.PR[, -basal], 2, sum)) > 0 & # no consumer without any prey
     min(apply(D.null[, -basal], 2, sum)) > 0 & min(apply(D.null[, -basal], 2, sum)) > 0) break 
  # for random web (null) these checkings are already embedded in the D.null function, but for scheme 14 with Diet.null2, it's possible that after simulate 500 time still come up with a web with isolated node - redraw community in this case. 
} # end of check-point (no isolated node) repeat-loop
#####