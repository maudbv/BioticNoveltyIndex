# Functions to simulate and visualize the Biotic Novelty Index
# written by Maud Bernard-Verdier and Conrad Schittko, February 2020

# This code comprises 5 functions, making use of the function BNI.calc:

# 1. simul.pool
# Generates a random pool of species and traits,
# with predetermined proportion of species in the three introduction status cateogries:
# natives(I; for indigenous), archaeobiota (A) and neobiota (N),
# and a given scenario of trait differences between species.
#
# 2. simul.comm
# Generates a random sets of communities from a given species pool
# with a predetermined proportion of neobiota species
#
# 3. simul.BNI
# Wrapper function for simul.comm, to generate multiple sets
# of random communities and calculate associated BNI, Rao and BNIs indices
#
# 4. plot.sim.BNI
# plotting function for a full graphical representation of the output of simul.BNI,
# with proportion of neophyte as the x-axis.
#
# 5. plot.simBNI.dates
# plotting function for sensitivity test on residence time, 
# with the residence time relevant to each scenario as the x-axis.
#________________________________________________________________________#




#________________________________________________________________________#

#### 1. FUNCTION simul.pool : Species POOL simulations   ####
#________________________________________________________________________#
#
# INPUT ##
# nsp             total number of species in the speices pool
# p.status        vector of proportion for natives, archaeo- and neobiota 
#                 in the regional pool. 
#                 Default is: c(nat = 0.70, arch = 0.15, neo = 0.15)
# t.arrival       Vector or minimal and maximal residence time (in number of years)
#                 to be modelled for each category (native/archaeo/neobiota). 
#                 Default is : c(nat.min= 8518,nat.max= 8518, 
#                               arch.min = 2756,arch.max = 2756,
#                               neo.min = 1, neo.max = 526)
# mean.traits     vector of mean values for the normal distributions from which
#                 to sample the trait values for each of the three categories: c(I,A,N).
#                 Default is c(0,0,0)
# sd.traits       vector of standard deviation for the normal distributions from which 
#                 to sample the trait values for each of the three categories.
#                 Default is c(1,1,1)
#
# OUTPUT : a list of elements:
# species.data    data frame of simulated introduction status and trait data for the nsp species
# comMat          simulated community matrix of dimension (nsites, nsp)
# traits          matrix of simulated trait data (for use in the BNI.calc function)
# nspecies          store variable nsp
# p.status          store variable p.status,
# t.arrival         store variable t.arrival,
# mean.traits       store variable mean.traits,
# sd.traits         store variable sd.traits
#_______________________________________________________________________

simul.pool <- function(nsp = 250,
                       p.status = c(nat = 0.70, arch = 0.15, neo = 0.15),
                       t.arrival = c(nat.min= 8518,nat.max= 8518,
                                     arch.min = 2756, arch.max = 2756,
                                     neo.min = 1, neo.max = 526),
                       mean.traits = c(0,0,0),
                       sd.traits = c(1,1,1)) {
  
  # Prepare species data matrix
  species.data <- data.frame(species = paste("sp", 1:nsp, sep = ""))
  
  # Assign species status (1 = native, 2 = archeobiota, 3 = neobiota), with 70% of natives
  species.data$status  <- sample(c(1, 2, 3), nsp, replace = T, prob = p.status)
  nat = species.data$status == 1 
  arc = species.data$status == 2 
  neo = species.data$status == 3 
  
  # Introduction: Assign years since introduction (YSI) based on the status
  species.data$YSI <- NA
  
  # Assign introduction years using uniform distributions between min and max:
  
  ## Natives
  species.data$YSI[nat] <- ceiling(runif(
    sum(nat),
    min = t.arrival[1],
    max = t.arrival[2]
  ))
  
  ## Archaeobiota 
  species.data$YSI[arc] <- ceiling(runif(
    sum(arc),
    min = t.arrival[3],
    max = t.arrival[4]
  ))
  
  ## Neobiota
  species.data$YSI[neo] <- ceiling(runif(
    sum(neo),
    min = t.arrival[5],
    max = t.arrival[6]
  ))
  
  
  ## Simulate random traits for each introduction status:
  
  # Trait T1 may be different accross species:
  species.data$T1[nat] = rnorm(sum(nat), mean = mean.traits[1] , sd = sd.traits[1])
  species.data$T1[arc] = rnorm(sum(arc), mean = mean.traits[2] , sd = sd.traits[2])
  species.data$T1[neo] = rnorm(sum(neo), mean = mean.traits[3] , sd = sd.traits[3])
  
  # Traits T2 and T3 are all sampled from the same N(0,1) distribution. 
  # These traits are only added to illustrate the use of a multivariate index
  species.data$T2 = rnorm(nsp, mean = 0 , sd = 1)  
  species.data$T3 =  rnorm(nsp, mean = 0 , sd = 1)
  
  
  # Return a list of elements defining the pool of species:
  pool <-
    list(
      species.data = species.data,
      traits = species.data[, c("T1", "T2","T3")],
      nspecies = nsp,
      p.status = p.status,
      t.arrival = t.arrival,
      mean.traits = mean.traits,
      sd.traits = sd.traits
    )
  
  return(pool)
}

#_______________________________________________________________________#

#### 2. FUNCTION simul.comm : assemble one random community from a species pool ####

#_______________________________________________________________________#
#
## INPUT ##
# pool        matrix of the regional pool of species, or output from simul.pool function
# SR          expected number of species in the community
# p.types     numeric vector of proportions for species in the three categories (I/A/N) 
#             Default is c(0.5, 0.25,0.25),
# simul.com   string, indicating method to dimulate the community:
#             "random" is a fully random lottery model;
#             "prop.neo" is weighted by the p.types vector;
#             Default is "random"
#
## Output : a list of elements with results of one set of simulated communities
# species.data    species.data,
# SR              Number of species in the community
# com = com       vector of community presence (1) and absence (0)
# simul.com       store the variable simul.com
# p.types         store the variable p.types
#_______________________________________________________________________

simul.comm <- function(species.pool,
                       SR = 25 ,
                       p.types = c(0.5, 0.25,0.25),
                       simul.com = c("random")) {
  
  
  if (class(species.pool) == "list")  species.pool <- species.pool$species.data 
  
  
  #### COMMUNITY ASSEMBLY
  
  # ASSEMBLY 0: totally random
  if (simul.com == "random") {
    com <- rep(0, nrow(species.pool))
    com[sample(nrow(species.pool), SR)] <- 1
    names(com) <- species.pool$species
  }
  
  # ASSEMBLY 2: proportion of neobiota is determined
  if (simul.com == "prop.neo") {
    com <- rep(0, nrow(species.pool))
    
    SRnat <- min (sum(species.pool$status == 1), floor(SR*p.types[1]), na.rm = T)
    SRarc <-  min (sum(species.pool$status == 2),floor(SR*p.types[2]), na.rm = T)
    SRneo <- min (sum(species.pool$status == 3), floor(SR*p.types[3]), na.rm = T)
    
    # native 
    if (SRnat>0)  com [sample(which(species.pool$status == 1), SRnat)] <- 1
    # archaeo
    if (SRarc>0)  com [sample(which(species.pool$status == 2), SRarc)] <- 1
    # neo
    if (SRneo>0) com [sample(which(species.pool$status == 3), SRneo)] <- 1
  }
  
  sim.com <-
    list(
      species.data = species.pool,
      SR = sum(com),
      com = com,
      simul.com = simul.com ,
      p.types = p.types
    )
  
  return(sim.com)
}

#### 2. FUNCTION simul.BNI:  Apply simulations and calculate BNI ####
# This funciton is a wrapper which applies functions simul.pool, simul.comm and BNI.calc.
# The funciton simulates a number of predefined scenarios for trait distributions,
# community assembly and species residence times, as presented in the manuscript
# by Schittko, Bernard-Verdier et al. (in review for GCB)
#
# INPUTS #
#
# scenario.traits   Trait scenario chosen among 6 options : 
#                      "fixed",
#                      "increasing mean", 
#                      "increasing variance",
#                      "lower to higher variance",
#                      "increasing mean and variance", 
#                      "increasing mean decreasing variance".
#                   Default is "fixed".
# scenario.date     Residence time (Date of arrival) scenario chosen among 4 options:
#                       "reference",
#                       "Archaeo range",
#                       "Neo recent range",
#                       "Neo threshold range".
#                   Default is "reference".
# simulation.com    Community assembly scenario chosen among 2 options: 
#                   "fixed" or "prop.neo". Default to "prop.neo".
# distance.method   Character string indicating the name of the distance method 
#                   to be used in BNI.calc. Must correspond to the methods in "dist" 
#                   of package "vegan". Default is "Euclidean".
# n.pools           Number of different species pool scenarios to be tested. Default is 20.
# nreps             Number of replicates for each species pool scenario. Default is 10.
# nb.species        Number of species to be simulated in the species pools. Default is 250.
# nb.sites          Number of community sites to be simulated. Default is 100.
# mean.richness     Species richness in each community. Default is 25.
# p.neo             Proportion of neophytes in each community
#                   (only valid if a fixed community scenario is used).
# proportion.status Numeric vector of length 3, giving the proportions of the three 
#                   different introduction categories of species. 
#                   Default is: c(nat = 0.70, arch = 0.15, neo = 0.15),
# time.arrival      Numeric vector of length 6, giving the min and max arrival time 
#                   of each species category. 
#                   Default is: c( nat.min = 8518,   nat.max = 8518,
#                                  arch.min = 2576,  arch.max = 2576,
#                                  neo.min = 1,      neo.max = 526)
# mean.tr           numeric vector of length 3, giving mean trait values
#                   for the three categories of species. Default is: c(0,0,0)
# sd.tr             numeric vector of length 3, giving SD trait values
#                   for the three categories of species. Default is: c(1,1,1)
#
# OUTPUT ##
## The function returns an array of simulation results of dimension nreps x npools.
## Each element of the list has the following structure:
#
# $index:   data.frame of n.site rows, with :
#             $BNI
#             $RaoQ
#             $BFI (Biotic index of Familiarity: BFI =  Rao's Q - BNI) 
#             $prop.neo: proportion or neobiota simulated
#             $hab : habitat type (1 to 5: defining % of neophyte expected)
#
# $species.data   data.frame of values for the simulated n.species : 
#                   $status:  (1:native, 2:archaeo, 3: neobiota),
#                   $YSI : years since introduction,
#                   $T1 to $T3 :  the three simulated traits 
#
# $simul.trait      scenario chosen for the trait simulation 
# $sd.trend         option chosen for the trait simulation when both mean and 
#                   SD of neobiota are different from natives 
#_______________________________________________________________________
simul.BNI <- function(scenario.traits = "fixed",
                      scenario.date = "reference",
                      simulation.com = "prop.neo",
                      distance.method =  "Euclidean",
                      n.pools = 20,
                      nreps = 10,
                      nb.species = 250,
                      nb.sites = 100,
                      mean.richness = 25,
                      proportion.status = c(nat = 0.70, arch = 0.15, neo = 0.15),
                      p.neo = 0.20,
                      p.arc = 0.10,
                      time.arrival = c(
                        nat.min = 8518,
                        nat.max = 8518,
                        arch.min = 2576,
                        arch.max = 2576,
                        neo.min = 1,
                        neo.max = 526
                      ),
                      mean.tr = c(0,0,0),
                      sd.tr = c(1,1,1)
                      ) {
#___________________________________________________________#    
# DATE scenarios across n.pools
#___________________________________________________________#    
    
    # Fixed dates across pools: (Scenario.date = "reference")
    time.mat <- matrix(time.arrival,
                       nrow = n.pools,
                       byrow = TRUE,
                       ncol = 6,
                       dimnames = list(
                         paste("pool",1: n.pools, sep ="_"),
                         names(time.arrival)
                       ) 
    )
    
    
    # Explore sensitivity to estimation of mean Archaeophyte arrivals
    if (scenario.date == "Archaeo range") {
      time.mat[,3] <- seq(time.arrival[3],
                          time.arrival[4], 
                          length.out = n.pools
      )
      time.mat[,4] <- time.mat[,3]
    }
    
    # explore sensitivity to the Archaeo/Neo threshold    
    if (scenario.date == "Neo threshold range") {
      time.mat[,5] <- time.arrival[5]
      time.mat[,6] <- seq(time.arrival[5],
                          time.arrival[6], 
                          length.out = n.pools
      )
    }
    
    # Explore sensitivity to recent arrivals
    if (scenario.date == "Neo recent range") {
      time.mat[,5] <- seq(time.arrival[5],
                          time.arrival[6], 
                          length.out = n.pools
      )
      time.mat[,6] <- time.arrival[6]
    }
    
    #___________________________________________________________#    
    # TRAIT scenarios across n.pools "scenario.traits" 
    #___________________________________________________________#    
    
    ## Start with no trait difference between introduction status:
    # => scenario.traits = "fixed"
    mean.tr.mat <- data.frame(matrix(mean.tr,
                                 nrow = n.pools,
                                 ncol= 3,
                                 byrow = TRUE,
                                 dimnames = list(
                                   paste("pool",1: n.pools, sep ="_"),
                                   c("I","A","N")
                                 )
    ))
    
    sd.tr.mat <- data.frame(matrix(sd.tr,
                               nrow = n.pools,
                               ncol= 3,
                               byrow = TRUE,
                               dimnames = list(
                                 paste("pool",1: n.pools, sep ="_"),
                                 c("I","A","N")
                               )
    ))
    
    # Scenario of increasing mean trait difference of neobiota 
    if (scenario.traits  == "increasing mean") {
      mean.tr.mat$N <- seq (0,10, length.out = n.pools)
    }
    # Scenario of increasing variance of neobiota 
    if (scenario.traits == "increasing variance") {
      sd.tr.mat$N <- seq (0,5, length.out = n.pools)
    }
    
    # Scenario of increasing variance of neobiota 
    if (scenario.traits == "lower to higher variance") {
      sd.tr.mat$N <- seq (0,5, length.out = n.pools)
      sd.tr.mat$I <- rep (2.5, n.pools)
      sd.tr.mat$A <- rep (2.5, n.pools)
    }
    # Scenario of increasing mean and variance of neobiota    
    if (scenario.traits == "increasing mean and variance") {
      mean.tr.mat$N <- seq (0,10, length.out = n.pools)
      sd.tr.mat$N <- seq (0,5, length.out = n.pools)
    }
    # Scenario of increasing mean and decreasing variance of neobiota  
    if (scenario.traits == "increasing mean decreasing variance") {
      mean.tr.mat$N <- seq (0,10, length.out = n.pools)
      sd.tr.mat$N <- seq (5,0, length.out = n.pools)
    }
    
    #___________________________________________________________#    
    # Proportion of neobiota in communities 
    #___________________________________________________________#    
    
    if (simulation.com =="prop.neo") {
      comm.types <- seq(0, 1, length.out = nb.sites)
    }
    if (simulation.com =="fixed") {
      comm.types <- rep(p.neo, nb.sites)
    }
   
    #___________________________________________________________#    
    # SIMULATIONS
    #___________________________________________________________#    
    
    ## simulate n.pools accordingt to the selected scenarios: 
    sim_list <- sapply(1:n.pools, function(i) {
      ## replicate nreps
     tmp <-  sapply(1:nreps, function(r){
        ## Build a random species pool
        spool <- simul.pool(nsp = nb.species,
                            p.status = proportion.status,
                            t.arrival = as.numeric(time.mat[i,]),
                            mean.traits = as.numeric(mean.tr.mat[i,]),
                            sd.traits = as.numeric(sd.tr.mat[i,]))
        
        ##  Build nsites communities from the pool :
        com.mat <- t(sapply(comm.types , function(k) {
          
          # Calculate the relative proportion of Archaeo and Natives:
          pA  <-  min(p.arc, 1-k)   # archaeobiota represent always p.arc = 10 % or less.
          ps <- c( (1-k-pA) , pA, k )
          
          # Simulate a community:
          scom <- simul.comm (spool$species.data,
                              SR = mean.richness ,
                              p.types = ps,
                              simul.com = simulation.com)
          
          return(scom$com)
        }, simplify = TRUE)
        )
      
        bni <- BNI.calc(com = com.mat,
                        trait.mat = spool$traits,
                        YSI = spool$species.data$YSI ,
                        dist.method = distance.method
                        )
        
        bni$index$prop.neo <- comm.types
        
        ## Add info on time of arrival settings
        bni$index <- cbind(bni$index, matrix(spool$t.arrival, byrow=TRUE,
                                             nrow = nrow(bni$index),
                                             ncol = 6,
                                             dimnames = list(
                                               NULL,
                                               c("nat.min" ,"nat.max",
                                                 "arch.min", "arch.max",
                                                 "neo.min", "neo.max")
                                               )
                                             )
        )
                           
        return(list(
            index = bni$index,
            com = com.mat,
            species.data = spool$species.data,
            scenario.traits = scenario.traits,
            simulation.com = simulation.com,
            scenario.date = scenario.date
            )
          )
      }, simplify = FALSE)
  }, simplify = TRUE)

    return(sim_list)
  }
#_______________________________________________________________________


#### 3. FUNCTION plot.simBNI :  Graphical output of BNI simulation ####
# Most appropriate for comparing Trait scenarios
#
# INPUTS #
# sim_list  an output object from simul.BNI
# type      vector of strings to select which combination 
#           of graphs to plot. Can be a combination of 
#           the following: c("Full", "BNI", "BNIs", "Rao").
#           Default to "Full".
# color.option     a character string indicating if the color gradient
#                  should vary with the different in trait mean or 
#                  SD of the neobiota in each simulation. 
#                  option = NULL colors will follow the mean variations.
# color.bar        logical: whether to plot a color bar legend
#
# which.bar        a character string indicating if two color 
#                  bars (SD and mean) should be represented together
#                  or in reverse order. 
#                  Default is NULL, only one bar will be plotted. 
#                 "inverse" is in reverse order. 
#                 "together" is in the same order.
# max.point       option to draw the maximum point of the LOESS curve,
#                 or the maximum of the simulated data. Default to "curve".
#
## Output : a four panel plot
#_______________________________________________________________________
plot.BNIsim <-
  function(sim_list,
           type = c("Full", "BNI", "BNIs", "Rao", "Traits"),
           color.option = c("mean"),
           max.point = "curve", 
           diff.col.max = TRUE,
           color.bar = TRUE, 
           which.bar = NULL) {
    
    require(fields)
    
    # Define color scheme  
    mean.seq = unlist(lapply(
      sim_list,
      FUN = function(x)
        mean(x$species.data$T1[x$species.data$status == 3])
    ))
    sd.seq = unlist(lapply(
      sim_list,
      FUN = function(x)
        sd(x$species.data$T1[x$species.data$status == 3])
    ))
    
    sd.nat.seq = unlist(lapply(
      sim_list,
      FUN = function(x)
        sd(x$species.data$T1[x$species.data$status !=3])
    ))
    
    
    allsd.seq = unlist(lapply(
      sim_list,
      FUN = function(x)
        sd(x$species.data$T1)
    ))
    
    
    col.seq = rep("#A1A1A1", length(sim_list))
    
    if (color.option == "mean") {
      sq <- floor((mean.seq - min(mean.seq)) / sd(mean.seq) *100) + 1
      col.seq = hcl.colors(max(sq))[sq]
    }
    
    if (color.option == "SD") {
      sq <- floor((sd.seq - min(sd.seq)) / sd(sd.seq) * 100) + 1
      col.seq = hcl.colors(max(sq))[sq]
    } 
    
    if (color.option == "allSD") {
      sq <- floor((allsd.seq - min(allsd.seq)) / sd(allsd.seq) * 100) + 1
      col.seq = hcl.colors(max(sq))[sq]
    } 
    
    col.seq = paste(col.seq, "50", sep = "")
    
      ## Plot only the trait values simulated
      if (type[1] == "Traits") {
        par(
          mfrow = c(1, 5),
          mar = c(1,0, 1, 0),
          oma = c(1, 4, 1, 1)
        )
        for (i in floor(seq(1,length(sim_list), length.out = 5)) ) {
          sim <- sim_list[[i]]
          boxplot(
            T1 ~ status,
            sim$species.data,
            names = levels(as.factor(sim$species.data$status)),
            ylim = c(-15, 15), 
            yaxt = "n"
          )
          if (i == 1) axis(2)
        }
      }
      
      ## BNI plots
      if (type[1] != "Traits") {
        
        # Select the plots to include
        nb.fig = length(type)
        if (type[1] == "Full")
          nb.fig = 4
        
        if (nb.fig == 1)
          par(
            mfrow = c(1, 1),
            mar = c(4, 4, 1, 2),
            oma = c(1, 1, 1, 1)
          )
        if (nb.fig == 2)
          par(
            mfrow = c(1, 2),
            mar = c(4, 4, 1, 2),
            oma = c(1, 1, 1, 1)
          )
        if (nb.fig == 3)
          par(
            mfrow = c(1, 3),
            mar = c(4, 4, 1, 2),
            oma = c(1, 1, 1, 1)
          )
        if (nb.fig == 4)
          par(
            mfrow = c(2, 2),
            mar = c(4, 4, 1, 2),
            oma = c(1, 1, 1, 1)
          )
        
        
        ## Illustrate max simulated trait differences
        if ("Full" %in%  type) {
          sim_max <- sim_list[[length(sim_list)]]
          boxplot(
            T1 ~ status,
            sim_max$species.data,
            names = c("Nat", "Arch", "Neo"),
            ylab = "Simulated trait values"
          )
        }
        
        ## BNI
        if ("BNI" %in% type | "Full" %in% type) {
          nd = data.frame(prop.neo = seq(0, 1, 0.01))
          # y.max <-
          #   max(unlist(lapply(
          #     sim_list,
          #     FUN = function(x)
          #       max(x$index$BNI)
          #   )))
          y.max <- 4
          
          max(c(sim_list[[1]]$index$BNI, sim_list[[20]]$index$BNI))
          plot(c(0, 1), c(0,  y.max),
               type = "n",  ann = F)
          title(xlab = "Proportion of neobiota", ylab = "BNI")
          
          
          for (i in 1:length(sim_list)) {
            lo <- loess(BNI ~ prop.neo, data = sim_list[[i]]$index, span = 0.9)
            plo <-predict(lo,  se = TRUE)
            lines(sim_list[[i]]$index$prop.neo, plo$fit, col = col.seq[i])
            col.max = "black"
            if (diff.col.max) {
              col.max = "firebrick"
              if (sd.seq[i] > sd.nat.seq[i]) col.max <- "black"
            }
            
            
            if(max.point =="curve"){
              mi <- which(plo$fit== max(plo$fit))
              points(seq(0, 1, 0.01)[mi],
                     plo$fit[mi],
                     col = col.max, cex = 0.8, pch = 20 )
            }
            
            if(max.point =="data"){
              mi <- which(sim_list[[i]]$index$BNI == max(sim_list[[i]]$index$BNI))
              points(sim_list[[i]]$index[mi,"prop.neo"],
                     sim_list[[i]]$index[mi, "BNI"], 
                     col = col.max, cex = 0.8, pch = 20 )
            }
          }
        }
        
        ## Rao
        if ("Rao" %in% type | "Full" %in% type) {
          nd = data.frame(prop.neo = seq(0, 1, 0.01))
          
          # y.max <-
          #   max(unlist(lapply(
          #     sim_list,
          #     FUN = function(x)
          #       max(x$index$RaoQ)
          #   )))
          y.max <- 4
          
          
          plot(c(0, 1), c(0, y.max),
               type = "n",  ann = F)
          title(xlab = "Proportion of neobiota", ylab = "Rao's Q")
          
          for (i in 1:length(sim_list)) {
            lo <- loess(RaoQ ~ prop.neo, data = sim_list[[i]]$index, span = 0.9)
            plo <- predict(lo, se = TRUE)
            lines(sim_list[[i]]$index$prop.neo, plo$fit, col = col.seq[i])
            
            col.max = "black"
            if (diff.col.max) {
              col.max = "firebrick"
              if (sd.seq[i] > sd.nat.seq[i]) col.max <- "black"
            }
            
            if(max.point =="curve"){
              mi <- which(plo$fit== max(plo$fit))
              points(seq(0, 1, 0.01)[mi],
                     plo$fit[mi],
                     col = col.max, cex = 0.8, pch = 20 )
            }
            
          }
        }
        
        
        
        
        ## BNI/Rao
        if ("BNIs" %in% type | "Full" %in% type) {
          nd = data.frame(prop.neo = seq(0, 1, 0.01))
          plot(c(0, 1), c(0, 1),
               type = "n",  ann = F)
          title(xlab = "Proportion of neobiota", ylab = "BNIs")
          
          for (i in 1:length(sim_list)) {
            BNIs <- (sim_list[[i]]$index$BNI / sim_list[[i]]$index$RaoQ)
            lo <- loess(BNIs ~ prop.neo, data = sim_list[[i]]$index, span = 0.6)
            plo <- predict(lo, se = TRUE)
            lines(sim_list[[i]]$index$prop.neo, plo$fit, col = col.seq[i])
          }
        }
        
        ## Draw legend bar
        
        if (color.bar) {
          if (color.option == "mean" & is.null(which.bar)) {
            sd.seq = unlist(lapply(
              sim_list,
              FUN = function(x)
                sd(x$species.data$T1[x$species.data$status == 3])
            ))
            image.plot(
              legend.only = TRUE,
              zlim = c(0, max(mean.seq)),
              col = hcl.colors(200),
              smallplot = c(0.70, 0.73, 0.4, 0.6),
              legend.lab = "neobiota mean",
              legend.line = 2,
              legend.cex = 0.7,
              border = "white",
              axis.args = list(cex.axis = 0.7)
            )
          }
          
          if (color.option == "SD" &  is.null(which.bar)) {
            image.plot(
              legend.only = TRUE,
              zlim = c(0, 5),
              col = hcl.colors(200),
              smallplot = c(0.70, 0.73, 0.4, 0.6),
              border = "white",
              axis.args = list(cex.axis = 0.7),
              legend.args=list( text = "neobiota SD", side=4, font=2, line=2, cex=0.7)
            )
          }
          
          if (!is.null(which.bar)) {
            my.image.plot(
              legend.only = TRUE,
              zlim = c(0, 10),
              col = hcl.colors(200),
              smallplot = c(0.70, 0.73, 0.4, 0.6),
              border = "white",
              axis.args = list(cex.axis = 0.7), axis.side = 4,
              legend.args=list( text = "neobiota mean", side=4, font=2, line=2, cex=0.7)
            )
            
            my.image.plot(
              legend.only = TRUE,
              col = hcl.colors(200),
              zlim = c(0, 10),
              smallplot = c(0.70, 0.73, 0.4, 0.6),
              border = "white",
              axis.args = list( cex.axis = 0.7),
              axis.lab = (if (which.bar == "inverse") 5:0 else 0:5),
              axis.at = seq(0, 10,2),
              axis.side = 2,
              legend.args=list( text = "neobiota SD", side=2, font=2, line=2, cex=0.7)
            )
          }
          
        }
      }
  }




#_______________________________________________________________________


#### 3. FUNCTION plot.dates.BNI :  Graphical output of BNI sensitivity to time ####

## Inputs ##
# sim_list  an output object from simul.BNI
#
## Output : a 2 panel plot
#_______________________________________________________________________
plot.sim.date.BNI <-function(sim_list ) {
  
  par(
    mfrow = c(1, 3),
    mar = c(4, 4, 1, 2),
    oma = c(1, 1, 1, 1)
  )
  

  if (sim_list[[1]]$scenario.date == "Archaeo range") {
    xaxs.lab = "Mean residence time of Archaeobiota"
    col.nam = "arch.min"
    legend.label = "Proportion of neobiota"
  } 
  
  if (sim_list[[1]]$scenario.date  == "Neo threshold range") {
    xaxs.lab = "Earliest non-native"
    col.nam = "neo.max"
    legend.label = "Proportion of non-native"
       } 
  
  if (sim_list[[1]]$scenario.date  == "Neo recent range") {
    xaxs.lab = "Most recent non-native"
    col.nam = "neo.min"
    legend.label = "Proportion of non-native"
  } 
  
 indices.mat <- do.call(
   args = sapply(1:length( sim_list),
    FUN = function(x) sim_list[[x]]$index,
    simplify = FALSE),
    what = rbind
   )
 
 indices.mat <- indices.mat[order(indices.mat$prop.neo),]

 indices.mat$perc.neo = 1 + ceiling(indices.mat$prop.neo  *100)
 col.seq = hcl.colors(100,palette ="Tealrose" )
 col.seq = paste( col.seq, "40", sep = "")
 col.seq = col.seq[indices.mat$perc.neo]
 
 ## BNI
  plot(range(indices.mat[,col.nam]), c(0,  max(indices.mat$BNI)),
       type = "n",  ann = F)
  title(xlab = xaxs.lab, ylab = "BNI")
  points(as.formula(paste('BNI ~ ', col.nam)),
         data = indices.mat,
         pch = 20,
         col = col.seq)
  
  ## BNI/Rao
  
  plot(range(indices.mat[,col.nam]), c(0,  max(indices.mat$BNIs)),
       type = "n",  ann = F)
  title(xlab = xaxs.lab, ylab = "BNIs")
  points(as.formula(paste('BNIs ~ ', col.nam)),
         data = indices.mat,
         pch = 20,
         col = col.seq)
  
  
  
plot.new()
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      col = hcl.colors(100,palette ="Tealrose" ),
      smallplot = c(0.20, 0.23, 0.3, 0.8),
      legend.lab = legend.label,
      legend.line = 2,
      legend.cex = 0.7,
      border = "white",
      axis.args = list(cex.axis = 0.7)
    )
  
}


