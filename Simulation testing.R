# Simulations testing the range of BNI, BNIs and Rao's Q values
# This script will reproduce the main figures in the Supplementary material 2
# of the manuscript by Schittko, Bernard-Verdier et al., 2020, GCB.

# Load the necessary functions:
source("FUNCTION BNI simulation functions.R")
source("BNI function.R")

# Load packages
library(vegan)
library(FD)

#___________________________________________________________#
#  Tested TRAIT scenarios                                ####
#___________________________________________________________#

# TRAIT SCENARIO 1: no trait difference between natives and Neobiota ####
simFixed <- simul.BNI(scenario.traits = "fixed",
                      n.pools = 20,
                      nreps = 10,
                       proportion.status = c(nat = 0.70, arch = 0.15, neo = 0.15))

plot.BNIsim(simFixed, type = c("Full"), color.option = "",color.bar = FALSE, diff.col.max = TRUE )

# TRAIT SCENARIO 2 : increasing trait differences between natives and aliens
simTraitDif <- simul.BNI( scenario.traits = "increasing mean",
                          n.pools = 20,
                          nreps = 10)

plot.BNIsim(simTraitDif, type = "Full",  diff.col.max = FALSE )


# TRAIT SCENARIO 3a : increasing trait variance of neobiota
simSDDif <- simul.BNI(scenario.traits = "increasing variance",
                      n.pools = 20,
                      nreps = 10)
plot.BNIsim(simSDDif , type = "Full", color.option = "SD")

# TRAIT SCENARIO 3b : increasing trait variance of neobiota
#                     when natives already have a high variance
simSDDif.natSD <- simul.BNI(scenario.traits = "lower to higher variance",
                            n.pools = 20,
                            nreps = 10)

plot.BNIsim(simSDDif.natSD , type = "Full", color.option = "SD")

# TRAIT SCENARIO 4a : increasing trait mean &s variance of neobiota
simbothDif <- simul.BNI(scenario.traits = "increasing mean and variance",
                        n.pools = 20,
                        nreps = 10)
plot.BNIsim(simbothDif , type = "Full",  color.option ="SD", which.bar = "together")


# TRAIT SCENARIO 4a :  increasing trait mean & decreasing variance of neobiota
simbothDif.inverse <- simul.BNI(scenario.traits = "increasing mean decreasing variance",  
                                n.pools = 20,
                                nreps = 10)
plot.BNIsim(simbothDif.inverse , type = "Full",   color.option ="mean", which.bar = "inverse")




#___________________________________________________________#
#  Sensitivity to DATE ESTIMATION scenarios              ####
#___________________________________________________________#

# DATE SCENARIO 2: Sensitivity to the estimation of mean Archaeobiota residence times: 
simArchDates <- simul.BNI(scenario.date = "Archaeo range",
                                 simulation.com = "prop.neo",
                                 time.arrival = c(nat.min= 8518,nat.max= 8518,
                                                  arch.min = 2000, arch.max = 4000,
                                                  neo.min = 100, neo.max = 100),
                                p.arc = 0.15,
                                n.pools = 50 , nreps = 10 )
 
 plot.sim.date.BNI(simArchDates)

 
# DATE SCENARIO 3: Effect of the most recent neobiota arrival (from 1 to 5000 years ago)
 # In scenario 3 and 4, we remove effect of acheobiota, and only test effect of
 # the max and min dates of introduction for non-natives.
 # We also give neobiota a high trait variance for more striking results.
 
 simRecent.Neo.only <- simul.BNI(scenario.traits = "fixed",
                            scenario.date = "Neo recent range",
                            simulation.com = "prop.neo",
                            # We give neos a high trait variance for more striking results:
                            sd.tr = c(1,1,2),   
                            time.arrival = c(nat.min= 8518,nat.max= 8518,
                                             # archaeo are reduced to natives:
                                             arch.min = 8518, arch.max = 8518,
                                             neo.min = 1, neo.max = 5000),
                            n.pools = 50 , nreps = 10)

 plot.sim.date.BNI( simRecent.Neo.only )
 
 
# DATE SCENARIO 4: Effect of the most ancient neobiota arrival (from 5000 to 1 years ago)
 simOlder.Neo.only <- simul.BNI(scenario.traits = "fixed",
                                 scenario.date = "Neo threshold range",
                                 simulation.com = "prop.neo",
                                  sd.tr = c(1,1,2),
                                 time.arrival = c(nat.min= 8518,
                                                  nat.max= 8518,
                                                  arch.min = 8518,
                                                  arch.max = 8518,
                                                  neo.min = 1,
                                                  neo.max = 5000),
                                 n.pools = 50 , nreps = 10)

 plot.sim.date.BNI( simOlder.Neo.only )
 
 
