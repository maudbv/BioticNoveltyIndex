# BIOTIC NOVELTY INDEX CALCULATION
#
# This R script is a function calculating the BNI, an index of biotic novelty 
# for biological communities. The function also calculates the 
# standardised version of the index (BNIs) and Rao's quadratic entropy.
#
# Requires R packages: FD, proxy
#
# written by Conrad Schittko and Maud Bernard-Verdier.
#
# Updated February 2020
#________________________________________________________________#
#
#
## INPUT ##
# com           numeric community matrix of species abundance 
#               of the format (sites x species)
# trait         numeric matrix of trait values which will be used
#               to calculate multitrait distances (species x traits)
# YSI           numeric vector of "years since introduction"
#               for each species in the region of interest (in years).
#               Must be in the same order as species in com and trait.
# dist.method   character string, name of a distance measure
#               accepted by the function "dist" in the R package vegan,
#               or "gower" which calculated the gower distance with package FD.
#               Default is "gower".
# scale.dist    logical value indicated if the distance matrix should be
#               scaled between 0 and 1. Default is FALSE.
# dist.mat      (optional) matrix of distance for all the species in "com" and "trait.mat"
#               ( !beware : Species must be in the same order as in com and YSI).
#               If provided, this will bypass the "trait.mat" and "dist.method".
#               It can be used to input phylogenetic distances.
#
#
## OUTPUT ##  a list of the following elements:
# index:      a data frame with the caluclated BNI, BNIs
#             and Rao's quadratic entropy for each site
# trait.mat   matrix of trait values provided in the input (if any)
# YSI         vector of year since introduction provided in the input
# rpi         vector of scaled residence times for each species
# cij         matrix of temporal coexistence coefficient for each pair of species
# t.dist      matrix of distances between pairs of species

#________________________________________________________________#


BNI.calc <- function(com,
                     YSI,
                     trait.mat = NULL,
                     dist.mat = NULL,
                     dist.method = "gower",
                     scale.dist = FALSE) {
  
  
  # Error if no distance matrix and no trait matrix are provided
  try(if(is.null(dist.mat) & is.null(trait.mat)) stop("No trait or distance matrix provided"))
  
  require(FD)
  require(proxy)
  
  # If only one community vector instead of a community matrix
  if (is.vector(com))
    com <- t(as.matrix(com))
  
  # Transform community matrix in relative abundances
  com <- com / rowSums(com)
   
  # If no distance matrix is provided:
  if (is.null(dist.mat)) {
  # Calculate pairwise trait distances "d(ij)"
  if (dist.method == "gower") {
    t.dist <- FD::gowdis(x = trait.mat) # Package FD
    t.dist <- as.matrix(t.dist)
  } else {
    t.dist <- dist(trait.mat, method = dist.method)
    t.dist <- as.matrix(t.dist)
  }
  }
  
  # If a distance matrix is already provided, no need to caluclate it
  # This may be used to calculate a phylogenetic version of the indices.
  if (!is.null(dist.mat)) t.dist = dist.mat
  
  
  # scale distances
  if (scale.dist)
    t.dist <- t.dist / max(t.dist, na.rm = TRUE)
  
  # calculate the scaled time of residence of each species
  rpi <-  (YSI) / (max(YSI))
  
  # calculate pairwise temporal coefficients of coexistence
  cij <- as.matrix(1 - dist(rpi, method = min))
  
  # combine the matrix of trait distance and the temporal coefficient
  # into one weight matrix : "d(ij) x c(ij)"
  dxc <- t.dist * cij
  dxc <- as.matrix(dxc)
  
  # Calculate the BNI as a cross product of two matrices:
  # community matrix x weight matrix
  BNI <- apply(com, 1, function(x)
    crossprod(x, dxc %*% x)) / 2
  
  # Calculate also Rao's quadratic entropy
  RaoQ <- apply(com, 1, function(x)
    crossprod(x, t.dist %*% x)) / 2
  
  # Calculate the standardized index BNIs
  BNIs <- BNI / RaoQ
  
  return(
    list(
      index = data.frame(BNI = BNI,
                         RaoQ = RaoQ,
                         BNIs = BNIs),
      trait.mat = trait.mat,
      YSI = YSI,
      rpi = rpi,
      cij = cij,
      t.dist = t.dist
    )
  )
}
