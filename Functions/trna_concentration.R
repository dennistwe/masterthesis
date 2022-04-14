#################################################################

# Functions to get equilibrium trna concentrations

#################################################################

# functions work in base R

#root.finder:
# root.finder solves the cubic polynomial and finds the roots for
# each starting concentration g0 of trna

# trna_eq:
# trna_eq takes root.finder as nested function and extracts the 
# smallest positive real roots as solutions for every starting 
# concentration g0

###################################################################
root.finder <- 
  function(constants = c(k1 = 0, k2 = 0, h0 = 0), g0){
  A <- constants[1] * constants[2]
  B <- 2*constants[1] * constants[2] * constants[3] - 
    constants[1] * constants[2] * g0 + constants[1]
  C <- constants[1] * constants[3] - constants[1] * g0 + constants[1]
  D <- -g0
  cube.solv <- polyroot(c(D,C,B,A))
  return(cube.solv)
}

trna_eq <- function(constants = c(k1 = 0, k2 = 0, h0 = 0), g0){
  b1 <- sapply(g0, function(x) root.finder(constants, x))
  find_zero_imaginary <- sapply(b1, function(x) round(Im(x),digits = 6) == 0)
  b2 <- b1[find_zero_imaginary]
  find_positive_real <- sapply(b2, function(x) Re(x) > 0)
  b3 <- b2[find_positive_real]
  b3 <- Re(b3)
 return(b3)
 
}

