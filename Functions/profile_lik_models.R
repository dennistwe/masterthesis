############################################################

# likelihood and profile likelihood functions

############################################################

ll.1to2.deg.add <- function(param, obs, g.curr, profile.par = NULL, ...){

  g0 <- g.curr
  if (is.null(profile.par)){
    Fpred <- bm1to2.deg.add(param["fhg"], param["k1"], g0)
  }  else if (all(profile.par == "fhg")){
    Fpred <- bm1to2.deg.add(..., param["k1"], g0)
  } else if (all(profile.par == "k1")){
    Fpred <- bm1to2.deg.add(param["fhg"], ..., g0)
  } else if (all(profile.par == c("fhg", "k1"))){
    Fpred <- bm1to2.deg.add(..., g0)
  }
  
  lnlik <- sum(dnorm(x = obs$que_fi, 
                     mean = Fpred, 
                     sd = param["sig"], 
                     log = TRUE))
  return(lnlik)
}
###########################
ll.1to2.deg <- function(param, obs, g.curr, profile.par = NULL, ...){
  g0 <- g.curr
  if (is.null(profile.par)){
    Fpred <- bm1to2.deg(param["fhg"], param["fhg2"], param["k1"], g0)
  } else if (all(profile.par == "fhg")){
    Fpred <- bm1to2.deg(..., param["fhg2"], param["k1"], g0)
  } else if (all(profile.par == "fhg2")){
    Fpred <- bm1to2.deg(param["fhg"], ..., param["k1"], g0)
  } else if (all(profile.par == "k1")){
    Fpred <- bm1to2.deg(param["fhg"], param["fhg2"], ..., g0)
  } else if (all(profile.par == c("fhg", "fhg2"))){
    Fpred <- bm1to2.deg(..., param["k1"], g0)
  } else if (all(profile.par == c("fhg", "k1"))){
    Fpred <- bm1to2.deg(..., param["fhg2"], g0)
  } else if (all(profile.par == c("fhg2", "k1"))){
    Fpred <- bm1to2.deg(param["fhg"], ..., g0)
  }
  
  lnlik <- sum(dnorm(x = obs$que_fi, 
                     mean = Fpred, 
                     sd = param["sig"], 
                     log = TRUE
                     )
               )
  return(lnlik)
}

###############################################################
ll.1to2.coop.add <- function(param, obs, g.curr, profile.par = NULL, ...){
  g0 <- g.curr
  if (is.null(profile.par)){
    Fpred <- bm1to2.coop.add(param["fhg"], param["k1"], param["k2"], g0)
  } else if (all(profile.par == "fhg")){
    Fpred <- bm1to2.coop.add(..., param["k1"], param["k2"], g0)
  } else if (all(profile.par == "k1")){
    Fpred <- bm1to2.coop.add(param["fhg"], ..., param["k2"], g0)
  } else if (all(profile.par == "k2")){
    Fpred <- bm1to2.coop.add(param["fhg"], param["k1"], ..., g0)
  } else if (all(profile.par == c("fhg", "k1"))){
    Fpred <- bm1to2.coop.add(..., param["k2"], g0)
  } else if (all(profile.par == c("fhg", "k2"))){
    Fpred <- bm1to2.coop.add(..., param["k1"], g0)
  } else if (all(profile.par == c("k1", "k2"))){
    Fpred <- bm1to2.coop.add(param["fhg"], ..., g0)
  }
  
  lnlik <- sum(dnorm(x = obs$que_fi, 
                     mean = Fpred, 
                     sd = param["sig"], 
                     log = TRUE))
  return(lnlik)
}

#######################################################################
ll.1to2.coop <- function(param, obs, g.curr, profile.par = NULL, ...){
  g0 <- g.curr
  if (is.null(profile.par)){
    Fpred <- bm1to2.coop(param["fhg"], param["fhg2"], param["k1"], 
                         param["k2"], g0)
  } else if (all(profile.par == "fhg")){
    Fpred <- bm1to2.coop(..., param["fhg2"], param["k1"], param["k2"], g0)
  } else if (all(profile.par == "fhg2")){
    Fpred <- bm1to2.coop(param["fhg"], ..., param["k1"], param["k2"], g0)
  } else if (all(profile.par == "k1")){
    Fpred <- bm1to2.coop(param["fhg"], param["fhg2"], ..., param["k2"], g0)
  } else if (all(profile.par == "k2")){
    Fpred <- bm1to2.coop(param["fhg"], param["fhg2"], param["k1"], ..., g0)
  } else if (all(profile.par == c("fhg", "fhg2"))){
    Fpred <- bm1to2.coop(..., param["k1"], param["k2"], g0)
  } else if (all(profile.par == c("fhg", "k1"))){
    Fpred <- bm1to2.coop(..., param["fhg2"], param["k2"], g0)
  } else if (all(profile.par == c("fhg", "k2"))){
    Fpred <- bm1to2.coop(..., param["fhg2"], param["k1"], g0)
  } else if (all(profile.par == c("fhg2", "k1"))){
    Fpred <- bm1to2.coop(param["fhg"], ..., param["k2"], g0)
  } else if (all(profile.par == c("fhg2", "k2"))){
    Fpred <- bm1to2.coop(param["fhg"], ..., param["k1"], g0)
  } else if (all(profile.par == c("k1", "k2"))){
    Fpred <- bm1to2.coop(param["fhg"], param["fhg2"], ..., g0)
  }
  
  lnlik <- sum(dnorm(x = obs$que_fi, 
                     mean = Fpred, 
                     sd = param["sig"], 
                     log = TRUE))
  return(lnlik)
}