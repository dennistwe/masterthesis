########################################################

# Functions to fit degenarative 1:2 binding models and
#                  cooperative 1:2 binding models

# Required packages: minpack.lm
########################################################
var_nls_deg <- function(..., start.conc){
  nls.model <- nlsLM(...)
  
  g.current <- start.conc
  k1.current <- as.numeric(coef(nls.model)["k1"])
  # max.iter <- 0
  sigma.new <- 2
  sigma.curr <- 1
  sigma.it <- c(sigma(nls.model) + 1, sigma(nls.model))
  results <- list()
  results$sigma.model[[1]] <- sigma(nls.model)   
  results$k1[[1]] <- k1.current
  results$trna.conc[[1]] <- g.current
  
  while(sigma.it[sigma.curr] >= sigma.it[sigma.new]){
    if(length(coef(nls.model)) == 2){ 
      model.loop <- nlsLM(que_fi ~ bm1to2.deg.add(fhg, k1, g.current),
                          data = fi,
                          start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                       k1 = k1.current))
    } else {
      model.loop <- nlsLM(que_fi ~ bm1to2.deg(fhg, fhg2, k1, g.current),
                          data = fi,
                          start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                       fhg2 = as.numeric(coef(nls.model)["fhg2"]),
                                       k1 = k1.current))
    }
    
    sigma.curr <- sigma.curr +1
    sigma.new <- sigma.new +1
    sigma.it <- c(sigma.it, sigma(model.loop))
    g.current <- trna_eq(constants = c(coef(model.loop)["k1"], 
                                       coef(model.loop)["k1"]/4, 
                                       5),
                         start.conc)
    k1.current <- as.numeric(coef(model.loop)["k1"])
    print(sigma(model.loop))
    
    results$sigma.model[[sigma.curr]] <- sigma(model.loop)
    results$k1[[sigma.curr]] <- k1.current
    results$trna.conc[[sigma.curr]] <- g.current
    
  }
  
  #selector <- which(results$sigma.model == min(results$sigma.model)) - 1 
  selector <- which(unlist(results$sigma.model) == min(unlist(results$sigma.model))) -1
  g.best <- unlist(results$trna.conc[selector])
  k1.best <- unlist(results$k1[selector])
  ###
  results.list <- list()
  if(length(coef(nls.model)) == 2){
    nls.model.best <- nlsLM(
      que_fi ~ bm1to2.deg.add(fhg, k1, g0 = g.best),
      data = fi,
      start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                   k1 = k1.best))
    
    results.list <- list(bm1to2.deg.add = nls.model.best,
                         g = g.best) 
    
  } else {
    nls.model.best <- nlsLM(que_fi ~ bm1to2.deg(fhg, fhg2, k1, g.best),
                            data = fi,
                            start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                         fhg2 = as.numeric(coef(nls.model)["fhg2"]),
                                         k1 = k1.best))
    results.list <- list(bm1to2.deg = nls.model.best,
                         g = g.best) 
    
  }
  
}


var_nls_coop <- function(..., start.conc){
  nls.model <- nlsLM(...)
  
  g.current <- start.conc
  k1.current <- as.numeric(coef(nls.model)["k1"])
  k2.current <- as.numeric(coef(nls.model)["k2"])
  # max.iter <- 0
  sigma.new <- 2
  sigma.curr <- 1
  sigma.it <- c(sigma(nls.model) + 1, sigma(nls.model))
  results <- list()
  results$sigma.model[[1]] <- sigma(nls.model)   
  results$k1[[1]] <- k1.current
  results$k2[[1]] <- k2.current
  results$trna.conc[[1]] <- g.current
  
  while(sigma.it[sigma.curr] >= sigma.it[sigma.new]){
    if(length(coef(nls.model)) == 3){ 
      model.loop <- nlsLM(que_fi ~ bm1to2.coop.add(fhg, k1, k2, g.current),
                          data = fi,
                          start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                       k1 = k1.current,
                                       k2 = k2.current))
    } else {
      model.loop <- nlsLM(que_fi ~ bm1to2.coop(fhg, fhg2, k1, k2, g.current),
                          data = fi,
                          start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                       fhg2 = as.numeric(coef(nls.model)["fhg2"]),
                                       k1 = k1.current,
                                       k2 = k2.current))
    }
    sigma.curr <- sigma.curr +1
    sigma.new <- sigma.new +1
    sigma.it <- c(sigma.it, sigma(model.loop))
    g.current <- trna_eq(constants = c(coef(model.loop)["k1"], 
                                       coef(model.loop)["k2"], 
                                       5),
                         start.conc)
    k1.current <- as.numeric(coef(model.loop)["k1"])
    k2.current <- as.numeric(coef(model.loop)["k2"])
    
    print(sigma(model.loop))
    results$sigma.model[[sigma.curr]] <- sigma(model.loop)
    results$k1[[sigma.curr]] <- k1.current
    results$k2[[sigma.curr]] <- k2.current
    results$trna.conc[[sigma.curr]] <- g.current
    
  }
  
  #selector <- which(results$sigma.model == min(results$sigma.model)) - 1 
  selector <- which(unlist(results$sigma.model) == min(unlist(results$sigma.model))) -1
  g.best <- unlist(results$trna.conc[selector])
  k1.best <- unlist(results$k1[selector])
  k2.best <- unlist(results$k2[selector])
  results.list <- list()
  
  if(length(coef(nls.model)) == 3){
    nls.model.best <- nlsLM(que_fi ~ bm1to2.coop.add(fhg, k1, k2, g.best),
                            data = fi,
                            start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                         k1 = k1.best,
                                         k2 = k2.best))
    results.list <- list(bm1to2.coop.add = nls.model.best,
                         g = g.best) 
  } else {
    nls.model.best <- nlsLM(que_fi ~ bm1to2.coop(fhg, fhg2, k1, k2, g.best),
                            data = fi,
                            start = list(fhg = as.numeric(coef(nls.model)["fhg"]),
                                         fhg2 = as.numeric(coef(nls.model)["fhg2"]),
                                         k1 = k1.best,
                                         k2 = k2.best))
    results.list <- list(bm1to2.coop = nls.model.best,
                         g = g.best) 
  }
}