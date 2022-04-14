#Model1
bm1to2.coop <- 
  function(fhg, fhg2, k1, k2, g0){
    5 * (fhg * k1 * g0 + fhg2 * k1 * k2 * g0^2)/
      (1 + k1 * g0 + k1 * k2 * g0^2)
  }

#Model2
bm1to2.coop.add <- 
  function(fhg, k1, k2, g0){
    5 * (fhg * (k1 * g0 + 2 * k1 * k2 *g0^2))/
      (1 + k1 * g0 + k1 * k2 * g0^2)
  }

#Model3
bm1to2.deg <- 
  function(fhg, fhg2, k1, g0){
    5 * (fhg * k1 * g0 + fhg2 * k1 *(k1 / 4) * g0^2)/
      (1 + k1 * g0 + k1 * (k1 / 4) * g0^2)
  }

#Model4
bm1to2.deg.add <- 
  function(fhg, k1, g0){
    5 * (fhg * (k1 * g0 + 2 * k1 * (k1 / 4) * g0^2))/
      (1 + k1 * g0 + k1 * (k1 / 4) * g0^2)
  }

#Model5
bm1to1 <- 
  function(fhg,k1,g0){
    fhg * (((g0 + 5 + (1 / k1)) - 
              sqrt((g0 + 5 + (1 / k1))^2 - 4 * 5 * g0)) / 2)
  }
