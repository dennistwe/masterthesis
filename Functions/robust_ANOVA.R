lincon<-function(x,con=0,tr=.2,alpha=.05,pr=FALSE){
  #
  #  A heteroscedastic test of d linear contrasts using trimmed means.
  #
  #  This version uses an improved method for computing the quantiles of a
  #  Studentized maximum modulus distriburtion
  # 
  #  The data are assumed to be stored in $x$ in list mode, a matrix
  #  or a data frame. If in list mode,
  #  length(x) is assumed to correspond to the total number of groups.
  #  It is assumed all groups are independent.
  #
  #  con is a J by d matrix containing the contrast coefficients that are used.
  #  If con is not specified, all pairwise comparisons are made.
  #
  #  Missing values are automatically removed.
  #
  #  pr=FALSE included to avoid errors using an earlier version of this function when
  #   dealing with two-way and higher designs
  #
  #  Adjusted p-values are based on the Studentized maximum modulus distribution with the 
  #  goal of controlling FWE
  #
  #  To apply the Kaiser-Bowden method, use the function kbcon
  #
  if(tr==.5) stop('Use the R function medpb to compare medians')
  if(is.data.frame(x)) x = as.matrix(x)
  if(is.matrix(x)) x <- listm(x)
  if(!is.list(x))stop('Data must be stored in a matrix or in list mode.')
  con<-as.matrix(con)
  J<-length(x)
  sam=NA
  h<-vector('numeric',J)
  w<-vector('numeric',J)
  xbar<-vector('numeric',J)
  for(j in 1:J){
    xx<-!is.na(x[[j]])
    val<-x[[j]]
    x[[j]]<-val[xx]  # Remove missing values
    sam[j]=length(x[[j]])
    h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
    # h is the number of observations in the jth group after trimming.
    w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
    xbar[j]<-mean(x[[j]],tr)
  }
  if(sum(con^2)==0){
    CC<-(J^2-J)/2
    psihat<-matrix(0,CC,9)
    dimnames(psihat)<-list(NULL,c('Group','Group','psihat','ci.lower','ci.upper',
                                  'p.value','Est.1','Est.2','adj.p.value'))
    test<-matrix(NA,CC,6)
    dimnames(test)<-list(NULL,c('Group','Group','test','crit','se','df'))
    jcom<-0
    for (j in 1:J){
      for (k in 1:J){
        if (j < k){
          jcom<-jcom+1
          test[jcom,3]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
          sejk<-sqrt(w[j]+w[k])
          test[jcom,5]<-sejk
          psihat[jcom,1]<-j
          psihat[jcom,2]<-k
          test[jcom,1]<-j
          test[jcom,2]<-k
          psihat[jcom,3]<-(xbar[j]-xbar[k])
          df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
          test[jcom,6]<-df
          psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
          psihat[jcom,7]=xbar[j]
          psihat[jcom,8]=xbar[k]
          crit=qsmm(1-alpha,CC,df)
          test[jcom,4]<-crit
          psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
          psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
          psihat[jcom,9]=1-psmm(test[jcom,3],CC,df)
        }}}}
  if(sum(con^2)>0){
    if(nrow(con)!=length(x)){
      stop('The number of groups does not match the number of contrast coefficients.')
    }
    CC=ncol(con)
    psihat<-matrix(0,ncol(con),6)
    dimnames(psihat)<-list(NULL,c('con.num','psihat','ci.lower','ci.upper',
                                  'p.value','adj.p.value'))
    test<-matrix(0,ncol(con),5)
    dimnames(test)<-list(NULL,c('con.num','test','crit','se','df'))
    df<-0
    for (d in 1:ncol(con)){
      psihat[d,1]<-d
      psihat[d,2]<-sum(con[,d]*xbar)
      sejk<-sqrt(sum(con[,d]^2*w))
      test[d,1]<-d
      test[d,2]<-sum(con[,d]*xbar)/sejk
      df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
      crit=qsmm(1-alpha,CC,df)
      test[d,3]<-crit
      test[d,4]<-sejk
      test[d,5]<-df
      psihat[d,3]<-psihat[d,2]-crit*sejk
      psihat[d,4]<-psihat[d,2]+crit*sejk
      psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
      psihat[d,6]=1-psmm(abs(test[d,2]),CC,df)
    }
  }
  list(n=sam,test=test,psihat=psihat)
}

smmcrit<-function(nuhat,C){
  #
  #  Determine the .95 quantile of the C-variate Studentized maximum
  #  modulus distribution using linear interpolation on inverse
  #  degrees of freedom
  #  If C=1, this function returns the .975 quantile of Student's t
  #  distribution.
  #
  if(C-round(C)!=0)stop("The number of contrasts, C, must be an  integer")
  if(C>=29)stop("C must be less than or equal to 28")
  if(C<=0)stop("C must be greater than or equal to 1")
  if(nuhat<2)stop("The degrees of freedom must be greater than or equal to 2")
  if(C==1)smmcrit<-qt(.975,nuhat)
  if(C>=2){
    C<-C-1
    m1<-matrix(0,20,27)
    m1[1,]<-c(5.57,6.34,6.89,7.31,7.65,7.93,8.17,8.83,8.57,
              8.74,8.89,9.03,9.16,9.28,9.39,9.49,9.59, 9.68,
              9.77,9.85,9.92,10.00,10.07,10.13,10.20,10.26,10.32)
    m1[2,]<-c(3.96,4.43,4.76,5.02,5.23,5.41,5.56,5.69,5.81,
              5.92,6.01,6.10,6.18,6.26,6.33,6.39,6.45,6.51,
              6.57,6.62,6.67,6.71,6.76,6.80,6.84,6.88, 6.92)
    m1[3,]<-c(3.38,3.74,4.01,4.20,4.37,4.50,4.62,4.72,4.82,
              4.89,4.97,5.04,5.11,5.17,5.22,5.27,5.32, 5.37,
              5.41,5.45,5.49,5.52,5.56,5.59,5.63,5.66,5.69)
    m1[4,]<-c(3.09,3.39,3.62,3.79,3.93,4.04,4.14,4.23,4.31,
              4.38,4.45,4.51,4.56,4.61,4.66,4.70,4.74,4.78,
              4.82,4.85,4.89,4.92,4.95,4.98,5.00,5.03,5.06)
    m1[5,]<-c(2.92,3.19,3.39,3.54,3.66,3.77,3.86,3.94,4.01,
              4.07,4.13,4.18,4.23,4.28,4.32,4.36,4.39,4.43,
              4.46,4.49,4.52,4.55,4.58,4.60,4.63,4.65,4.68)
    m1[6,]<-c(2.80,3.06,3.24,3.38,3.49,3.59,3.67,3.74,3.80,
              3.86,3.92,3.96,4.01,4.05,4.09,4.13,4.16,4.19,
              4.22,4.25,4.28,4.31,4.33,4.35,4.38,4.39,4.42)
    m1[7,]<-c(2.72,2.96,3.13,3.26,3.36,3.45,3.53,3.60,3.66,
              3.71,3.76,3.81,3.85,3.89,3.93,3.96,3.99, 4.02,
              4.05,4.08,4.10,4.13,4.15,4.18,4.19,4.22,4.24)
    m1[8,]<-c(2.66,2.89,3.05,3.17,3.27,3.36,3.43,3.49,3.55,
              3.60,3.65,3.69,3.73,3.77,3.80,3.84,3.87,3.89,
              3.92,3.95,3.97,3.99,4.02,4.04,4.06,4.08,4.09)
    m1[9,]<-c(2.61,2.83,2.98,3.10,3.19,3.28,3.35,3.41,3.47,
              3.52,3.56,3.60,3.64,3.68,3.71,3.74,3.77,3.79,
              3.82,3.85,3.87,3.89,3.91,3.94,3.95, 3.97,3.99)
    m1[10,]<-c(2.57,2.78,2.93,3.05,3.14,3.22,3.29,3.35,3.40,
               3.45,3.49,3.53,3.57,3.60,3.63,3.66,3.69,3.72,
               3.74,3.77,3.79,3.81,3.83,3.85,3.87,3.89,3.91)
    m1[11,]<-c(2.54,2.75,2.89,3.01,3.09,3.17,3.24,3.29,3.35,
               3.39,3.43,3.47,3.51,3.54,3.57,3.60,3.63,3.65,
               3.68,3.70,3.72,3.74,3.76,3.78,3.80,3.82,3.83)
    m1[12,]<-c(2.49,2.69,2.83,2.94,3.02,3.09,3.16,3.21,3.26,
               3.30,3.34,3.38,3.41,3.45,3.48,3.50,3.53,3.55,
               3.58,3.59,3.62,3.64,3.66,3.68,3.69,3.71,3.73)
    m1[13,]<-c(2.46,2.65,2.78,2.89,2.97,3.04,3.09,3.15,3.19,
               3.24,3.28,3.31,3.35,3.38,3.40,3.43,3.46,3.48,
               3.50,3.52,3.54,3.56,3.58,3.59,3.61,3.63,3.64)
    m1[14,]<-c(2.43,2.62,2.75,2.85,2.93,2.99,3.05,3.11,3.15,
               3.19,3.23,3.26,3.29,3.32,3.35,3.38,3.40,3.42,
               3.44,3.46,3.48,3.50,3.52,3.54,3.55,3.57,3.58)
    m1[15,]<-c(2.41,2.59,2.72,2.82,2.89,2.96,3.02,3.07,3.11,
               3.15,3.19,3.22,3.25,3.28,3.31,3.33,3.36,3.38,
               3.39,3.42,3.44,3.46,3.47,3.49,3.50,3.52,3.53)
    m1[16,]<-c(2.38,2.56,2.68,2.77,2.85,2.91,2.97,3.02,3.06,
               3.09,3.13,3.16,3.19,3.22,3.25,3.27,3.29,3.31,
               3.33,3.35,3.37,3.39,3.40,3.42,3.43,3.45,3.46)
    m1[17,]<-c(2.35,2.52,2.64,2.73,2.80,2.87,2.92,2.96,3.01,
               3.04,3.07,3.11,3.13,3.16,3.18,3.21,3.23,3.25,
               3.27,3.29,3.30,3.32,3.33,3.35,3.36,3.37,3.39)
    m1[18,]<-c(2.32,2.49,2.60,2.69,2.76,2.82,2.87,2.91,2.95,
               2.99,3.02,3.05,3.08,3.09,3.12,3.14,3.17, 3.18,
               3.20,3.22,3.24,3.25,3.27,3.28,3.29,3.31,3.32)
    m1[19,]<-c(2.29,2.45,2.56,2.65,2.72,2.77,2.82,2.86,2.90,
               2.93,2.96,2.99,3.02,3.04,3.06,3.08,3.10, 3.12,
               3.14,3.16,3.17,3.19,3.20,3.21,3.23,3.24,3.25)
    m1[20,]<-c(2.24,2.39,2.49,2.57,2.63,2.68,2.73,2.77,2.79,
               2.83,2.86,2.88,2.91,2.93,2.95,2.97,2.98, 3.01,
               3.02,3.03,3.04,3.06,3.07,3.08,3.09,3.11,3.12)
    if(nuhat>=200)smmcrit<-m1[20,C]
    if(nuhat<200){
      nu<-c(2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,30,40,60,200)
      temp<-abs(nu-nuhat)
      find<-order(temp)
      if(temp[find[1]]==0)smmcrit<-m1[find[1],C]
      if(temp[find[1]]!=0){
        if(nuhat>nu[find[1]]){
          smmcrit<-m1[find[1],C]-
            (1/nu[find[1]]-1/nuhat)*(m1[find[1],C]-m1[find[1]+1,C])/
            (1/nu[find[1]]-1/nu[find[1]+1])
        }
        if(nuhat<nu[find[1]]){
          smmcrit<-m1[find[1]-1,C]-
            (1/nu[find[1]-1]-1/nuhat)*(m1[find[1]-1,C]-m1[find[1],C])/
            (1/nu[find[1]-1]-1/nu[find[1]])
        }
      }}
  }
  smmcrit
}


winvar <- function(x,tr=.2){
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr*n)+1
  itop <- length(x)-ibot+1
  xbot <- y[ibot]
  xtop <- y[itop]
  y <- ifelse(y<=xbot,xbot,y)
  y <- ifelse(y>=xtop,xtop,y)
  wv <- var(y)
  wv
}


listm<-function(x){
  #
  # Store the data in a matrix or data frame in a new
  # R variable having list mode.
  # Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
  #
  if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
  y<-list()
  for(j in 1:ncol(x))y[[j]]<-x[,j]
  y
}

#m2l=listm
#matrix2list=listm

qsmm<-function(q, r, nu) {
  #r=number of comparisons
  if (!is.finite(nu)) 
    return(qnorm(1 - 0.5 * (1 - q^(1/r))))
  res = uniroot(function(c, r, nu, q) {
    psmm(c, r = r, nu = nu) - q
  },
  c(0, 100), r = r, nu = nu, q = q)
  res$root
}


psmm = function(x, r, nu) {
  res = integrate(psmm.x, 0, Inf, c = x, r = r, nu = nu)
  res$value
}


psmm.x=function(x, c, r, nu) {
  snu = sqrt(nu)
  sx = snu * x
  lgx = log(snu) - lgamma(nu/2) + (1 - nu/2) * log(2) + 
    (nu - 1) * log(sx) + (-sx^2/2)
  exp(r * log(2 * pnorm(c * x) - 1) + lgx)
}
