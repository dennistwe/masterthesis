####..........generic functions and objects....................####

prof.seq <- function(lower, upper, len = 100){
  seq(lower, upper, length = len)
}
single.prof.plot <- function(model.list, mle.obj, theta){
  prof.tbl <- tibble(loglik = unlist(model.list$lik.profile[theta]),
                     theta = unlist(model.list$profile.vals[theta])
                     )
  ci95 <- mle.obj$value - (qchisq(0.95,1)/2)
  prof.tbl <- prof.tbl %>% 
    mutate(ci_bound = loglik >= ci95)
  
  ggplot(prof.tbl, aes(y = loglik, x = theta, color = ci_bound))+
    geom_segment(aes(x = theta, 
                     y = 0, 
                     xend = theta, 
                     yend = loglik))+
    geom_point()+
    scale_colour_manual(values = c("grey60", "black"))+
    theme(legend.position = "none")+
    geom_vline(xintercept = mle.obj$par[theta],
               color = "#0033FF", size = 1.1, lty = 2)+
    geom_hline(yintercept = ci95,
               color = "#CC0000", size = 1.1)+
    xlab(names(
      model.list$profile.vals[which(names(model.list$profile.vals) == theta)])
    )
}
contour.plot <- function(model.list, mle.obj, theta){
  conf95 <- mle.obj$value - qchisq(0.95,2)/2
  grid.vals <- data.frame(expand.grid(
    unlist(model.list$profile.vals[theta[1]]), 
    unlist(model.list$profile.vals[theta[2]])))
  prof.vals <- data.frame(model.list$lik.profile[paste(theta[1], 
                                                       theta[2], 
                                                       sep = ".")])
  names(prof.vals) <- "loglik"
  prof.contour <- cbind(grid.vals, prof.vals)
  
  ggplot(prof.contour, aes(y = Var1, 
                           x = Var2, 
                           z = loglik))+
    theme_minimal()+
    geom_tile(aes(fill = loglik))+
    stat_contour(color = "black")+
    geom_point(aes(y = mle.obj$par[theta[1]],
                   x = mle.obj$par[theta[2]]),
               color = "#0033FF")+
    scale_fill_continuous(low="grey70",high="grey15")+
    stat_contour(breaks=c(conf95), color = "#CC0000")+
    ylab(theta[1])+
    xlab(theta[2])
}
profile.vals <- list()
lik.profile <- list()
bm1to1_nls <- readRDS("Output/model1.rds")
bm1to2.deg.add_nls <- readRDS("Output/model2.rds")
bm1to2.deg_nls <- readRDS("Output/model3.rds")
bm1to2.coop.add_nls <- readRDS("Output/model4.rds")
bm1to2.coop_nls <- readRDS("Output/model5.rds")
fi <- readRDS("Output/fi.rds")

bm1to2.deg.add_nls$nls.coef <- c(coef(bm1to2.deg.add_nls$bm1to2.deg.add),
                             sig = sigma(bm1to2.deg.add_nls$bm1to2.deg.add))
bm1to2.deg_nls$nls.coef <- c(coef(bm1to2.deg_nls$bm1to2.deg),
                             sig = sigma(bm1to2.deg_nls$bm1to2.deg))
bm1to2.coop.add_nls$nls.coef <- c(coef(bm1to2.coop.add_nls$bm1to2.coop.add),
                                  sig = sigma(bm1to2.coop.add_nls$bm1to2.coop.add))
bm1to2.coop_nls$nls.coef <- c(coef(bm1to2.coop_nls$bm1to2.coop),
                              sig = sigma(bm1to2.coop_nls$bm1to2.coop))
mytheme_axes_box =
  theme_bw()+
  theme(axis.ticks.length=unit(.07, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 1))+
  theme(axis.text.y = element_text(color="black",size = 11))+
  theme(axis.text.x = element_text(color="black",size = 11))+
  theme(axis.title.y = element_text(size = 13,face = "bold"))+
  theme(axis.title.x = element_text(size = 13,face = "bold"))+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.border = element_rect(size = 1, color = "black"))

####..........maximum likelihood estimation....................####
####bm1to2.deg.add####
mle.deg.add <- optim(bm1to2.deg.add_nls$nls.coef,
                     ll.1to2.deg.add,
                     obs = fi,
                     g.curr = bm1to2.deg.add_nls$g,
                     control = list(fnscale=-1)
                     )
####bm1to2.deg####
mle.deg <- optim(bm1to2.deg_nls$nls.coef,
                 ll.1to2.deg,
                 obs = fi,
                 g.curr = bm1to2.deg_nls$g,
                 control = list(fnscale=-1))

####bm1to2.coop.add####
mle.coop.add <- optim(bm1to2.coop.add_nls$nls.coef,
                 ll.1to2.coop.add,
                 obs = fi,
                 g.curr = bm1to2.coop.add_nls$g,
                 control = list(fnscale=-1))
####bm1to2.coop####
mle.coop <- optim(bm1to2.coop_nls$nls.coef,
                      ll.1to2.coop,
                      obs = fi,
                      g.curr = bm1to2.coop_nls$g,
                      control = list(fnscale=-1))
####..........Single param profiles bm1to2.deg.add.............####
####k1#####
bm1to2.deg.add_nls$profile.vals$k1 <- prof.seq(0.05, 2)

bm1to2.deg.add_nls$lik.profile$k1 <- sapply(
  bm1to2.deg.add_nls$profile.vals$k1, function(k1) 
    optim(bm1to2.deg.add_nls$nls.coef,
          ll.1to2.deg.add,
          obs = fi,
          g.curr = bm1to2.deg.add_nls$g,
          profile.par = "k1",
          k1 = k1,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_deg1 <- single.prof.plot(bm1to2.deg.add_nls, mle.deg.add, "k1")+
  mytheme_axes_box+
  ylab(TeX("$LogLik$",
           bold = TRUE))+
  xlab(TeX("$\\beta_{12}$ $\\[ µm^{-2} \\]$",
           bold = TRUE))

bm1to2.deg.add_tbl <- tibble(
  loglik_vals = bm1to2.deg.add_nls$lik.profile$k1,
  param_vals = bm1to2.deg.add_nls$profile.vals$k1
)

ci95 <- mle.deg.add$value - (qchisq(0.95, 1) / 2)

bm1to2.deg.add_tbl <- bm1to2.deg.add_tbl %>% 
  mutate(ci_bound = loglik_vals >= ci95)

ggplot(bm1to2.deg.add_tbl, 
       aes(x = param_vals, 
           y = loglik_vals,
           color = ci_bound))+
  geom_segment(aes(x = param_vals, 
                   y = 0, 
                   xend = param_vals, 
                   yend = loglik_vals))+
  geom_point()+
  scale_colour_manual(values = c("grey60", "black"))+
  mytheme_axes_grids+
  theme(legend.position = "none")+
  ylab(TeX("$PL_{i}(p)$",
           bold = TRUE))+
  xlab(TeX("$\\theta_{i} = p_{i}$",
           bold = TRUE))+
  geom_hline(yintercept = ci95,
             color = "#CC0000", size = 1)+
  annotate("text", 
           x = 1.3, 
           y = 28.8,
           label = TeX("$log(L(\\hat{\\theta}_{nls}))$ - $\\frac{\\chi^2_{(0.95,1)}}{2}$",
                       bold = TRUE))
  
####fhg#####################################################

bm1to2.deg.add_nls$profile.vals$fhg <- prof.seq(0.01, 0.2)
 
bm1to2.deg.add_nls$lik.profile$fhg <- sapply(
  bm1to2.deg.add_nls$profile.vals$fhg, function(fhg) 
    optimParallel(bm1to2.deg.add_nls$nls.coef,
          ll.1to2.deg.add,
          obs = fi,
          g.curr = bm1to2.deg.add_nls$g,
          profile.par = "fhg",
          fhg = fhg,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_deg2 <- single.prof.plot(bm1to2.deg.add_nls, mle.deg.add, "fhg")+
  mytheme_axes_box+
  xlab(TeX("$F_{\\Delta HG}$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  ylab(TeX("LogLik",
           bold = TRUE))

####.............profile contour bm1to2.deg.add...............#####
####fhg.k1####
grid.deg.add.fhg.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.deg.add_nls$profile.vals$fhg), 
    unlist(bm1to2.deg.add_nls$profile.vals$k1))))
                                              
bm1to2.deg.add_nls$lik.profile$fhg.k1 <- sapply(
  grid.deg.add.fhg.k1, function(v) 
    optimParallel(bm1to2.deg.add_nls$nls.coef,
          ll.1to2.deg.add,
          obs = fi,
          g.curr = bm1to2.deg.add_nls$g,
          profile.par = c("fhg","k1"),
          fhg = v[1],
          k1 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_deg3 <- contour.plot(bm1to2.deg.add_nls, 
             mle.deg.add, 
             theta = c("fhg", "k1"))+
  ylab(TeX("$F_{\\Delta HG}$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  xlab(TeX("$\\beta_{12}$ $\\[ µm^{-2} \\]$",
           bold = TRUE))+
  mytheme_axes_box


####.............Single param profiles bm1to2.deg.............####
####k1####
bm1to2.deg_nls$profile.vals$k1 <- prof.seq(0.05, 2)
bm1to2.deg_nls$lik.profile$k1 <- sapply(
  bm1to2.deg_nls$profile.vals$k1, function(k1) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = "k1",
          k1 = k1,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.deg_nls, mle.deg, "k1")

####fhg####
bm1to2.deg_nls$profile.vals$fhg <- prof.seq(0.01, 0.3)
bm1to2.deg_nls$lik.profile$fhg <- sapply(
  bm1to2.deg_nls$profile.vals$fhg, function(fhg) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = "fhg",
          fhg = fhg,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.deg_nls, mle.deg, "fhg")

####fhg2####
bm1to2.deg_nls$profile.vals$fhg2 <- prof.seq(0.01, 0.3)
bm1to2.deg_nls$lik.profile$fhg2 <- sapply(
  bm1to2.deg_nls$profile.vals$fhg2, function(fhg2) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = "fhg2",
          fhg2 = fhg2,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.deg_nls, mle.deg, "fhg2")

####..............profile contours bm1to2.deg.................####
####fhg.fhg2####
grid.deg.fhg.fhg2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.deg_nls$profile.vals$fhg), 
    unlist(bm1to2.deg_nls$profile.vals$fhg2))))

bm1to2.deg_nls$lik.profile$fhg.fhg2 <- sapply(
  grid.deg.fhg.fhg2, function(v) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = c("fhg","fhg2"),
          fhg = v[1],
          fhg2 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.deg_nls, mle.deg, theta = c("fhg", "fhg2"))

####fhg.k1####
grid.deg.fhg.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.deg_nls$profile.vals$fhg), 
    unlist(bm1to2.deg_nls$profile.vals$k1))))

bm1to2.deg_nls$lik.profile$fhg.k1 <- sapply(
  grid.deg.fhg.k1, function(v) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = c("fhg","k1"),
          fhg = v[1],
          k1 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.deg_nls, mle.deg, theta = c("fhg", "k1"))

####fhg2.k1####
grid.deg.fhg2.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.deg_nls$profile.vals$fhg2), 
    unlist(bm1to2.deg_nls$profile.vals$k1))))

bm1to2.deg_nls$lik.profile$fhg2.k1 <- sapply(
  grid.deg.fhg2.k1, function(v) 
    optim(bm1to2.deg_nls$nls.coef,
          ll.1to2.deg,
          obs = fi,
          g.curr = bm1to2.deg_nls$g,
          profile.par = c("fhg2","k1"),
          fhg2 = v[1],
          k1 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.deg_nls, mle.deg, theta = c("fhg2", "k1"))

####..........Single param profiles bm1to2.coop.add.............####
####fhg####
bm1to2.coop.add_nls$profile.vals$fhg <- prof.seq(0.06, 0.32)

bm1to2.coop.add_nls$lik.profile$fhg <- sapply(
  bm1to2.coop.add_nls$profile.vals$fhg, function(fhg) 
    optim(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = "fhg",
          fhg = fhg,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.coop.add_nls, mle.coop.add, "fhg")

####k1####
bm1to2.coop.add_nls$profile.vals$k1 <- prof.seq(0.08, 0.4)

bm1to2.coop.add_nls$lik.profile$k1 <- sapply(
  bm1to2.coop.add_nls$profile.vals$k1, function(k1) 
    optim(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = "k1",
          k1 = k1,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.coop.add_nls, mle.coop.add, "k1")

####k2####
bm1to2.coop.add_nls$profile.vals$k2 <- prof.seq(0.001, 0.12)

bm1to2.coop.add_nls$lik.profile$k2 <- sapply(
  bm1to2.coop.add_nls$profile.vals$k2, function(k2) 
    optim(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = "k2",
          k2 = k2,
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.coop.add_nls, mle.coop.add, "k2")

####..........profile contours bm1to2.coop.add.................####
####fhg.k1####
grid.coop.add.fhg.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop.add_nls$profile.vals$fhg), 
    unlist(bm1to2.coop.add_nls$profile.vals$k1))))

bm1to2.coop.add_nls$lik.profile$fhg.k1 <- sapply(
  grid.coop.add.fhg.k1, function(v) 
    optim(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = c("fhg","k1"),
          fhg = v[1],
          k1 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop.add_nls, mle.coop.add, theta = c("fhg", "k1"))

####fhg.k2####
grid.coop.add.fhg.k2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop.add_nls$profile.vals$fhg), 
    unlist(bm1to2.coop.add_nls$profile.vals$k2))))

bm1to2.coop.add_nls$lik.profile$fhg.k2 <- sapply(
  grid.coop.add.fhg.k2, function(v) 
    optim(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = c("fhg","k2"),
          fhg = v[1],
          k2 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop.add_nls, mle.coop.add, theta = c("fhg", "k2"))


####k1.k2####
grid.coop.add.k1.k2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop.add_nls$profile.vals$k1), 
    unlist(bm1to2.coop.add_nls$profile.vals$k2))))

bm1to2.coop.add_nls$lik.profile$k1.k2 <- sapply(
  grid.coop.add.k1.k2, function(v) 
    optimParallel(bm1to2.coop.add_nls$nls.coef,
          ll.1to2.coop.add,
          obs = fi,
          g.curr = bm1to2.coop.add_nls$g,
          profile.par = c("k1","k2"),
          k1 = v[1],
          k2 = v[2],
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop.add_nls, mle.coop.add, theta = c("k1", "k2"))

####............Single param profiles bm1to2.coop.............####
####fhg####
bm1to2.coop_nls$profile.vals$fhg <- prof.seq(0.001, 0.017)
bm1to2.coop_nls$lik.profile$fhg <- sapply(
  bm1to2.coop_nls$profile.vals$fhg, function(fhg) 
    optim(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = "fhg",
          fhg = fhg,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_coop1 <- single.prof.plot(bm1to2.coop_nls, mle.coop, "fhg")+
  mytheme_axes_box+
  xlab(TeX("$F_{\\Delta HG}$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  ylab(TeX("LogLik",
           bold = TRUE))

####fhg2####
bm1to2.coop_nls$profile.vals$fhg2 <- prof.seq(0.15, 0.25)
bm1to2.coop_nls$lik.profile$fhg2 <- sapply(
  bm1to2.coop_nls$profile.vals$fhg2, function(fhg2) 
    optim(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = "fhg2",
          fhg2 = fhg2,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.coop_nls, mle.coop, "fhg2")

####k1####
bm1to2.coop_nls$profile.vals$k1 <- prof.seq(1000000, 5000000)
bm1to2.coop_nls$lik.profile$k1 <- sapply(
  bm1to2.coop_nls$profile.vals$k1, function(k1) 
    optim(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = "k1",
          k1 = k1,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_coop2 <- single.prof.plot(bm1to2.coop_nls, mle.coop, "k1")+
  mytheme_axes_box+
  xlab(TeX("$K_1$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  ylab(TeX("LogLik",
           bold = TRUE))+
  theme(axis.text.x = element_text(angle = 25))

####k2####
bm1to2.coop_nls$profile.vals$k2 <- prof.seq(0.06, 0.2)
bm1to2.coop_nls$lik.profile$k2 <- sapply(
  bm1to2.coop_nls$profile.vals$k2, function(k2) 
    optim(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = "k2",
          k2 = k2,
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

single.prof.plot(bm1to2.coop_nls, mle.coop, "k2")

####............profile contours bm1to2.coop...................####
####fhg.fhg2####
grid.coop.fhg.fhg2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$fhg), 
    unlist(bm1to2.coop_nls$profile.vals$fhg2))))

bm1to2.coop_nls$lik.profile$fhg.fhg2 <- sapply(
  grid.coop.fhg.fhg2, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("fhg","fhg2"),
          fhg = v[1],
          fhg2 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop_nls, mle.coop, theta = c("fhg", "fhg2"))

####fhg.k1####
grid.coop.fhg.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$fhg), 
    unlist(bm1to2.coop_nls$profile.vals$k1))))

bm1to2.coop_nls$lik.profile$fhg.k1 <- sapply(
  grid.coop.fhg.k1, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("fhg","k1"),
          fhg = v[1],
          k1 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

profile_coop3 <- contour.plot(bm1to2.coop_nls, 
             mle.coop, 
             theta = c("fhg", "k1"))+
  mytheme_axes_box+
  ylab(TeX("$F_{\\Delta HG}$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  xlab(TeX("$K_1$ $\\[ µm^{-1} \\]$",
           bold = TRUE))+
  theme(axis.text.x = element_text(angle = 25))

####fhg.k2####
grid.coop.fhg.k2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$fhg), 
    unlist(bm1to2.coop_nls$profile.vals$k2))))

bm1to2.coop_nls$lik.profile$fhg.k2 <- sapply(
  grid.coop.fhg.k2, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("fhg","k2"),
          fhg = v[1],
          k2 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop_nls, mle.coop, theta = c("fhg", "k2"))

####fhg2.k1####
grid.coop.fhg2.k1 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$fhg2), 
    unlist(bm1to2.coop_nls$profile.vals$k1))))

bm1to2.coop_nls$lik.profile$fhg2.k1 <- sapply(
  grid.coop.fhg2.k1, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("fhg2","k1"),
          fhg2 = v[1],
          k1 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop_nls, mle.coop, theta = c("fhg2", "k1"))

####fhg2.k2####
grid.coop.fhg2.k2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$fhg2), 
    unlist(bm1to2.coop_nls$profile.vals$k2))))

bm1to2.coop_nls$lik.profile$fhg2.k2 <- sapply(
  grid.coop.fhg2.k2, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("fhg2","k2"),
          fhg2 = v[1],
          k2 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop_nls, mle.coop, theta = c("fhg2", "k2"))

####k1.k2####
grid.coop.k1.k2 <- data.frame(
  t(expand.grid(unlist(
    bm1to2.coop_nls$profile.vals$k1), 
    unlist(bm1to2.coop_nls$profile.vals$k2))))

bm1to2.coop_nls$lik.profile$k1.k2 <- sapply(
  grid.coop.k1.k2, function(v) 
    optimParallel(bm1to2.coop_nls$nls.coef,
          ll.1to2.coop,
          obs = fi,
          g.curr = bm1to2.coop_nls$g,
          profile.par = c("k1","k2"),
          k1 = v[1],
          k2 = v[2],
          method = "L-BFGS-B",
          lower = c(0.001,0.001),
          control=list(fnscale = -1))$val)

contour.plot(bm1to2.coop_nls, mle.coop, theta = c("k1", "k2"))


ggarrange(profile_deg1, profile_coop1, 
          profile_deg2, profile_coop2, 
          profile_deg3, profile_coop3,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)
