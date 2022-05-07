library(bookdown)
library(minpack.lm)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(ggsci)
library(latex2exp)
library(AICcmodavg)
library(vtreat)
library(kableExtra)
library(broom)

source("Functions/binding_models.R", local = knitr::knit_global())
source("Functions/ggplot_theme.R", local = knitr::knit_global())
bm1to1_nls <- readRDS("Output/model1.rds")
bm1to2.deg.add_nls <- readRDS("Output/model2.rds")
bm1to2.deg_nls <- readRDS("Output/model3.rds")
bm1to2.coop.add_nls <- readRDS("Output/model4.rds")
bm1to2.coop_nls <- readRDS("Output/model5.rds")
fi <- readRDS("Output/fi.rds")

k_fold_cross <- function(data,
                         model.list, 
                         model.fun = c("bm1to1", 
                                       "bm1to2.deg.add", 
                                       "bm1to2.deg", 
                                       "bm1to2.coop.add",
                                       "bm1to2.coop"),
                         nFolds){
  select.model <- match.arg(model.fun)
  nRow <- nrow(data)
  model.formula <- formula(model.list[[select.model]])
  model.coef <- as.list(coef(model.list[[select.model]]))
  model.list$cross_valid <- list()
  #############################################################
  # Set up list of all possible split combinations
  fold_comb <- combn(nRow, nRow / nFolds)
  sublist_app <- map(seq_len(ncol(fold_comb)), ~fold_comb[,.])
  
  split_list <- rep(list(list(train = vector("integer", length = 9L),
                              app = vector("integer", length = 3L))), 
                    times = ncol(fold_comb))
  
  split_list <- split_list %>% 
    map2(sublist_app, ~list(.x$app <- .y, 
                            .x$train <- seq(nRow)[- .y])) 
  
  for(i in seq_along(1:ncol(fold_comb))){
    names(split_list[[i]]) <- c("app", "train")
  }
  
  predict.vals <- vector(mode = "double", length = nRow / nFolds)
  
  if (select.model != "bm1to1"){
    data$trna <- model.list$g
    names(data) <- c("g.best", "que_fi")
  } 
  
  for(i in seq_along(1:ncol(fold_comb))){
    split <- split_list[[i]] 
    model <- nlsLM(model.formula, 
                   data = data[split$train,],
                   start = model.coef)
    predict.vals <- predict(model, newdata = data[split$app,])
    
    obs.vals <- data[split$app,"que_fi"]
    m <- length(predict.vals)
    fold_rmse <- sqrt(sum((obs.vals - predict.vals)^2) / m)
    model.list$cross_valid[[i]] <- fold_rmse
  }
  return(model.list)
}


fi <- fi[,1:2]
bm1to1_nls <- k_fold_cross(data = fi,
                           model.list = bm1to1_nls,
                           model.fun = c("bm1to1"),
                           nFolds = 4)

bm1to2.deg_nls <- k_fold_cross(
  data = fi, 
  model.list = bm1to2.deg_nls, 
  model.fun = c("bm1to2.deg"),
  nFolds = 4)

bm1to2.deg.add_nls <- k_fold_cross(
  data = fi, 
  model.list = bm1to2.deg.add_nls, 
  model.fun = c("bm1to2.deg.add"),
  nFolds = 4)

bm1to2.coop.add_nls <- k_fold_cross(
  data = fi, 
  model.list = bm1to2.coop.add_nls, 
  model.fun = c("bm1to2.coop.add"),
  nFolds = 4)

bm1to2.coop_nls <- k_fold_cross(
  data = fi, 
  model.list = bm1to2.coop_nls, 
  model.fun = c("bm1to2.coop"),
  nFolds = 4)

cross_valid.ls <- list(bm1to1_nls,
                       bm1to2.deg.add_nls,
                       bm1to2.deg_nls,
                       bm1to2.coop.add_nls,
                       bm1to2.coop_nls)
names(cross_valid.ls) <- c("bm1to1", 
                           "bm1to2.deg.add", 
                           "bm1to2.deg", 
                           "bm1to2.coop.add",
                           "bm1to2.coop")

cross_error <- cross_valid.ls %>% 
  map("cross_valid") %>% 
  flatten()

model_sigma <- cross_valid.ls %>% 
  map(`[[`, 1) %>% 
  map_dbl(~ sigma(.x))

nRow <- 12
nFolds <- 4
fold_comb <- combn(nRow, nRow / nFolds)
sublist_app <- map(seq_len(ncol(fold_comb)), ~fold_comb[,.])
sublist_app_chr <- sublist_app %>% 
  map(`[`, c(1,2,3)) %>%
  map_chr(~paste0(.x, collapse = "-")) 

compCV.tbl <- tibble(model = factor(rep(names(cross_valid.ls), each = 220), 
                                    levels = names(cross_valid.ls)),
                     RMSE = unlist(cross_error),
                     fold = rep(c(1:220), times = 5),
                     fold_points = rep(sublist_app_chr, times = 5))

compCV_bar.tbl <- compCV.tbl %>% 
  group_by(model) %>% 
  summarise(cv_RMSE = mean(RMSE)) %>% 
  ungroup() %>% 
  mutate(nls_RMSE = model_sigma)

###############Plots#####################
posn.d <- position_dodge(0.9)
plotA <- ggplot(compCV.tbl, aes(y = RMSE, 
                                x = model))+
  geom_jitter(width = 0.3, color = "grey50", alpha = 0.6)+
  
  #stat_summary(fun.data = "mean_sdl", 
   #            fun.args = list(mult = 1), 
    #           geom = "errorbar", width = 0, 
     #          position = posn.d, size = 1.1)+
  stat_summary(fun = mean,
               geom = "pointrange",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = posn.d,
               shape = 23,
               size = 1)+
  stat_summary(fun.y = mean, geom = "point", 
               position = posn.d,
               fill = "red", size = 3.5, shape = 23)+
  mytheme_axes_grids+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 15))

plotB <- compCV_bar.tbl %>% 
  pivot_longer(cols = - model, 
               names_to = "Prediction",
               values_to = "RMSE") %>% 
  mutate(Prediction = factor(Prediction,
                             levels = c("nls_RMSE", "cv_RMSE"))) %>%
  
  ggplot(aes(x = model, y = RMSE, fill = Prediction))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_brewer(palette = "Set1")+
  mytheme_axes_grids+
  theme(legend.position = c(0.8, 0.9),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 15))

ggarrange(plotA, plotB, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

