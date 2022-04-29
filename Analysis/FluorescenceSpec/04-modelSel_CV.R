library(vtreat)
library(dplyr)
library(tidyr)
library(purrr)
#install.packages("Hmisc")
library(Hmisc)

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

saveRDS(compCV.tbl, file = "Output/compCV.rds")
#write.csv(compCV.tbl, "Output/compCV.csv")

compCV_bar.tbl <- compCV.tbl %>% 
  group_by(model) %>% 
  summarise(cv_RMSE = mean(RMSE)) %>% 
  ungroup() %>% 
  mutate(nls_RMSE = model_sigma) %>% 
  mutate(diff = cv_RMSE - nls_RMSE)

###############Plots#####################
posn.d <- position_dodge(0.9)
plotA <- ggplot(compCV.tbl, aes(y = RMSE, 
                                x = model))+
  geom_jitter(width = 0.3, color = "grey50", alpha = 0.6)+
  stat_summary(fun.y = mean, geom = "point", 
               position = posn.d,
               fill = "red", size = 3.5, shape = 23)+
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "errorbar", width = 0, 
               position = posn.d, size = 1.1)+
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
  theme(legend.position = c(0.9, 0.9),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 15))

ggplot(compCV.tbl, aes(x = RMSE, fill = model))+
  geom_histogram(bins = 25)+
  facet_wrap(~model)

cvBox <- ggplot(compCV.tbl, aes(y = RMSE, 
                       x = model, 
                       fill = model))+
  geom_boxplot()+
  mytheme_axes_box+
  scale_fill_uchicago()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15,
                                   vjust = 0.6))
# outlier plot
compCV_out <- compCV.tbl %>%   
  group_by(model) %>% 
  filter(RMSE > quantile(RMSE)[4] + 1.5 * IQR(RMSE)) 

outPlot <- compCV_out %>% 
  ggplot(aes(y = fold_points, x = model))+
  geom_point()+
  ylab("test fold")+
  mytheme_axes_box+
  theme(axis.text.x = element_text(angle = 15,
                                   vjust = 0.6))
  
ggarrange(cvBox, outPlot,
          labels = c("A","B"),
          nrow = 1, ncol = 2)


is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

compCV.out.tbl <- compCV.tbl %>% 
  group_by(model) %>% 
  mutate(is_outlier = ifelse(is_outlier(RMSE), RMSE, as.numeric(NA)))
compCV.out.tbl$fold_points[which(is.na(compCV.out.tbl$is_outlier))] <- as.numeric(NA)


pA <- ggplot(compCV.out.tbl, aes(y=RMSE, x=model, fill = model)) + 
  geom_boxplot() + 
  geom_text(aes(label=fold_points),
            na.rm=TRUE,
            hjust = -0.2,
            angle = 10,
            size = 3.2,
            check_overlap = TRUE)+
  mytheme_axes_box+
  scale_fill_uchicago()+
  theme(axis.text.x = element_text(angle = 15,
                                   vjust = 0.6),
        legend.position = "none")

trim <- 0.2
n <- 220
cutoff <- floor(n * trim) 

compCV.trim.tbl<- compCV.tbl %>% 
  group_by(model) %>% 
  arrange(model,RMSE) %>% 
  mutate(index = c(1:220)) %>% 
  filter(index > cutoff & index <= 220 - cutoff) %>%
  dplyr::select(-index) %>% 
  ungroup()

pB <- compCV.trim.tbl %>% 
  ggplot(aes(y = RMSE, x = model, fill = model))+
  geom_boxplot()+
  mytheme_axes_box+
  scale_fill_uchicago()+
  theme(axis.text.x = element_text(angle = 25,
                                   vjust = 0.6))



