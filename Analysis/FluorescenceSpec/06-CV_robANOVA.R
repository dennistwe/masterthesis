# robust anova
source("Functions/robust_ANOVA.R")
# bring tidy tibble into list mode
compCV.ls <- compCV.tbl %>% 
  dplyr::select(c(model, RMSE)) %>% 
  mutate(index = rep(c(1:220), times = 5)) %>% 
  pivot_wider(names_from = model, values_from = RMSE) %>% 
  dplyr::select(-index) %>% 
  as.list()


robustComp <- lincon(compCV.ls, custom_cont)

ci.rob.df <- data.frame(lower = robustComp$psihat[,3],
                        upper = robustComp$psihat[,4],
                        contrast = c("c1","c2","c3","c4"))

ci.B <- ci.rob.df %>% 
  pivot_longer(-contrast,names_to = "bound", values_to = "margin") %>% 
  mutate(coefficient = rep(robustComp$psihat[,2], each = 2)) %>% 
  ggplot(aes(x = margin, y = contrast))+
  geom_point(size = 2)+
  geom_line(size = 1.0)+
  geom_point(aes(x = coefficient), color = "blue", size = 3, shape = 18)+
  geom_vline(xintercept = 0, lty = 2, size = 1.0,
             color = "red")+
  theme_bw()+
  theme(axis.ticks.length=unit(.07, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 1))+
  theme(axis.text.y = element_text(color="black",size = 12))+
  theme(axis.text.x = element_text(color="black",size = 12))+
  theme(axis.title.y = element_text(size = 14,face = "bold"))+
  theme(axis.title.x = element_text(size = 14,face = "bold"))+
  theme(panel.grid.minor.x = element_line(color = "grey80"),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(size = 1.2, color = "black"))


####anova on trimmed data set######
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

compCV.trim.tbl %>% 
  ggplot(aes(y = RMSE, x = model, fill = model))+
  geom_boxplot()
#contrasts(compCV.trim.tbl$model) <- custom_cont
augModel1.trim.lm <- lm(RMSE ~ model, data = compCV.trim.tbl)
augModel1.trim.aov <- aov(augModel1.trim.lm)
cont_summary1.trim <- summary(augModel1.trim.aov,
                           split = list(
                             model = list(
                               c1 = 1, c2 = 2, c3 = 3, c4 = 4)
                           ))


psihat <- robustComp$psihat[,2]
csquare <- apply(custom_cont, 2, function(x) sum(x^2))
psihat / csquare

resid.B <- augModel1.trim.lm %>% 
  augment() %>% 
  ggplot(aes(sample = .resid))+
  stat_qq(aes(color = model), alpha = 0.7, size = 1.3)+
  stat_qq_line(color = "black", size = 1.3)+
  scale_color_uchicago()+
  mytheme_axes_box+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10,
                                   face = "bold"),
        legend.position = c(0.25,0.72))


contrast.ci <- function(model.lm, level = .05){
  ci_margin <- function(x){
    omega <- sum(x^2 / n)
    err <- sqrt(fcrit * msw * omega)
    err
  }
  rss <- resid(model.lm)
  rssDF <- model.lm$df.residual
  msw <- t(rss) %*% rss  / rssDF
  fcrit <- qf(level, 1, rssDF, lower.tail = FALSE)
  model.coef <- coef(model.lm)[-1]
  #cont_mat <- model.matrix(model.lm)[,-1]
  #N <- dim(cont_mat)[1]
  n <- 220
  csquare <- apply(custom_cont, 2, function(x) sum(x^2))
  ci_mat <- matrix(NA, length(model.coef), 2L,
                   dimnames = list(names(model.coef), c("lower", "upper")))
  cont_se <- apply(custom_cont, 2, ci_margin)
  ci_mat[,1] <- model.coef * csquare - cont_se
  ci_mat[,2] <- model.coef * csquare + cont_se
  ci_mat
}

#####rob.ci.df#####
rob.ci.df <- data.frame(Source = c("C1", "C2", "C3", "C4"),
                        psihat = c(robustComp$psihat[,2]),
                        df = robustComp$test[,5],
                        tstat = robustComp$test[,2],
                        se = robustComp$test[,4],
                        p = robustComp$psihat[,6])
