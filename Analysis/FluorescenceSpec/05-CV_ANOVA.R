library(dplyr)
library(tidyr)
library(ggplot2)
#remove.packages("Hmisc")
#install.packages("Hmisc")
#library(Hmsic)
library(MASS)
#crossval <- read.csv("Output/compCV.csv")
#crossval <- crossval[,-1]
compCV.tbl <- compCV.tbl%>% 
  mutate(fold_class = ifelse(grepl("1-2", fold_points),
                                  "outlier",
                                  "normal"),
         fold_class = factor(fold_class, levels = c("normal", 
                                                    "outlier")))

#custom_cont <- matrix(c(-4,1,1,1,1,
 #                       0,1,1,1,-3,
  #                      0,-1,-1,2,0,
   #                     0,-1,1,0,0), 
   #                  nrow = 5)
custom_cont <- matrix(c(-4/5,1/5,1/5,1/5,1/5,
                        0,1/4,1/4,1/4,-3/4,
                        0,-1/3,-1/3,2/3,0,
                        0,-1/2,1/2,0,0), 
                      nrow = 5)
colnames(custom_cont) <- c("c1", "c2", "c3", "c4")
rownames(custom_cont) <- c("bm1to1", "bm1to2.deg.add", 
                      "bm1to2.deg", "bm1to2.coop.add",
                      "bm1to2.coop")
###################################
# test for trimming sequence and taking mean vs mean(seq, trim)
cutoff <- 220 * 0.2 
floor(cutoff)
crossval %>% 
  group_by(model) %>% 
  arrange(model,RMSE) %>% 
  mutate(index = c(1:220)) %>% 
  filter(index > cutoff & index <= 220 - cutoff) %>% 
  summarise(reps = length(RMSE),
            avg = mean(RMSE))  
  
crossval %>% 
  group_by(model) %>%
  summarise(avg = mean(RMSE, trim = 0.2))

####################LinearReg and ANOVA##############################
contrasts(compCV.tbl$model) <- custom_cont
#contrasts(compCV.tbl$fold_class) <- c(1,-1)
compModel.lm <- lm(RMSE ~ 1, data = compCV.tbl)
augModel1.lm <- lm(RMSE ~ model, data = compCV.tbl)
#augModel2.lm <- lm(RMSE ~ fold_class + model, data = compCV.tbl)
#augModel3.lm <- lm(RMSE ~ fold_class + fold_class:model + model, data = compCV.tbl)
# this model not meaningful, investigates whether the CONTRASTS are different
#between the outlier class and the normal class, not a valuable information
augModel1.aov <- aov(augModel1.lm)
#augModel2.aov <- aov(augModel2.lm)
#augModel3.aov <- aov(augModel3.lm)

cont_summary1 <- summary(augModel1.aov,
        split = list(
          model = list(
            c1 = 1, c2 = 2, c3 = 3, c4 = 4)
        ))

# Residual Plot
augment(augModel.lm) %>% 
  ggplot(aes(y = RMSE, x = model))+
  geom_point()+
  geom_point(
    aes(y = .fitted, x = model), 
    color = "red",
    size = 3,
    shape = 12)+
  theme_classic()
#qq-plot
resid.B <- 
  augModel1.lm %>% 
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

resid.A <- augModel1.lm %>% 
  augment() %>% 
  ggplot(aes(x = model, y = .resid, fill = model))+
  geom_boxplot()+
  ylab("residuals")+
  mytheme_axes_box+
  scale_fill_uchicago()+
  theme(axis.text.x = element_blank(),
        legend.position = "none")
  
ggarrange(resid.A, resid.B,
          labels = c("A","B"),
          ncol = 2, nrow = 1)
####calculations of RSS, MSE, F for individual contrasts#####
num <- sum(xbar*custom_cont[,4])^2
denom <- sum((custom_cont[,4]^2)/sam)

ss <- apply(custom_cont, 2, function(x) (sum(xbar*x)^2) / (sum((x^2)/sam)))
(sum(xbar*custom_cont[,4]))/(sum(custom_cont[,4]^2)) # model params

anova(lm(RMSE ~ 1, data = crossval), augModel.lm)

# getting SSA from sum formula
ssa <- crossval %>% 
  group_by(model) %>% 
  mutate(ss = (RMSE - mean(RMSE))^2) %>% 
  summarise(sse = sum(ss)) %>% 
  summarise(msw = sum(sse)) %>% 
  pull(msw)

####defining PRE####
pre <- function(sum.aov){
  nRow <- nrow(sum.aov[[1]][2])
  ssr <- sum.aov[[1]][2][-nRow,]
  sseA <- sum.aov[[1]][2][nRow,]
  pre <- ssr / (ssr + sseA)
  return(pre)
}
  
#####defining ANOVA summary table####
aov.compModel <- tidy(anova(compModel.lm))
aov.df <- data.frame(Source = c("Model", "C1", "C2", "C3", "C4",
                                "Residuals", "Total"),
                     df = unlist(c(cont_summary1[[1]][1], 
                                   aov.compModel[1,"df"])),
                     SS = unlist(c(cont_summary1[[1]][2], 
                                       aov.compModel[1,"sumsq"])),
                     MS = unlist(c(cont_summary1[[1]][3], 
                                   aov.compModel[1,"meansq"])),
                     Fstat = unlist(c(cont_summary1[[1]][4],NA)),
                     p = unlist(c(cont_summary1[[1]][5],NA)),
                     PRE = c(pre(cont_summary1), NA, NA))
aov.df %>% 
  mutate(across(where(is.numeric), ~ round(., digits = 4)))
rownames(aov.df) <- NULL


# Finding outliers
crossval_ext <- read.csv("Output/compCV.csv")
crossval_ext <- crossval_ext[,-1]
crossval_ext <- crossval_ext %>% 
  mutate(model = factor(model, levels = c("bm1to1", 
                                          "bm1to2.deg.add", 
                                          "bm1to2.deg", 
                                          "bm1to2.coop.add",
                                          "bm1to2.coop")))

crossval_out <- crossval %>%   
group_by(model) %>% 
  filter(RMSE > quantile(RMSE)[4] + 1.5 * IQR(RMSE)) 

# use this plot finally
crossval_out %>% 
ggplot(aes(y = fold_points, x = model))+
  geom_point()

#####plot for CI based on contrast.ci fun#####
source("Functions/confint.R")
ci <- contrast.ci(augModel1.lm)

cont.ci.df <- data.frame(lower = ci[,1],
                      upper = ci[,2],
                      contrast = c("c1","c2","c3","c4"))

linconts <- diag(apply(custom_cont, 2, 
                       function(x) augModel1.lm$coefficients[-1] * sum(x^2)))
cont.ci.df %>% 
  pivot_longer(-contrast,names_to = "bound", values_to = "margin") %>% 
  mutate(coefficient = rep(linconts, each = 2)) %>% 
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

