library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ggsci)
library(repurrrsive)
library(broom)
library(minpack.lm)
library(AICcmodavg)
library(ggpubr)
library(latex2exp)
library(caret)

source("Functions/read_fi.R")
source("Functions/trna_concentration.R")
source("Functions/binding_models.R")
source("Functions/nls_reg_models.R")
source("Functions/ggplot_theme.R")
source("Functions/profile_lik_models.R")

fi_raw <- read_fi_raw(data.path = "RawData/FluorescenceSpec",
                       file.index =  c(2,4,5),
                       trna.conc = 30)

####...................data preparation.......................####
fi_raw_list <- fi_raw %>% 
  dplyr::select(-trna) %>% 
  as.list()

fi_baseline <- fi_raw_list %>% 
  map(~last(.x))

fi <- map2_df(fi_raw_list, fi_baseline, function(.x, .y) (.y - .x) / .y) %>% 
  rowwise() %>% 
  mutate(que_fi = mean(c(replica1,
                         replica2,
                         replica3)),
         stdev = sd(c(replica1,
                      replica2,
                      replica3))) %>% 
  ungroup() %>% 
  mutate(trna = fi_raw$trna) %>% 
  dplyr::select(trna, que_fi, stdev) %>% 
  filter(trna != 0)

saveRDS(fi, "Output/fi.rds")
####....................nls modeling..........................####
####bm1to1####
bm1to1_nls <- list(bm1to1 = 
                     nls(que_fi ~ bm1to1(fhg, k1, g0 = trna),
                     data = fi,
                     start = list(fhg = 0.164223,
                                   k1 = 0.280638)))
saveRDS(bm1to1_nls, file = "Output/model1.rds")                  
####bm1to2.deg.add####
bm1to2.deg.add_nls <- var_nls_deg(que_fi ~ bm1to2.deg.add(fhg, k1, trna),
                     data = fi,
                     start = list(fhg = 0.0944,
                                  k1 = 0.133),
        start.conc = fi$trna)
saveRDS(bm1to2.deg.add_nls, file = "Output/model2.rds")
####bm1to2.deg####
bm1to2.deg_nls <- var_nls_deg(que_fi ~ bm1to2.deg(fhg, fhg2, k1, trna),
                           data = fi,
                           start = list(fhg = 1.25,
                                        fhg2 = -8.0,
                                        k1 = 0.00355),
                   start.conc = fi$trna)
saveRDS(bm1to2.deg_nls, file = "Output/model3.rds")
####bm1to2.coop.add####
bm1to2.coop.add_nls <- var_nls_coop(que_fi ~ bm1to2.coop.add(fhg, k1, k2, trna),
                         data = fi,
                         start = list(fhg = 0.087,
                                      k1 = 5.9,
                                      k2 = 0.84),
                         start.conc = fi$trna)
saveRDS(bm1to2.coop.add_nls, file = "Output/model4.rds")
####bm1to2.coop####
bm1to2.coop_nls <- var_nls_coop(que_fi ~ bm1to2.coop(fhg, fhg2, k1, k2, trna),
                         data = fi,
                         start = list(fhg = 0.087,
                                      fhg2 = 0,
                                      k1 = 5.9,
                                      k2 = 0.84),
                         start.conc = fi$trna)
saveRDS(bm1to2.coop_nls, file = "Output/model5.rds")

####.......................plots.............................####
nls_reg <- list()
nls_reg$bm1to1_nls <- bm1to1_nls
nls_reg$bm1to1_nls$g <- fi$trna
nls_reg$bm1to2.deg.add_nls <- bm1to2.deg.add_nls
nls_reg$bm1to2.deg_nls <- bm1to2.deg_nls
nls_reg$bm1to2.coop.add_nls <- bm1to2.coop.add_nls
nls_reg$bm1to2.coop_nls <- bm1to2.coop_nls

summary(nls_reg)

nls_reg_pred <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map(~ predict(.x))

nls_reg_resid <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map(~ resid(.x))

nls_reg_trna <- nls_reg %>% 
  map(`[[`, 2)

nls_reg_names <- rep(c(names(nls_reg)), each = 12)
que_fi <- rep(fi$que_fi, times = 5)
stdev <- rep(fi$stdev, times = 5)

model_predict_tbl<- tibble(fi_pred = nls_reg_pred,
                           fi_resid = nls_reg_resid,
                           trna = nls_reg_trna) %>% 
  unnest(cols = c("fi_pred", "fi_resid", "trna")) %>% 
  mutate(model = nls_reg_names,
         que_fi = que_fi,
         stdev = stdev)
  

plot1 <- ggplot(model_predict_tbl, aes(x = trna, 
                                       y = que_fi, 
                                       color = model))+
  geom_point(size = 2.5)+
  geom_line(aes(y = fi_pred), size = 1)+
  geom_errorbar(aes(ymin = que_fi - stdev,
                    ymax = que_fi + stdev),
                    width = 0.3,
                    size = 0.8,
                color = "grey20",
                alpha = 0.3)+
  ylab(TeX("FI Quenching $\\[ \\% \\]$", bold = TRUE))+
  xlab(expression(bold(tRNA[Bs]^{Phe}~"[µM]")))+
  scale_color_uchicago()+
  mytheme_axes_grids+
  theme(legend.title = element_blank(),
        legend.position = c(0.7,0.35))
  
  

plot3 <- fi_raw %>% 
  pivot_longer(-trna, names_to = "Replicate", values_to = "fi") %>%
  mutate(fi = fi / 10^5) %>% 
  ggplot(aes(x = trna, y = fi, color = Replicate))+
  geom_point(size = 2.5)+
  ylab(TeX("FI $\\div 10^5$",
           bold = TRUE))+
  xlab(expression(bold(tRNA[Bs]^{Phe}~"[µM]")))+
  scale_color_nejm()+
  mytheme_axes_grids+
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.8))
  
ggarrange(plot3, plot1, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


plot2 <- ggplot(model_predict_tbl, aes(sample = fi_resid, color = model))+
  stat_qq(size = 2)+
  stat_qq_line(size = 1)+
  facet_grid(rows = vars(model))+
  scale_color_brewer(palette = "Dark2")

####Table#####################################################
deg <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map_dbl(~ nrow(fi) - length(coef(.x)))

rmse <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map_dbl(~sigma(.x))

fhg <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map_dbl(~ coef(.x)["fhg"])

fhg2 <- nls_reg %>% 
  map(`[[`, 1) %>% 
  map_dbl(~ coef(.x)["fhg2"])


nls_results_tbl <-
  as_tibble(matrix(nrow = 5, ncol = 8), 
            .name_repair = "minimal")
names(nls_results_tbl) <- c("Model", 
                            "df",
                            "RMSE",
                            "Fhg",
                            "Fhg2",
                            "beta12",
                            "K1",
                            "K2")
nls_results_tbl <- nls_results_tbl %>% 
  mutate(beta12 = as.double(beta12),
         K1 = as.double(K1),
         K2 = as.double(K2), 
         Model = c("bm1to1", 
                   "bm1to2.deg.add", 
                   "bm1to2.deg", 
                   "bm1to2.coop.add",
                   "bm1to2.coop"),
         df = deg,
         RMSE = rmse,
         Fhg = fhg,
         Fhg2 = fhg2)  
  



nls_results_tbl[2:3, "beta12"] <- c(coef(bm1to2.deg.add_nls$bm1to2.deg.add)["k1"],
                                    coef(bm1to2.deg_nls$bm1to2.deg)["k1"])

nls_results_tbl[4:5, "K1"] <- c(coef(bm1to2.coop.add_nls$bm1to2.coop.add)["k1"],
                               coef(bm1to2.coop_nls$bm1to2.coop)["k1"])
nls_results_tbl[4:5, "K2"] <- c(coef(bm1to2.coop.add_nls$bm1to2.coop.add)["k2"],
                                coef(bm1to2.coop_nls$bm1to2.coop)["k2"])

nls_results_tbl <- nls_results_tbl %>% 
  mutate_if(is.double, signif, digits = 3)

