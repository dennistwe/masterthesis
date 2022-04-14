####..........model performance analysis by aicc..............####
model.list <- list(bm1to1 = bm1to1_nls,
                   bm1to2.deg.add = bm1to2.deg.add_nls,
                   bm1to2.deg = bm1to2.deg_nls,
                   bm1to2.coop.add = bm1to2.coop.add_nls,
                   bm1to2.coop = bm1to2.coop_nls)

aicc_score <- model.list %>%
  map(`[[`, 1) %>% 
  map_dbl(~ AICc(.x)) 

k <- model.list %>% 
  map(`[[`, 1) %>% 
  map_dbl(~ length(coef(.x))) 
  

model.aicc <- tibble(model = names(model.list),
                     k = k +1,
                     aicc = aicc_score)

model.aicc <- model.aicc %>% 
  mutate(del_aicc = aicc - min(aicc),
         aicc_wt = (exp(-del_aicc/2)/sum(exp(-del_aicc/2)))*100) %>% 
  arrange(del_aicc) %>% 
  mutate(cum_wt = cumsum(aicc_wt))

model.aicc.best <- model.aicc[1,]

model.aicc.reality <- model.aicc %>% 
  select(model, aicc)
model.aicc.reality[6,1] <- c("reality")
model.aicc.reality[6,2] <- c(-60)

only.reality <- model.aicc.reality[6,]
only.best.model <- model.aicc.reality %>% 
  filter(model == "bm1to2.deg.add")

#aictab(model.list)

delta.aicc.plot <- model.aicc %>% 
  ggplot(aes(y = del_aicc, x = model))+
  geom_point(size = 3.5, 
             shape = 23, fill = "black",
             color = "black")+
  geom_point(data = model.aicc.best,
             fill = "royalblue4",
             size = 3.5, 
             shape = 23,)+
  coord_polar()+
  theme_minimal()+
  theme(axis.ticks.length=unit(.07, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 1))+
  theme(axis.text.y = element_text(color="black",size = 12))+
  theme(axis.text.x = element_text(color="black",size = 9,face = "bold"))+
  theme(axis.title.y = element_text(size = 14,face = "bold"))+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_line(colour = "grey70",
                                        size = 0.8))+
  ylab(TeX("$\\Delta AIC_c$", bold = TRUE))


aicc.plot <- model.aicc.reality %>% 
  ggplot(aes(y = aicc, x = model))+
  geom_point(size = 3.5, 
             shape = 23, fill = "black",
             color = "black")+
  geom_point(data = only.reality,
             fill = "red3",
             size = 3.5, 
             shape = 23,)+
  geom_point(data = only.best.model,
             fill = "royalblue4",
             size = 3.5, 
             shape = 23,)+
  coord_polar()+
  theme_minimal()+
  theme(axis.ticks.length=unit(.07, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 1))+
  theme(axis.text.y = element_text(color="black",size = 12))+
  theme(axis.text.x = element_text(color="black",size = 9,face = "bold"))+
  theme(axis.title.y = element_text(size = 14,face = "bold"))+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_line(colour = "grey70",
                                        size = 0.8))+
  ylab(TeX("$AIC_c$", bold = TRUE))

ggarrange(aicc.plot, delta.aicc.plot, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
