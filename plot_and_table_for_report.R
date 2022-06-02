# plot and produce tables for the methods
library(tidyverse)
library(viridis)

diets <- read.csv('../output/goa_pprey_matrix.csv')

# plot
diets <- diets %>%
  mutate(tmp = substr(name, 6, (nchar(name)-3))) %>%
  rowwise() %>%
  mutate(Prey_stage = ifelse(substr(tmp,1,1) == '1', 1, ifelse(substr(tmp,1,1) == '2', 2, NA)),
         Pred_stage = ifelse(substr(tmp,nchar(tmp),nchar(tmp)) == '1', 1, ifelse(substr(tmp,nchar(tmp),nchar(tmp)) == '2', 2, NA)),
         Pred_name = ifelse(substr(tmp,1,1) %in% c('1','2'), substr(tmp,2,4), tmp)) %>%
  select(Pred_name, Pred_stage, Prey_stage, KWT:DRsed) 

# reconstruct invertebrate diets
diets_inv_tmp <- diets %>% filter(is.na(Pred_stage))

diets_inv <- diets_inv_tmp[rep(seq_len(nrow(diets_inv_tmp)), each = 4), ]
diets_inv$Pred_stage <- rep(c(1,2), nrow(diets_inv_tmp)*2)
diets_inv$Prey_stage <- rep(c(1,1,2,2), nrow(diets_inv_tmp))

# and now bring together
diets <- rbind(diets %>% filter(!is.na(Pred_stage)), diets_inv)

diets_long <- diets %>%
  pivot_longer(-(Pred_name:Prey_stage), names_to = 'Prey_name', values_to = 'Prop') %>%
  mutate(Prop = Prop * 10 * 100,
         Stage = paste0('Prey', Pred_stage, ':Predator', Prey_stage)) 

# attach long names
fg <- read.csv('../data/GOA_Groups.csv') %>% select(Code, Name)

diets_long <- diets_long %>%
  left_join(fg, by = c('Pred_name' = 'Code')) %>%
  left_join(fg, by = c('Prey_name' = 'Code')) %>%
  select(Name.x, Name.y, Stage, Prop) %>%
  rename(Pred_name = Name.x, Prey_name = Name.y) %>%
  mutate(Prop = na_if(Prop, 0)) %>%
  filter(!is.na(Prey_name))

p <- diets_long %>% ggplot()+
  geom_tile(aes(x = Prey_name, y = Pred_name, fill = Prop), color = 'darkgrey')+
  scale_fill_viridis()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1.05, vjust = 1, size = 11),
        axis.text.y = element_text(size = 11))+
  facet_wrap(~Stage)+
  labs(x = "Prey", y = "Predator",
       fill = "Ingested\nprey (%)")

p

ggsave('complete_dietmatrix.png', p, width = 19, height = 18, units = 'in', dpi = 300)

# really hard to read as a plot...

# make a table

diet_tab <- diets %>%
  left_join(fg, by = c('Pred_name' = 'Code')) %>%
  mutate(Pred_stage = str_replace(Pred_stage, '1', 'J'),
         Pred_stage = str_replace(Pred_stage, '2', 'A'),
         Prey_stage = str_replace(Prey_stage, '1', 'J'),
         Prey_stage = str_replace(Prey_stage, '2', 'A')) %>%
  select(Name, Pred_stage, Prey_stage, KWT:DRsed) %>%
  set_names(c('Predator','Predator stage','Prey stage'), 
            fg$Name, 'Carrion_sediments', 'Detritus_labile_sediment', 'Detritus_refractory_sediment') %>%
  mutate(across(where(is.numeric), ~ .x*10*100))

# write as csv
write.csv(diet_tab, 'diets_table.csv', row.names = F)
