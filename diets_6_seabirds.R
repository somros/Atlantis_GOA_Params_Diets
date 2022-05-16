# Alberto Rovellini
# April 4 2022

# Diets Part 6: Diet preferences for seabirds from Aydin et al. (2007)
# which in turn took it from Hunt et al. (2000)
# this is the starting point for seabird diets (pPREY matrix). Updating with more recent datasets would be good - or using recent data 
# to validate the realized diets

library(tidyverse)
library(data.table)
library(readxl)
library(viridis)

dat <- read_excel('../data/Diets GOA from Aydin et al. (2007).xlsx', sheet = 'Seabirds_calculations' ) 
fg <- read.csv('../data/GOA_Groups.csv') %>% select(Code,Name)

key <- dat[c(1,2),] %>% 
  t() %>% 
  data.frame() %>% 
  drop_na() %>% 
  set_names(c('Code','Biomass')) %>% 
  mutate(Species = rownames(.), Biomass = as.numeric(Biomass))

dat1 <- dat[-c(1,2),] %>%
  select(-Biomass) %>%
  mutate(across(c(where(is.character), -Atlantis_prey), as.numeric)) %>%
  group_by(Atlantis_prey) %>%
  summarize(across(everything(), sum)) %>% # add up the proportions of diet on the same prey for the same pred
  ungroup() %>%
  pivot_longer(cols=-Atlantis_prey, names_to = 'Predator', values_to = 'Prop') %>%
  left_join(key, by = c('Predator'='Species')) %>% # add the biomass of the predator species
  select(Code, Predator, Biomass, Atlantis_prey, Prop) %>%
  group_by(Code, Atlantis_prey) %>%
  mutate(meanprop = weighted.mean(Prop,Biomass)/100) %>% # mean of the prop for a pred group based on the propos in the species that compose it and their biomasses used as weights
  ungroup() %>%
  select(Code, Atlantis_prey, meanprop) %>%
  distinct() %>%
  set_names(c('Pred','Prey','Prop'))

# # view this
# dat1 %>% 
#   left_join(fg, by = c('Pred'='Code')) %>%
#   left_join(fg, by = c('Prey'='Code')) %>%
#   select(Name.x,Name.y,Prop) %>%
#   set_names('Predator','Prey','Proportion') %>%
#   ggplot()+
#   geom_tile(aes(x=Prey,y=Predator,fill=Proportion),color='darkgrey')+
#   scale_fill_viridis()+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1.05, vjust = 1, size = 14),
#         axis.text.y = element_text(size = 14))+
#   labs(title = "Proportion of predator diet", 
#        x = "Prey", y = "Predator",
#        fill = "Proportion of ingested prey")
  
# now expand to all prey codes, order, duplicate stages, and write out
# this does not capture ontogenetic preferences and at that point gape size becomes especially important
# keep an eye on the diet - if too much adult groundfish is eaten we will need to modify

pprey_seabirds <- expand.grid(Pred = unique(dat1$Pred), Prey = fg$Code) %>%
  left_join(dat1, by = c('Pred','Prey')) %>%
  mutate(Prop = replace_na(Prop, 0)) %>%
  arrange(factor(Pred, levels = fg$Code)) %>%
  mutate(Prop = Prop/10) %>% # divide by 10 as starting point
  pivot_wider(names_from = Prey, values_from = Prop) %>%
  slice(rep(1:n(), each = 4)) %>% # replicate 4 times (no ontogenetic preferences for now)
  mutate(DCsed = 0,
         DLsed = 0,
         DRsed = 0,
         PreyStage = rep(c(1,1,2,2), 4), # write stages 
         PredStage = rep(c(1,2,1,2), 4),
         name = paste0('pPREY',PreyStage,Pred,PredStage, ' ',(length(fg$Code)+3))) %>%
  select(name, KWT:DRsed)

# write out
write.csv(pprey_seabirds, '../output/seabirds_pprey.csv', row.names = F)
