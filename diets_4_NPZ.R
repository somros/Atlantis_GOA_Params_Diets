# Alberto Rovellini
# March 25 2022

# Diets Part 4: Diet preferences for plankton from NPZ

# This code takes plankton diet preferences from Coyle et al. (2019) - which is an NPZ model for the GOA 
# The code writes out the rows of the pprey matrix for ZM, ZS, and EUP. Those will then be processed by write_pprey.R

library(tidyverse)

plankton <- read.csv('../data/plankton.csv')
fg <- read.csv('../data/GOA_Groups.csv') %>% select(Code,Name)

# for each predator from NPZ, turn those numbers in proportions

plankton <- plankton %>%
  group_by(Predator) %>%
  mutate(Prop = Value/sum(Value)) %>%
  ungroup()

# group by atlantis predator and take averages of the proportions of NPZ prey belonging to the same Atlantis prey

plankton1 <- plankton %>%
  group_by(Predator, Predator_Atlantis, Prey_Atlantis) %>%
  summarise(Prop = sum(Prop)) %>% ungroup()

plankton2 <- expand.grid(unique(plankton1$Predator), unique(plankton1$Prey_Atlantis)) %>%
  set_names(c('Predator','Prey_Atlantis')) %>%
  arrange(Predator) %>%
  left_join((plankton1 %>% select(-Predator_Atlantis)), by = c('Predator','Prey_Atlantis')) %>%
  left_join((plankton1 %>% select(Predator, Predator_Atlantis) %>% distinct()), by = 'Predator') %>%
  mutate(Prop = replace_na(Prop, 0))

plankton3 <- plankton2 %>%
  group_by(Predator_Atlantis, Prey_Atlantis) %>%
  summarise(Prop = mean(Prop)) %>%
  ungroup()

# view
plankton3 %>%
  ggplot(aes(y=Prop,x=Predator_Atlantis,fill=Prey_Atlantis))+
  geom_bar(position = 'stack',stat='identity')

# we are missing macrozooplankton information here, but we can get those from Kerim's tech memo (they were quite generic)
# and then we can calibrate

# Euphausiids align OK to Kerim's model

# write a pprey matrix with rows for EUP, ZM, ZS
# first, change names to codes in plankton3
plankton4 <- plankton3 %>%
  left_join(fg, by = c('Predator_Atlantis'='Name')) %>%
  left_join(fg, by = c('Prey_Atlantis'='Name')) %>%
  set_names(c('predname','preyname','prop','pred','prey')) %>%
  select(pred,prey,prop)

# make pprey
pprey_plankton <- expand.grid(unique(plankton4$pred), fg$Code) %>% 
  set_names(c('pred','prey')) %>%
  left_join(plankton4, by = c('pred','prey')) %>%
  mutate(prop = replace_na(prop, 0),
         prop = prop/10) %>% # initial estimate like in CalCUrrent
  pivot_wider(names_from = prey, values_from = prop) %>%
  mutate(DCsed = 0, 
         DLsed = 0, 
         DRsed = 0, #add columns for detritus
         name = paste0('pPREY',pred,' ',81)) %>% 
  select(name,KWT:DRsed)
  
# write out
write.csv(pprey_plankton, '../output/plankton_pprey.csv', row.names = F)
