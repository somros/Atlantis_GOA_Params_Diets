# This computes plankton diets 

library(tidyverse)

plankton <- read.csv('../data/plankton.csv')

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
  summarise(Prop = mean(Prop))

# view
plankton3 %>%
  ggplot(aes(y=Prop,x=Predator_Atlantis,fill=Prey_Atlantis))+
  geom_bar(position = 'stack',stat='identity')
