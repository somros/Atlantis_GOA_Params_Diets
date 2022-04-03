# Alberto Rovellini
# April 1 2022

# Diets Part 5: Diet preferences for marine mammals and sharks from Aydin et al. (2007)

library(tidyverse)
library(readxl)
library(data.table)

# read in all data
fg <- read.csv('../data/GOA_Groups.csv') %>%
  select(Name,Code)

groups <- c('PIN','WHT','KWT','KWR','WHG','WHH','WHB','DOL','SHD','SHP')

datlist <- list()

for (i in 1:length(groups)){
  this_group <- groups[i]
  
  if(this_group=='PIN'){
    
    dat <- read_excel('../data/Diets GOA from Aydin et al. (2007).xlsx', sheet = this_group ) %>%
      select(Predator, Prey_Atlantis, Proportion)
    dat <- dat[-1,]
    dat <- dat %>% group_by(Predator, Prey_Atlantis) %>% 
      summarise(Prop = sum(Proportion)) %>% 
      ungroup()
    
  } else {
    
    dat <- read_excel('../data/Diets GOA from Aydin et al. (2007).xlsx', sheet = this_group ) %>%
      select(Predator, Prey_Atlantis, Proportion) %>% 
      group_by(Predator, Prey_Atlantis) %>% 
      summarise(Prop = sum(Proportion)) %>% 
      ungroup()
    
  }
  
  datlist[[i]] <- dat
}

dat <- rbindlist(datlist) %>% set_names(c('Predator','Prey','Prop'))

pprey_mammals_sharks <- expand.grid(Predator=groups, Prey = fg$Code) %>%
  left_join(dat, by = c('Predator','Prey')) %>%
  arrange(Predator) %>%
  mutate(Prop = replace_na(Prop, 0),
         Prop = Prop/100/10) %>% # starting point is prop / 10
  pivot_wider(names_from = Prey, values_from = Prop) %>%
  arrange(factor(Predator, levels = fg$Code)) %>%
  slice(rep(1:n(), each = 4)) %>% # replicate 4 times (no ontogenetic preferences for now)
  mutate(DCsed = 0, # add sediment detritus cols
         DLsed = 0,
         DRsed = 0,
         PreyStage = rep(c(1,1,2,2), length(groups)), # write stages 
         PredStage = rep(c(1,2,1,2), length(groups)),
         name = paste0('pPREY',PreyStage,Predator,PredStage, ' ',(length(fg$Code)+3))) %>%
  select(name, KWT:DRsed)

# write out
write.csv(pprey_mammals_sharks, '../output/mammals_sharks_pprey.csv', row.names = F)
