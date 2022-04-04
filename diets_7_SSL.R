# Alberto Rovellini
# April 4 2022

# Diets Part 7: Diet preferences for Steller sea lions

# The values are averages from Aydin et al (2007), Trites et al. (2007), Tollit et al (2015, 2017)
# Key thing to remember is that the point of those (and other) studies is that SSL diet varies widely between regions (W-E), season, sex, and pregnancy status
# We do need some average for the PPREY matrix, but we will need to check that the realized diets make sense in space and time
# For the availability (PPREY) matrix, we take several averages: W and C in Aydin et al. (2007), over seasons for Trites et al. (2007),
# in time and space for Tollit et al. (2015, 2017 - although the authors had calculated those), and between EC-V and EC-F for Tollit
# et al. (2015). We then average across studies.
# Also values have different meanings in these studies: Aydin et al. (2007) reports consumption, Trites frequency of occurrence, and 
# Tollit energetic content (and many other metrics). So this is to be taken as a starting point - precise calculations would likely
# not be worth the squeeze since the PPREY is constant in time and space and has to be calibrated

library(tidyverse)
library(data.table)
library(readxl)
library(viridis)

dat <- read_excel('../data/Diets GOA from Aydin et al. (2007).xlsx', sheet = 'SSL_calculations', range = 'K1:P28')
fg <- read.csv('../data/GOA_Groups.csv') %>% select(Code,Name)

dat1 <- dat %>% 
  select(Prey_Atlantis,Mean) %>%
  set_names('Prey','Prop') %>%
  mutate(Pred = 'SSL')

# view
# dat1 %>%
#   left_join(fg, by = c('Pred'='Code')) %>%
#   left_join(fg, by = c('Prey'='Code')) %>%
#   select(Name.x,Name.y,Prop) %>%
#   set_names('Predator','Prey','Proportion') %>%
#   ggplot()+
#   geom_bar(aes(x=Prey,y=Proportion), stat = 'identity')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1.05, vjust = 1, size = 14),
#         axis.text.y = element_text(size = 14))+
#   labs(title = "Percentage of predator diet")

pprey_ssl <- expand.grid(Pred = unique(dat1$Pred), Prey = fg$Code) %>%
  left_join(dat1, by = c('Pred','Prey')) %>%
  mutate(Prop = replace_na(Prop, 0),
         Prop = Prop/100/10) %>% # turn percentage to proportion and then divide by 10
  pivot_wider(names_from = Prey, values_from = Prop) %>%
  slice(rep(1:n(), each = 4)) %>% # replicate 4 times (no ontogenetic preferences for now)
  mutate(DCsed = 0,
         DLsed = 0,
         DRsed = 0,
         PreyStage = c(1,1,2,2), # write stages 
         PredStage = c(1,2,1,2),
         name = paste0('pPREY',PreyStage,Pred,PredStage, ' ',(length(fg$Code)+3))) %>%
  select(name, KWT:DRsed)

# write out
write.csv(pprey_ssl, '../output/ssl_pprey.csv', row.names = F)
