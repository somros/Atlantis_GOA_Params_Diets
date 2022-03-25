# Alberto Rovellini
# March 25 2022

# Diets part 3: Diet preferences for invertebrates from Ecopath

# This code reads invertebrate diets information from Aydin et al. (2007) and wriets out the lines of the PPREY matrix
# for the inverts

# TODO: drop mesozooplankton and EUP as we are taking those from NPZ

library(tidyverse)
library(viridis)

dat <- read.csv("../data/inverts.csv")

# turn to long format
dat1 <- dat %>% pivot_longer(cols = BB:ZM, names_to = "Prey_code", values_to = "Percent_diet")

# group by atlantis code and take average of percentages - need to weight by biomass though
counter <- dat1 %>% 
  select(Pred_code,Biomass..t.km2.) %>% 
  distinct() %>% group_by(Pred_code) %>% 
  summarise(Biom_ag = sum(Biomass..t.km2.)) %>% 
  ungroup()

dat2 <- dat1 %>%
  mutate(Prop_diet = Percent_diet/100) %>%
  left_join(counter, by = "Pred_code") %>%
  mutate(Weighted_prop = Biomass..t.km2.*Prop_diet) %>%
  group_by(Pred_code,Biom_ag,Prey_code) %>%
  summarise(Diet_comp = sum(Weighted_prop)/Biom_ag) %>%
  ungroup() %>%
  distinct()

# join with longer names
ag <- read.csv('../data/GOA_Groups.csv', fileEncoding = 'UTF-8-BOM')
ag <- ag %>% select(Code, Name) %>% set_names('Code', 'Name')

dat3 <- dat2 %>% 
  left_join(ag, by = c('Pred_code'='Code')) %>% 
  left_join(ag, by = c('Prey_code'='Code')) %>%
  select(Name.x,Name.y,Diet_comp) %>%
  set_names(c("Predator","Prey","Proportion")) 

dat3 %>% group_by(Predator) %>% summarise(check = sum(Proportion)) %>% pull(check) %>% mean()

# plot

p <- dat3 %>%
  mutate(Proportion = na_if(Proportion,0)) %>%
  ggplot()+
  geom_tile(aes(x = Prey, y = Predator, fill = Proportion), color = "darkgrey")+
  scale_fill_viridis()+
  #scale_fill_gradient2(low = "blue", high = "red", midpoint = .5, na.value="grey")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1.05, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))+
  labs(title = "Proportion of predator diet", 
       x = "Prey", y = "Predator",
       fill = "Proportion of ingested prey")

p

ggsave("inverts.png",p,width = 12, height = 8)

# Prepare PPREY for inverts -----------------------------------------------

dat4 <- dat2 %>% select(-Biom_ag)

dat5 <- expand.grid(unique(dat4$Pred_code), ag$Code) %>%
  as.data.frame() %>%
  set_names(c('Pred_code','Prey_code')) %>% 
  full_join(dat4, by = c('Pred_code','Prey_code')) %>%
  mutate(Diet_comp = replace_na(Diet_comp, 0)) %>%
  mutate(Diet_comp = Diet_comp/10) %>% # as a starting ballpark, as it was done in CalCurrent to start
  pivot_wider(names_from = 'Prey_code', values_from = 'Diet_comp') %>%
  arrange(factor(Pred_code, levels = ag$Code)) %>%
  filter(Pred_code != 'EUP', Pred_code != 'ZM') %>% # dropping EUP and ZM because we do those from NPZ
  mutate(name = paste0('pPREY', Pred_code, ' ', (length(ag$Code)+3))) %>% # +3 because of the detrital sediment
  select(name, KWT:DR) %>%
  mutate(DCsed = 1e-09,   # add sediment columns
         DLsed = 1e-09,
         DRsed = 1e-09)

# write out
write.csv(dat5, '../output/inverts_pprey.csv', row.names = F)
