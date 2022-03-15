# 2/17/2022 Alberto Rovellini
# This code takes output from the diet code and writes out a template for the PPREY matrix for Atlantis GOA

library(tidyverse)
library(data.table)

# Read data ---------------------------------------------------------------

goa_groups <- read.csv('../data/GOA_Groups.csv', fileEncoding = 'UTF-8-BOM')
goa_key <- goa_groups %>% select(Code,Name)
goa_codes <- goa_groups %>% pull(Code) # fix DF
numgroups <- length(goa_codes)

# data from REEM bottom trawl groundfish data
reem_files <- list.files('../output/REEM/',pattern='.csv',full.names = T)
reem_list <- list()

# REEM --------------------------------------------------------------------
for (i in 1:length(reem_files)){
  this_reem <- read.csv(reem_files[i])
  
  # establish what combination of prey-pred the file contains
  this_name <- list.files('../output/REEM/',pattern='.csv',full.names = F)[i]
  
  this_pprey <- case_when( 
    grepl('pred_J_prey_J',this_name) ~ 'pPREY1XXX1',
    grepl('pred_A_prey_J',this_name) ~ 'pPREY1XXX2',
    grepl('pred_J_prey_A',this_name) ~ 'pPREY2XXX1',
    grepl('pred_A_prey_A',this_name) ~ 'pPREY2XXX2'
  )
  
  # change names to 3-letter codes because we need those in the pPREY matrix
  this_reem <- this_reem %>%
    left_join(goa_key, by = c('Pred_at_name'='Name')) %>%
    left_join(goa_key, by = c('Prey_at_name'='Name')) %>%
    select(Code.x,Code.y,pcW_mean) %>%
    set_names(c('pred','prey','value'))
  
  # complete the combinations to add all functional groups
  all_codes <- expand.grid(goa_codes,c(goa_codes,"DCsed","DLsed","DRsed")) %>%
    set_names(c('pred','prey'))
  
  all_codes <- all_codes[order(match(all_codes$pred,goa_codes)),] # order them like the groups in Group.csv
  
  all_pprey <- all_codes %>%
    left_join(this_reem, by=c('pred','prey')) %>%
    mutate(value=value/10, # first approximation to get pPREY from diet proportion preferences
           value=replace_na(value,1e-9)) # using a small number as fillvalue for now - TODO change this
  
  all_pprey_wide <- all_pprey %>% pivot_wider(names_from = prey, values_from = value)
  
  all_pprey_wide <- all_pprey_wide %>% mutate(pred=paste0(gsub('XXX.*','',this_pprey),
                                                          pred,
                                                          gsub('.*XXX','',this_pprey),
                                                          ' ',
                                                          numgroups+3)) %>% #+3 because of the detrital sediment
    set_names(c('name',colnames(all_pprey_wide)[-1]))
  
  reem_list[[i]] <- all_pprey_wide
}

pprey <- rbindlist(reem_list)

# sort it in the correct order for readability
pprey <- pprey %>% mutate(predcode = substr(name,7,(nchar(name)-4)),
                          preystage = substr(name,6,6),
                          predstage = substr(name,(nchar(name)-3),(nchar(name)-3))) %>%
  arrange(factor(predcode,levels = goa_codes),preystage,predstage) %>%
  select(-predstage,-preystage)

# remove all inverts, plankton, detritus, bacteria
inverts <- goa_groups %>% 
  filter(GroupType %in% setdiff(goa_groups$GroupType, c('FISH','MAMMAL','BIRD','SHARK'))) %>%
  pull(Code)

pprey_verts <- pprey %>% filter(predcode %in% setdiff(goa_codes,inverts))

# Invertebrates -----------------------------------------------------------
pprey_inverts <- read.csv('../output/inverts_pprey.csv')

pprey_inverts <- pprey_inverts %>%
  mutate(predcode = substr(name, 6, 8))

# drop blanks from predcode column
pprey_inverts$predcode <- gsub(' ', '', pprey_inverts$predcode)

# check what inverts we are missing from Ecopath
inverts_that_eat <- goa_groups %>% 
  filter(GroupType %in% setdiff(goa_groups$GroupType, c('FISH','MAMMAL','BIRD','SHARK','PHYTOBEN','LG_PHY','SM_PHY','SED_BACT','PL_BACT','CARRION','LAB_DET','REF_DET'))) %>%
  pull(Code)

missing <- setdiff(inverts_that_eat, pprey_inverts$predcode) # missing "BG" "BO" "ZM" "ZS"

# add missing inverts

pprey_missing_inverts <- data.frame('name'=paste0('pPREY', missing, ' ', numgroups+3),
                                    matrix(1e-09, nrow = 1, ncol = numgroups+3)) %>%
  mutate(predcode = substr(name, 6, 7)) %>%
  set_names(colnames(pprey_inverts))

# Bind all -------------------------------------------------------

# now bind verts and inverts and write them out as a csv
pprey_all <- rbind(pprey_verts, pprey_inverts, pprey_missing_inverts) %>%
  select(-predcode)

write.csv(pprey_all,'../output/goa_pprey_matrix.csv',row.names = F)
