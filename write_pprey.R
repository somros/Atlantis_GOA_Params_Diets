# 2/17/2022 Alberto Rovellini

# Diets Part 8: bring it all together and write PPREY

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
           value=replace_na(value,0)) # assuming that NA is 0
  
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

# drop: 
# salmon and herring as we will add them with GOAIERP
# mammals and sharks as we will de them from other sources
# seabirds

to_drop <- c('SCH','SCO','SPI','SCM','SSO','HER',
             'KWT','WHT','KWR','DOL','WHG','WHH','WHB','SSL','PIN',
             'BDF','BDI','BSF','BSI',
             'SHD','SHP')

pprey_verts <- pprey_verts %>%
  filter(predcode %in% setdiff(goa_codes, to_drop)) 

# GOAIERP salmon and herring ----------------------------------------------

pprey_goaierp <- read.csv('../output/goaierp_pprey.csv') %>%
  mutate(predcode = substr(name,7,(nchar(name)-4)))

# Invertebrates -----------------------------------------------------------
pprey_inverts <- read.csv('../output/inverts_pprey.csv')

pprey_inverts <- pprey_inverts %>%
  mutate(predcode = substr(name, 6, 8))

# drop blanks from predcode column
pprey_inverts$predcode <- gsub(' ', '', pprey_inverts$predcode)

# drop EUP from these because we enter them from NPZ
pprey_inverts <- pprey_inverts %>%
  filter(predcode != 'EUP')

# check what inverts we are missing from Ecopath
inverts_that_eat <- goa_groups %>% 
  filter(GroupType %in% setdiff(goa_groups$GroupType, c('FISH','MAMMAL','BIRD','SHARK','PHYTOBEN','LG_PHY','SM_PHY','SED_BACT','PL_BACT','CARRION','LAB_DET','REF_DET'))) %>%
  pull(Code)

setdiff(inverts_that_eat, pprey_inverts$predcode)
# missing "BG"  "BO"  "EUP" "ZM"  "ZS", but we do the last 3 from NPZ 
missing <- c('BG','BO')

# add missing inverts - we will need to fill these manually in PPREY

pprey_missing_inverts <- data.frame('name'=paste0('pPREY', missing, ' ', numgroups+3),
                                    matrix(0, nrow = 1, ncol = numgroups+3)) %>%
  mutate(predcode = substr(name, 6, 7)) %>%
  set_names(colnames(pprey_inverts))

# Plankton from NPZ -------------------------------------------------------
pprey_plankton <- read.csv('../output/plankton_pprey.csv') 

pprey_plankton <- pprey_plankton %>%
  mutate(predcode = substr(name, 6, 8))

# drop blanks from predcode column
pprey_plankton$predcode <- gsub(' ', '', pprey_plankton$predcode)

# Mammals and sharks ------------------------------------------------------
pprey_mammals_and_sharks <- read.csv('../output/mammals_sharks_pprey.csv') 

pprey_mammals_and_sharks <- pprey_mammals_and_sharks %>%
  mutate(predcode = substr(name, 7, 9))

# drop blanks from predcode column
pprey_mammals_and_sharks$predcode <- gsub(' ', '', pprey_mammals_and_sharks$predcode)

# Seabirds ----------------------------------------------------------------
pprey_seabirds <- read.csv('../output/seabirds_pprey.csv') 

pprey_seabirds <- pprey_seabirds %>%
  mutate(predcode = substr(name, 7, 9))

# drop blanks from predcode column
pprey_seabirds$predcode <- gsub(' ', '', pprey_seabirds$predcode)

# Steller sea lions -------------------------------------------------------
pprey_ssl <- read.csv('../output/ssl_pprey.csv') 

pprey_ssl <- pprey_ssl %>%
  mutate(predcode = substr(name, 7, 9))

# drop blanks from predcode column
pprey_ssl$predcode <- gsub(' ', '', pprey_ssl$predcode)

# Bind all -------------------------------------------------------

# now bind everything together, reorder them, and and write them out as a csv
pprey_all <- rbind(pprey_verts, 
                   pprey_goaierp, 
                   pprey_inverts, 
                   pprey_missing_inverts, 
                   pprey_plankton,
                   pprey_mammals_and_sharks,
                   pprey_seabirds,
                   pprey_ssl) %>%
  arrange(factor(predcode, levels = goa_codes))

# # as a cheat, still keep a small value for birds, fish, and large sharks until we address their diets (probably manually)
# # TODO: fix this
# these <- c((goa_groups %>% filter(GroupType %in% c('BIRD','MAMMAL')) %>% pull(Code)), c('SHP','SHD'))
# pprey_all[pprey_all$predcode %in% these,2:82] <- 1e-9
# 
pprey_all <- pprey_all %>% select(-predcode)

# we are missing capelin, sandlace, BG, BO
# From Aydin et al. (2007): "Sand lance diet composition information was unavailable for the GOA, so a similar diet composition of 90% euphausiids and 10% copepods was applied to this group as was assumed for all species of small pelagic forage fish"
# This applies to both. We do the same here except we add large zooplankton (not a group in Ecopath)
# and redistribute, as a first pass, 80% EUP, 10% ZL, 10% ZM
# For benthic grazers assume: 30% PL, 30% MA, 20% BB, 10% DLsed, 10% DCsed
# For meiobenthos assume: 70% BB, 10% DLsed, 10% DCsed, 10% DRsed
#Capelin
pprey_all[grep('CAP',pprey_all$name),]$EUP <- 0.8/10
pprey_all[grep('CAP',pprey_all$name),]$ZL <- 0.1/10
pprey_all[grep('CAP',pprey_all$name),]$ZM <- 0.1/10
#Sandlance
pprey_all[grep('SAN',pprey_all$name),]$EUP <- 0.8/10
pprey_all[grep('SAN',pprey_all$name),]$ZL <- 0.1/10
pprey_all[grep('SAN',pprey_all$name),]$ZM <- 0.1/10
#BG
pprey_all[grep('BG',pprey_all$name),]$PL <- 0.3/10
pprey_all[grep('BG',pprey_all$name),]$MA <- 0.3/10
pprey_all[grep('BG',pprey_all$name),]$BB <- 0.2/10
pprey_all[grep('BG',pprey_all$name),]$DLsed <- 0.1/10
pprey_all[grep('BG',pprey_all$name),]$DCsed <- 0.1/10
#BO
pprey_all[grep('BO',pprey_all$name),]$BB <- 0.7/10
pprey_all[grep('BO',pprey_all$name),]$DLsed <- 0.1/10
pprey_all[grep('BO',pprey_all$name),]$DCsed <- 0.1/10
pprey_all[grep('BO',pprey_all$name),]$DRsed <- 0.1/10

write.csv(pprey_all,'../output/goa_pprey_matrix.csv',row.names = F)
