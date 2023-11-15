# looking at detailed diets
# This output is really large (0.3 GB per time step it seems), so produce it sparingly
# Good for checking diet proportions after 1 year
# Good for checking consumption after one year
# Things to look at: diet proportions (as bar chart) from different runs compared to the "true" PPREY
# Consumption of each prey as their eaten / biomass
# Consumption of each predator as consumed / biomass
library(tidyverse)

diet1 <- read.table("../../output_files/data/out_1410/outputGOA01410_testDetailedDietCheck.txt",sep=" ",header=T)
diet2 <- read.table("../../output_files/data/out_1411/outputGOA01411_testDetailedDietCheck.txt",sep=" ",header=T)
diet3 <- read.table("../../output_files/data/out_1412/outputGOA01412_testDetailedDietCheck.txt",sep=" ",header=T)
dietBase <- read.table("../../output_files/data/out_1413/outputGOA01413_testDetailedDietCheck.txt",sep=" ",header=T)

# need age at maturity
agemat <- read.table("../data/age_mat.txt", sep = " ")
agemat <- agemat %>% select(1,4) %>% set_names("Code", "agemat") %>% mutate(Code = str_replace(Code, "_age_mat", ""))

# drop spatial elements, aggregate age classes into 
# agemat likely wrong - check

clean_diet <- function(detailedDiet, run){
  
  dietStage <- detailedDiet %>%
    slice_max(Time) %>% # for longer files may need to cut this
    group_by(Predator, Cohort) %>%
    summarise(across(KWT:DR, sum)) %>%
    pivot_longer(-c(Predator, Cohort), names_to = "Prey", values_to = "mt") %>%
    left_join(agemat, by = c("Predator"="Code")) %>%
    mutate(Stage = ifelse(Cohort < agemat, "Juvenile", "Adult")) %>%
    mutate(Stage = replace_na(Stage, "Adult")) %>%
    select(Predator, Stage, Prey, mt) %>%
    group_by(Predator, Stage, Prey) %>%
    summarize(mt = sum(mt)) %>%
    mutate(Run = run)
  return(dietStage)
  
}

dietStage1 <- clean_diet(diet1, "Top_preds")  
dietStage2 <- clean_diet(diet2, "Groundfish")  
dietStage3 <- clean_diet(diet3, "Forage") 
dietStageBase <- clean_diet(dietBase, "Base")

# remove large files and collect
rm(list = c("diet1", "diet2", "diet3", "dietBase"))
gc()

# bind
dietStage <- rbind(dietStageBase, dietStage1, dietStage2, dietStage3)

# get proportions
dietStage <- dietStage %>%
  group_by(Run, Predator, Stage) %>%
  mutate(tot = sum(mt, na.rm = T)) %>%
  ungroup() %>%
  mutate(prop = mt / tot) %>%
  select(-tot)

# introduce true diets from original PPREY matrix
repAtlmat <- "../../build_init_prm_10COHORTS/data/GOAbioparam_test.prm" #true diets
#setwd("/home/atlantis/calcurrenttestruns/calcurrenttestrunssr04/calcurrenttestruns/AtlantisCalCurrV4_S4_0124/outputFolder")
FG_tab_Atl <- read.csv("../data/GOA_Groups.csv")

no_age <- c(FG_tab_Atl$GroupType[FG_tab_Atl$NumStages==1], "DCsed", "DLsed", "DRsed")
bio_prm <- readLines(repAtlmat)

# identify and index the PPREY matrix from the PRM file
diets_start <- grep("pPREY1KWT1", bio_prm) # flag of the first line in the PRM - change to first species
pprey_ind <- which(startsWith(x=bio_prm, "pPREY")==T)
diets_end <- max(pprey_ind)+2
DM_prm <- bio_prm[seq(diets_start,diets_end)]
names_pprey <- bio_prm[pprey_ind]
val_pprey <- bio_prm[pprey_ind+1]

FG <- gsub(" ","",unique(gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\   .*","",names_pprey))))) # consumers
names_pprey_age <- unique(gsub("pPREY","",gsub('[[:digit:]]+', '',  gsub("\\   .*","", grep("pPREY1", names_pprey, value=T))))) # age-structured consumers

DM_to_format <- t(
  sapply(seq(1, length(val_pprey)),
         function(x){
           vec <- unlist(strsplit(val_pprey[x], " "))
           return(vec)
         })
)

for_order_FG <- FG_tab_Atl[,1] # groups from group file
colnames(DM_to_format) <- c(for_order_FG,c("DCsed","DLsed","DRsed"))

formatted_DM_original <- DM_to_format %>%
  as_tibble() %>%
  cbind(label=gsub("\\ .*","",names_pprey)) %>%
  cbind(Predator=gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\ .*","",names_pprey)))) %>%
  cbind(PreyAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1)%in%c(1,2),
                            substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1),
                            "1")) %>%  
  cbind(PredatorAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))))%in%c(1,2), 
                                substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey)))),
                                "2")
  ) %>% 
  mutate(PredatorAgeClass=ifelse(PredatorAgeClass==1,"Juvenile", "Adult"),
         PreyAgeClass=ifelse(PreyAgeClass==2,"Adult", "Juvenile")) %>%
  dplyr::select(c("label","Predator", "PreyAgeClass", "PredatorAgeClass", colnames(DM_to_format)))

# this produces a long table with the initial PPREY values
true_diets <- formatted_DM_original %>%
  pivot_longer(5:ncol(formatted_DM_original), values_to="diet_init", names_to="Prey") %>%
  mutate(diet_init= as.numeric(diet_init)) %>%
  group_by(Predator, PredatorAgeClass, Prey) %>%
  summarize(diet_init = sum(diet_init)) %>% # drop preystage
  rename(Stage = PredatorAgeClass) %>%
  group_by(Predator, Stage) %>%
  mutate(tot = sum(diet_init)) %>%
  ungroup() %>%
  mutate(prop = diet_init / tot)

# Proportional diets comparisons ------------------------------------------

namekey <- FG_tab_Atl %>% dplyr::select(Code, LongName)

propdiet <- dietStage %>%
  dplyr::select(Predator, Stage, Prey, prop, Run) %>%
  rbind(true_diets %>% dplyr::select(Predator, Stage, Prey, prop) %>% mutate(Run = "True")) %>%
  left_join(namekey, by = c("Predator"="Code")) %>%
  rename(PredatorName = LongName) %>%
  left_join(namekey, by = c("Prey"="Code")) %>%
  rename(PreyName = LongName)

# t1 <- propdiet %>% filter(Predator == "SHD")

# list groups
top_preds <- FG_tab_Atl %>%
  filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>%
  pull(Code)
tier3 <- c("ATF","COD","POL","ATF","POP","SBF","FFS","FFD","RFS","HAL","REX","FHS")
forage <- c("HER","EUL","CAP","SAN","FOL") # FOL is wrong change to FOS

fg_to_plot <- c(top_preds, tier3, forage)

# plot
for(i in 1:length(fg_to_plot)){
  
  this_pred <- fg_to_plot[i]
  
  if(this_pred %in% top_preds) {
    this_run <- "Top_preds"
  } else if (this_pred %in% tier3) {
    this_run <- "Groundfish"
  } else {
    this_run <- "Forage"
  }
  
  plot_prop <- propdiet %>%
    mutate(Run = factor(Run, levels = c("Base","True","Top_preds","Groundfish","Forage"))) %>%
    filter(Run %in% c("Base","True",this_run)) %>%
    filter(Predator == this_pred)%>%
    filter(prop > 0.01) %>% #filter out traces
    ggplot()+
    geom_bar(aes(x = Stage, y = prop, fill = PreyName), stat = "identity", position = "stack")+
    theme_bw()+
    ggtitle(paste("Predator:", propdiet %>% filter(Predator == this_pred) %>% pull(PredatorName) %>% unique()))+
    facet_wrap(~Run)
  plot_prop
  # save
  ggsave(paste0("figures/",this_pred,".png"), plot_prop, width = 8, height = 6)
  
}
  

# Check consumption -------------------------------------------------------

# read in total biomass
# these will be different for each run (Base and the scenarios)
biomBase <- read.table("../../output_files/data/out_1413/outputGOA01413_testAgeBiomIndx.txt",sep=" ",header=T)
biom1 <- read.table("../../output_files/data/out_1410/outputGOA01410_testAgeBiomIndx.txt",sep=" ",header=T)
biom2 <- read.table("../../output_files/data/out_1411/outputGOA01411_testAgeBiomIndx.txt",sep=" ",header=T)
biom3 <- read.table("../../output_files/data/out_1412/outputGOA01412_testAgeBiomIndx.txt",sep=" ",header=T)

clean_biom <- function(biomage, run){
  
  this_biomage <- biomage %>% 
    slice_max(Time) %>% 
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(agemat, by = "Code") %>%
    mutate(Stage = ifelse(Age < agemat, "Juvenile", "Adult")) %>%
    mutate(Stage = replace_na(Stage, "Adult")) %>%
    group_by(Code, Stage) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    mutate(Run = run)
  return(this_biomage)
}

biomStage1 <- clean_biom(biom1, "Top_preds")  
biomStage2 <- clean_biom(biom2, "Groundfish")  
biomStage3 <- clean_biom(biom3, "Forage") 
biomStageBase <- clean_biom(biomBase, "Base")

# bind
biomStage <- rbind(biomStageBase, biomStage1, biomStage2, biomStage3)

# predator consumption
pred_cons <- dietStage %>%
  group_by(Run, Predator, Stage) %>%
  summarize(eaten_mt = sum(mt)) %>%
  left_join(biomStage, by = c("Run", "Predator"="Code", "Stage")) %>%
  mutate(annual_consumption = eaten_mt / biomass_mt)

# No. These numbers are really small. This highlights one of a few things:
# Consumption is way too low (they do not eat)
# The units in the detailedDietCheck file are not mt
# The output is unreliable
# I think consumption is indeed low in this model, but it should not be negligible like these numbers suggest.

