# comparing calibrated diet, true diet (as pre-calibration PPREY), and diet to calibrate
# Comparing the proportions of the true diets to the proportions in the dietcheck output
# Doing so as prop_new/prop_true

library("ggplot2")
library("knitr")
library("dplyr")
library("tidyr")
library(RColorBrewer)
library("readxl")
library(tidyr)
library(tibble)
library(reshape2)
library(here)

#############################################
# true diet:

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
  #    cbind(PreyAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\   .*","",names_pprey)),1,1)%in%c(1,2),
  #                                 substr(gsub("pPREY", "", gsub("\\   .*","",names_pprey)),1,1),
  #                                 "2")) %>%  
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

# tis produces a long table with the initial PPREY values
# These are not necessarily proportional diets
formatted_DM_original_long <- formatted_DM_original %>%
  pivot_longer(5:ncol(formatted_DM_original), values_to="diet_init", names_to="Prey") %>%
  group_by(label, Predator, PredatorAgeClass, Prey) %>%
  mutate(diet_init= as.numeric(diet_init)) %>%
  summarize(diet_init=sum(diet_init, na.rm=T)) # this does not seem to do anything

# These should add up to 0.1 (or close) for each "row"
# formatted_DM_original_long %>%
#   group_by(label) %>%
#   summarize(check = sum(diet_init)) %>%
#   pull(check) %>%
#   hist()

# Add prey stage
formatted_DM_original_long <- formatted_DM_original_long %>%
  mutate(PreyAgeClass = as.numeric(substr(label, 6,6))) %>%
  mutate(PreyAgeClass = replace_na(PreyAgeClass, 2)) %>%
  mutate(PreyAgeClass = ifelse(PreyAgeClass == 1, "Juvenile","Adult"))

# Clean, aggregate without life stages and turn into proportions
true_diet <- formatted_DM_original_long %>%
  group_by(Predator, PredatorAgeClass, Prey) %>%
  summarize(diet_init = sum(diet_init, na.rm = T)) %>% # add up pprey values on juv and adult prey
  group_by(Predator, PredatorAgeClass) %>%
  mutate(tot = sum(diet_init, na.rm = T)) %>%
  mutate(diet_prop_true = diet_init / tot) %>%
  ungroup() %>%
  dplyr::select(Predator, PredatorAgeClass, Prey, diet_prop_true)

# # check that they add up to 1
# true_diet %>% group_by(Predator, PredatorAgeClass) %>% summarize(check = sum(diet_prop_true)) %>% pull(check)
# true_diet %>% filter(Predator == "BDF", Prey == "ATF")

###################################################################################
# run that we modified - how bad was it before?
pred_diets <- read.table("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1401/outputGOA01401_testDietCheck.txt", header=T)
# read in age at maturity - mind that this is not likely to align with the age at which the ontogenetic shift in diets occurs
age_mat <- read.table("../data/age_mat.txt", sep = ' ') %>%
  select(V1, V4) %>%
  purrr::set_names(c('Code','age')) %>%
  mutate(Code = gsub('_age_mat','',Code))

formatted_DM_pred_long <- pred_diets %>%
  filter(Time==365) %>% # can consider using later time steps
  pivot_longer(6:ncol(pred_diets), values_to="diet_prop", names_to="Prey") %>%
  left_join(age_mat, by = c('Predator'='Code')) %>%
  mutate(PredatorAgeClass=ifelse(!Predator %in% no_age,
                                 ifelse(Cohort < age, "Juvenile", "Adult"), "Adult")) %>% 
  group_by(Time, Predator, PredatorAgeClass, Prey) %>% 
  summarize(diet_prop_est=sum(diet_prop, na.rm=T)) %>% # summing contributions to cohorts into life stages
  group_by(Time, Predator, PredatorAgeClass) %>%
  mutate(diet_prop_tot=sum(diet_prop_est, na.rm=T)) %>% 
  mutate(diet_prop_est=diet_prop_est/diet_prop_tot) %>% # rescales the contributions by the total, to make sure there are proportional (do we want this though?)
  ungroup() %>%
  dplyr::select(Predator, PredatorAgeClass, Prey, diet_prop_est) %>%
  drop_na() # CAREFUL: dropping all non-consumers from PPREY

# check
# make sure you subset to the species of interest
# some examples:
top_preds <- FG_tab_Atl %>%
  filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>%
  pull(Code)
# tier3 <- c("ATF","COD","POL","ATF","POP","SBF","FFS","FFD","RFS","HAL","REX","FHS")
# forage <- c("HER","EUL","CAP","SAN","FOL")

todo <- top_preds

check_old <- true_diet %>%
  filter(Predator %in% todo) %>%
  left_join(formatted_DM_pred_long %>%
              filter(Predator %in% todo),
            by = c("Predator", "PredatorAgeClass", "Prey")) %>%
  filter(Predator %in% todo) %>%
  filter(diet_prop_true > 0.01) %>%
  filter(diet_prop_est > 0.01) %>%
  group_by(Predator, PredatorAgeClass) %>%
  mutate(check = diet_prop_est / diet_prop_true) %>%
  ungroup() %>%
  distinct()

hist(check_old$check, breaks = 80)

# new run to check

# model output (the latest that we want to fix)
pred_diets <- read.table("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1408/outputGOA01408_testDietCheck.txt", header=T)

formatted_DM_pred_long <- pred_diets %>%
  filter(Time==365) %>% # can consider using later time steps
  pivot_longer(6:ncol(pred_diets), values_to="diet_prop", names_to="Prey") %>%
  left_join(age_mat, by = c('Predator'='Code')) %>%
  mutate(PredatorAgeClass=ifelse(!Predator %in% no_age,
                                 ifelse(Cohort < age, "Juvenile", "Adult"), "Adult")) %>% 
  group_by(Time, Predator, PredatorAgeClass, Prey) %>% 
  summarize(diet_prop_est=sum(diet_prop, na.rm=T)) %>% # summing contributions to cohorts into life stages
  group_by(Time, Predator, PredatorAgeClass) %>%
  mutate(diet_prop_tot=sum(diet_prop_est, na.rm=T)) %>% 
  mutate(diet_prop_est=diet_prop_est/diet_prop_tot) %>% # rescales the contributions by the total, to make sure there are proportional (do we want this though?)
  ungroup() %>%
  dplyr::select(Predator, PredatorAgeClass, Prey, diet_prop_est) %>%
  drop_na() # CAREFUL: dropping all non-consumers from PPREY

# check
check_new <- true_diet %>%
  left_join(formatted_DM_pred_long, by = c("Predator", "PredatorAgeClass", "Prey")) %>%
  filter(Predator %in% todo) %>%
  distinct() %>%
  filter(diet_prop_true > 0.01) %>%
  filter(diet_prop_est > 0.01) %>%
  group_by(Predator, PredatorAgeClass) %>%
  mutate(check = diet_prop_est / diet_prop_true) %>%
  ungroup() 

hist(check_new$check, breaks = 80)
