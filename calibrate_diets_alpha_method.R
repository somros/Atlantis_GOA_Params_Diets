#########################################
###### Calibrating the diet matrix - method 2
#########################################
library(tidyverse)

#############################################
# Alberto Rovellini
## Starting from version with amplified consumption, run 1689
## Here we use as target (i.e., "true" diets) the original, pre-calibration PRM file

source("diet_calibration_functions_V2.R")

# which runs are we working on?
runs_folder <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/"
run_to_calibrate <- 1689

my_alpha <- 0.75

# which groups?
FG_tab_Atl <- read.csv("../data/GOA_Groups.csv")
for_order_FG <- FG_tab_Atl[,1] # groups from group file

no_age <- c(FG_tab_Atl$Code[FG_tab_Atl$NumStages==1], "DCsed", "DLsed", "DRsed")
top_preds <- FG_tab_Atl %>%
  filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>%
  pull(Code)
tier3 <- c("ATF","COD","POL","ATF","POP","SBF","FFS","FFD","RFS","HAL","REX","FHS")
forage <- c("HER","EUL","CAP","SAN","FOL")

todo <- tier3 # pick some groups, or do whole food web at once

# Read target diets -------------------------------------------------------
# For GOA, this is a pre-calibration PPREY matrix, but it could be other form of data
repAtlmat <- "../../build_init_prm_10COHORTS/data/GOAbioparam_test.prm"
bio_prm <- readLines(repAtlmat)

# Get original diet matrix data
formatted_DM_original_long <- format_diet_matrix(
  prm_file = repAtlmat,
  for_order_FG = for_order_FG,
  scale_factor = 10
)

# Diet output from run to calibrate ---------------------------------------
## Starting from version with amplified consumption, run 1689

# model output (the latest that we want to fix)
pred_diets <- read.table(paste0(runs_folder,"out_",run_to_calibrate,"/outputGOA0",run_to_calibrate,"_testDietCheck.txt"), header=T)
# read in age at maturity - this is what the model uses for the ontogenetic shift in diets
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
  group_by(Time, Predator, Prey, PredatorAgeClass) %>% 
  summarize(diet_prop_est=sum(diet_prop, na.rm=T)) %>% # summing contributions to cohorts into life stages
  group_by(Time, Predator, PredatorAgeClass) %>%
  mutate(diet_prop_tot=sum(diet_prop_est, na.rm=T)) %>% 
  mutate(diet_prop_est=diet_prop_est/diet_prop_tot) %>% # rescales the contributions by the total, to make sure there are proportional (do we want this though?)
  ungroup() %>%
  dplyr::select(Predator,PredatorAgeClass,Prey,diet_prop_est)

colnames(formatted_DM_pred_long) <- c("Predator","PredStage","Prey","Value")

# IMPORTANT: this last step above rescales the PPREY contributions so that they sum up to 1.
# In the GOA model we needed to increase individual PPREY values, and the total did not sum up to 1
# Rescaling PPREY entries this way may help with proportional diet contributions, but if may also deflate consumption
# The PPREY matrix is defined as an availability matrix, and there is no real requirement that the values entered are proportional

# check props here:
# tt <- formatted_DM_pred_long %>%
#   group_by(Predator, PredatorAgeClass) %>%
#   summarise(check  = sum(diet_prop_est))


# PPREY matrix to modify --------------------------------------------------
# read in new diet matrix (make this a function as it is the same workflow used above)
latestAtlmat <- paste0(runs_folder,"out_",run_to_calibrate,"/GOAbioparam_test.prm")
bio_prm_latest <- readLines(latestAtlmat)

# Read PPREY to modify, keep prey stages
formatted_DM_latest_long_age <- format_diet_matrix(
  prm_file = latestAtlmat,
  for_order_FG = for_order_FG,
  keep_age_specific = TRUE
)

# Read PPREY to modify, group across prey stages
formatted_DM_latest_long <- format_diet_matrix(
  prm_file = latestAtlmat,
  for_order_FG = for_order_FG
)

# store labels
# For the latest diet matrix, keeping age-specific data and getting names
names_pprey <- format_diet_matrix(
  prm_file = latestAtlmat,
  for_order_FG = for_order_FG,
  return_names_only = TRUE
)

# Apply calibration function ----------------------------------------------
calibrated_diet_result <- calibrate_diet_matrix_grouped(
  target_diet = formatted_DM_original_long,
  input_diet = formatted_DM_latest_long,
  output_diet = formatted_DM_pred_long,
  selected_predators = tier3,
  alpha = my_alpha
)

# isolate the multipliers from the results
these_mult <- calibrated_diet_result[[1]]

# bring pprey to calibrate back in
corr_diets <- formatted_DM_latest_long_age %>%
  left_join(these_mult, by = c("PredatorCODE"="Predator", "PredatorAgeClass"="PredStage","Prey"="Prey")) %>%
  mutate(pprey_new = diet_value * Value,
         check = pprey_new / diet_value) %>%
  mutate(pprey_new = ifelse(pprey_new>=1,0.99,pprey_new)) # make sure they are <1

# turn into a table to be transformed into the pprey matrix
newpprey <- corr_diets  %>%
  mutate(pprey_new = replace_na(pprey_new, 0)) %>%
  dplyr::select(label, Prey, pprey_new) %>%
  pivot_wider(values_from=pprey_new, names_from=Prey) %>%
  dplyr::select(label,all_of(for_order_FG),DCsed,DLsed,DRsed)

# pull the modified groups from here
pprey_file <- paste0('calibrated_tier3_',run_to_calibrate,my_alpha,'.txt')
file.create(pprey_file)

for(i in 1:length(names_pprey)){
  
  this_pprey_name <-  gsub('(81).*','\\1',names_pprey[i])
  this_label <- gsub('\\  81.*','',names_pprey[i])
  
  if(this_label %in% unique(newpprey$label)){
    
    cat(this_pprey_name, file=pprey_file, append=TRUE,'\n')
    cat(unlist(newpprey %>% filter(label == this_label) %>% select(-label)), file=pprey_file, append=TRUE, '\n')
    
  } else {
    
    cat(this_pprey_name, file=pprey_file, append=TRUE,'\n')
    cat(val_pprey[i], file=pprey_file, append=TRUE, '\n')
    
  }
  
}
