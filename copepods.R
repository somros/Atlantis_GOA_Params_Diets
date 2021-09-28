# Check how many of the copepod entries are unidentified

library(data.table)
library(tidyverse)

# In the stomach content data, prey are aggregated. To get down to copepod species, let us use the NODC codes and get the names from the key

# first, get key and subset to entries containing the worth "copepod"

load("data/nodc_matching_key.Rdata")

# check that nodc are unique

test <- nodc_matching_key %>% select(Nodc) %>% distinct() %>% nrow()
nrow(nodc_matching_key)-test # threre are nine duplicates. Which ones are these?

dups <- nodc_matching_key[nodc_matching_key %>% select(Nodc) %>% duplicated(),]

chech_dups <- nodc_matching_key %>% filter(Nodc %in% dups$Nodc)

# there are a few duplicated NODC, but overall not a problem.
# now identify in the key everything that refers to a copepod. Keep it generic

copepods <- nodc_matching_key %>% filter(grepl("copepod", Name.x, ignore.case = TRUE))

# now read in the stomach content data from https://apps-afsc.fisheries.noaa.gov/refm/reem/webdietdata/dietdataintro.php

stomach <- read.csv("data/1_Gulf Of Alaska_RawPP.csv")

# subset to the Prey_nodc of the copepods. join with the key for the prey names

prey_is_copepod <- stomach %>% filter(Prey_nodc %in% copepods$Nodc)
prey_is_copepod <- prey_is_copepod %>% select(Prey_nodc) %>% left_join(copepods, by =c("Prey_nodc" = "Nodc")) 

# how many of the copepod entries are actually in the data?
prey_is_copepod %>% select(Name.x) %>% distinct() %>% pull()

# these do not appear
setdiff(copepods %>% select(Name.x) %>% pull(), prey_is_copepod %>% select(Name.x) %>% distinct() %>% pull())

# level of taxonomic ID in data is low. Where are all the records in the key from? Perhaps from the Prey length data set

percentages <- prey_is_copepod %>% group_by(Name.x) %>% tally() %>% ungroup() %>% mutate(perc = n/sum(n)*100)

# plot this

ggplot(data = prey_is_copepod, aes(x = Name.x))+
  geom_bar(aes(y = (..count..)/sum(..count..)))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# in the gut content data, 66% of the entries are "Calanoida", which could be Neocalanus or Pseudocalanus.
# Other most abundant occurrences are "large" and "medium" copepods. 
# With this data it will not be possible to parameterize diets on large and small copepods.

