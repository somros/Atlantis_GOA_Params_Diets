# Code to explore Atlantis diet interactions from biol.prm and init.nc
# Conceptually similar to Javier's ReactiveAtlantis code, but more detail
# For each predator, it takes PPREY matrix and identifies the "favorite" prey items as prescribed in the preferences
# For each prey items (by ontogenetic stage of the predator), it takes:
# horizontal dists, vertical dists, gape size from biol.prm, weight at age from init.nc
# The function returns: by stage: horizontal overlap (map and coefficient), vertical overlap; by age class: gape size limitations

library(tidyverse)
library(rbgm)
library(maps)
library(mapdata)
library(tidync)
library(sf)
library(viridis)
library(ggh4x)
library(RColorBrewer)
library(patchwork)
library(kableExtra)

# Read data ---------------------------------------------------------------

# geometry
bgm <- read_bgm("../data/GOA_WGS84_V4_final.bgm")
goa_sf <- box_sf(bgm)
st_crs(goa_sf) <- st_crs(attr(goa_sf$geometry, "crs")$proj) 

# functional groups
grps <- read.csv("../data/GOA_Groups.csv")
fg <- grps %>% pull(Code) # groups from group file
verts <- grps %>% filter(GroupType %in% c("MAMMAL","BIRD","SHARK","FISH")) %>% pull(Code)
inverts <- setdiff(fg, verts)

# utility group tables
vertebrate_groups <- grps %>% filter(GroupType%in%c("FISH","SHARK","BIRD","MAMMAL")) %>% mutate(BiomassType="vertebrate")
plankton_groups <- grps %>% filter(GroupType %in% c("PWN",'CEP','LG_ZOO','MED_ZOO','SM_ZOO','LG_PHY','SM_PHY')) %>%
  mutate(BiomassType="plankton")
bottom_groups <- grps %>% filter(GroupType %in% c("MOB_EP_OTHER",'SED_EP_FF','SED_EP_OTHER','PHYTOBEN')) %>%
  mutate(BiomassType="2D")
other_groups <- grps %>% filter(GroupType %in% c("LG_INF","MICROPHTYBENTHOS","SED_BACT","PL_BACT","SM_INF","CARRION","LAB_DET","REF_DET"))%>%
  mutate(BiomassType="other")
biomass_groups <- bind_rows(vertebrate_groups,plankton_groups,bottom_groups,other_groups)

# add to grps df
grps <- grps %>% left_join(biomass_groups)

# age at maturity info
agemat <- read.table("../data/age_mat.txt")
agemat <- agemat %>% mutate(species = substr(V1, 1, 3)) %>% rename(agemat = V2) %>% dplyr::select(species, agemat)

# prm of run to look at
bio_prm <- readLines("../../output_files/data/out_1328/GOAbioparam_test.prm")
# init.nc of run to look at
init_file <- "../../../AtlantisGOA_Base/GOA_cb_summer.nc"
init <- tidync(init_file)
init_nc <- ncdf4::nc_open(init_file)

# Turn PRM and INIT to data frames  ---------------------------------------
# turn pprey matrix to a data frame to work with
no_age <- c(grps$GroupType[grps$NumStages==1], "DCsed", "DLsed", "DRsed")

# identify and index the PPREY matrix from the PRM file
diets_start <- grep("pPREY1KWT1", bio_prm) # flag of the first line in the PRM - change to first species
pprey_ind <- which(startsWith(x=bio_prm, "pPREY")==T)
diets_end <- max(pprey_ind)+2
DM_prm <- bio_prm[seq(diets_start,diets_end)]
names_pprey <- bio_prm[pprey_ind]
val_pprey <- bio_prm[pprey_ind+1]

DM_to_format <- t(
  sapply(seq(1, length(val_pprey)),
         function(x){
           vec <- unlist(strsplit(val_pprey[x], " "))
           return(vec)
         })
)

colnames(DM_to_format) <- c(fg,c("DCsed","DLsed","DRsed"))

pprey_mat <- DM_to_format %>%
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
  mutate(PredatorAgeClass=ifelse(PredatorAgeClass==1,"juvenile", "adult"),
         PreyAgeClass=ifelse(PreyAgeClass==2,"adult", "juvenile")) %>%
  dplyr::select(c("label","Predator", "PreyAgeClass", "PredatorAgeClass", colnames(DM_to_format)))

# create dataframe of horizontal distributions
s <- data.frame()
for(i in 1:length(fg)){
  
  print(paste("doing s for", fg[i]))
  
  this_fg <- fg[i]
  
  # find the handles for the S1-4 parameter lines
  if(this_fg %in% verts) {
    s_pars <- c(
      paste0(paste0("F", this_fg, "_S"),1:4," 109"),
      paste0(paste0("F", this_fg, "_S"),1:4,"juv 109")
    )
  } else {
    s_pars <- paste0(paste0("F", this_fg, "_S"),1:4," 109")
  }
  
  s_lines <- sapply(s_pars, function(s_par) grep(s_par, bio_prm)) 
  if(!is.numeric(s_lines)) {next}
  s_vec_lines <- s_lines + 2 # +2 because the s vector is 2 lines below in GOA prm, in other models may be different
  
  # read names and values
  s_df <- data.frame()
  for(j in 1:length(s_lines)){
    s_mat <- matrix(as.numeric(unlist(str_split(bio_prm[s_vec_lines[j]], " "))), nrow = 1)
    s_df_tmp <- data.frame(s_mat) %>% `colnames<-`(0:108)
    s_df_tmp <- cbind("name" = names(s_lines[j]), s_df_tmp)
    s_df <- rbind(s_df, s_df_tmp)
  }
  
  # add columns for stage and season
  s_df <- s_df %>%
    mutate(species = this_fg,
           stage = ifelse(grepl("juv", name), 
                          "juvenile", 
                          "adult"),
           season = as.numeric(ifelse(grepl("juv", name), 
                                      substr(name, (nchar(name)-7), (nchar(name)-7)),
                                      substr(name, (nchar(name)-4), (nchar(name)-4))))) %>%
    dplyr::select(-name)
  # pivot
  s_df_long <- s_df %>%
    pivot_longer(cols = -c(species, stage, season), names_to = "b", values_to = "s") %>%
    mutate(b = as.numeric(b))
  
  # add to df
  s <- rbind(s, s_df_long)
}

# create dataframe for vertical distributions
# there are day and night distributions for juveniles and adults, so 4 vectors each for verts, 2 each for inverts
v <- data.frame()
for(i in 1:length(fg)){
  
  print(paste("doing v for", fg[i]))
  
  this_fg <- fg[i]
  
  # find the handles for the VERTday and VERTnight parameter lines
  if(this_fg %in% verts) {
    v_pars <- c(
      paste0(paste0("VERTday_", this_fg),1:2,"  6"),
      paste0(paste0("VERTnight_", this_fg),1:2,"  6")
    )
  } else {
    v_pars <- c(
      paste0("VERTday_", this_fg,"  6"),
      paste0("VERTnight_", this_fg,"  6")
    )
  }
  
  v_lines <- sapply(v_pars, function(v_par) grep(v_par, bio_prm)) 
  if(!is.numeric(v_lines)) {next}
  v_vec_lines <- v_lines + 2 # +2 because the v vector is 2 lines below in GOA prm, in other models may be different
  
  # read names and values
  v_df <- data.frame()
  for(j in 1:length(v_lines)){
    v_mat <- matrix(as.numeric(unlist(str_split(bio_prm[v_vec_lines[j]], " "))), nrow = 1)
    v_df_tmp <- data.frame(v_mat) %>% `colnames<-`(6:1) # 1 is the surface here
    v_df_tmp <- cbind("name" = names(v_lines[j]), v_df_tmp)
    v_df <- rbind(v_df, v_df_tmp)
  }
  
  # add columns for stage and day/night
  v_df <- v_df %>%
    mutate(species = this_fg,
           stage = ifelse(grepl("1", name), 
                          "juvenile", 
                          "adult"),
           time = ifelse(grepl("day", name), 
                         "day", 
                         "night")) %>%
    dplyr::select(-name)
  # pivot
  v_df_long <- v_df %>%
    pivot_longer(cols = -c(species, stage, time), names_to = "z", values_to = "v") %>%
    mutate(z = as.numeric(z))
  
  # add to df
  v <- rbind(v, v_df_long)
  
}

# create a dataframe with gape size information
g <- data.frame()
for(i in 1:length(fg)){
  
  print(paste("doing g for", fg[i]))
  
  this_fg <- fg[i]
  
  # find the handles for the KLP_ and KUP parameter lines
  
  g_pars <- c(
    paste0("KLP_", this_fg),
    paste0("KUP_", this_fg)
  )
  
  g_lines <- sapply(g_pars, function(g_par) grep(g_par, bio_prm)) 
  if(!is.numeric(g_lines)) {next}
  
  # read names and values
  g_df <- data.frame()
  for(j in 1:length(g_lines)){
    this_g_line <- unlist(str_split(bio_prm[g_lines[j]], "   "))
    this_g_name <- this_g_line[1]
    this_g_val <- as.numeric(this_g_line[2])
    g_df_tmp <- data.frame("name" = this_g_name, "g" = this_g_val)
    g_df <- rbind(g_df, g_df_tmp)
  }
  
  # add columns for species and whether it is upper or lower limit
  g_df <- g_df %>%
    mutate(species = this_fg,
           limit = ifelse(grepl("KLP", name), 
                          "low", 
                          "high")) %>%
    dplyr::select(-name)
  
  # add to df
  g <- rbind(g, g_df)
  
}

# This function takes a predator and its favorite 5 prey by ontogenetic stage (per PPREY matrix) and makes some plots:
# Horizontal overlap (as the lowest value between s of predator and prey) - how do we handle stages and seasons?
# vertical overlap
# gape size overlap by age class (take size of pred, size of prey )
# the whole concept of favorite prey items needs some thought
# 

diet_from_init <- function(predator) {
  
  print(paste("Extracting info for", predator))
  
  predator_name <- grps %>% filter(Code == predator) %>% pull(Name) # need this for the netcdf variables
  
  # pull pprey for this predator
  this_pprey <- pprey_mat %>% 
    filter(Predator == predator) %>% 
    dplyr::select(-label) %>%
    pivot_longer(cols = -c(Predator, PredatorAgeClass, PreyAgeClass), names_to = "Prey", values_to = "pprey") %>%
    mutate(pprey = as.numeric(pprey)) %>%
    filter(Prey %in% verts)
  
  # make a plot of the predation according to pprey
  d_plot <- this_pprey %>%
    group_by(Predator, PredatorAgeClass, PreyAgeClass) %>%
    mutate(tot = sum(pprey),
           prop = pprey / tot) %>%
    filter(prop > 0) %>%
    arrange(Predator, PredatorAgeClass, PreyAgeClass, desc(prop)) %>%
    mutate(cumprop = cumsum(prop)) %>%
    ungroup() %>%
    filter(cumprop < .9) %>%
    ggplot()+
    geom_bar(aes(x = PreyAgeClass, y = prop, fill = Prey), stat = "identity")+
    theme_bw()+
    labs(x="",y="prop from PPREY",title = paste0(predator, "'s PPREY preferences (proportional)"))+
    facet_grid(~PredatorAgeClass)
  
  d_file <- paste("diets_from_init", predator, "d.png", sep="/")

  ##########################################################################################
  # HORIZONTAL OVERLAP
  # pull horizontal distributions
  # of the predator
  s_pred <- s %>% filter(species == predator)
  
  # need to do this starting from predator life stage
  pred_stages <- s_pred %>% pull(stage) %>% unique()
  s_overlap <- data.frame()
  for(i in 1:length(pred_stages)){
    
    fav_prey_df <- this_pprey %>%
      filter(PredatorAgeClass == pred_stages[i]) %>%
      group_by(Predator, PredatorAgeClass, PreyAgeClass) %>%
      arrange(desc(pprey)) %>%
      filter(pprey > 0) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    fav_prey <- fav_prey_df %>% pull(Prey) %>% unique()
    
    if(length(fav_prey)==0) {next}
    
    # pull dists of the favorite prey for this life stage of the predator
    s_prey <- s %>% filter(species %in% fav_prey)
    
    for(j in 1:length(fav_prey)){ # prey items eaten by this life stage of this predator
      
      prey_stages <- fav_prey_df %>% filter(Prey == fav_prey[j]) %>% pull(PreyAgeClass) %>% unique()
      
      for(k in 1:length(prey_stages)){ # life stages of this prey for the life stage of the predator
        
        seasons <- s_pred %>% filter(stage == pred_stages[i]) %>% pull(season) %>% unique()
        
        for (l in 1:length(seasons)) {
          
          s_pred_tmp <- s_pred %>% 
            filter(stage == pred_stages[i], season == seasons[l]) %>%
            rename(pred = species, pred_stage = stage, s_pred = s)
          
          s_prey_tmp <- s_prey %>% filter(species == fav_prey[j], stage == prey_stages[k], season == seasons[l]) %>%
            rename(prey = species, prey_stage = stage, s_prey = s)
          
          # join
          s_overlap_tmp <- s_pred_tmp %>%
            left_join(s_prey_tmp, by = c("season","b")) %>%
            rowwise() %>%
            mutate(lo = min(c(s_pred, s_prey))) %>%
            ungroup() %>%
            mutate(overlap = sum(lo)) %>% # a value of 0 means no overlap, 1 means complete overlap (identical dists, expected for cannibalism)
            dplyr::select(season, b, pred, pred_stage, prey, prey_stage, lo, overlap)
          
          # now bind into one dataframe for plotting
          s_overlap <- rbind(s_overlap, s_overlap_tmp)
          
        }
      }
    }
    
  }
  
  # if at this stage s_overlap has no rows it means that this group eats no vertebrate prey and we should break the function
  if(nrow(s_overlap)==0){stop("This group eats no predators")}
  
  # drop na's and adjust the stages for facetting
  s_overlap <- s_overlap %>%
    drop_na(prey) %>%
    mutate(pred_stage = paste("predator", pred_stage, sep = ":"),
           prey_stage = paste("prey", prey_stage, sep = ":"))
  
  # make a plot per prey per season
  p <- list()
  for(i in 1:length(unique(s_overlap$prey))){
    
    this_prey <- unique(s_overlap$prey)[i]
    
    for(j in 1:length(unique(s_overlap$season))){
      
      this_season <- unique(s_overlap$season)[j]
      
      idx <- (i -1) * length(unique(s_overlap$season)) + j # make index for the list of plots
      
      p[[idx]] <- goa_sf %>% 
        dplyr::select(box_id) %>%
        full_join(s_overlap, by = c("box_id" = "b")) %>%
        filter(prey == unique(s_overlap$prey)[i], season == unique(s_overlap$season)[j]) %>%
        ggplot()+
        geom_sf(aes(fill = lo)) +
        scale_fill_viridis() +
        theme_bw() +
        geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=round(overlap,2)))+
        labs(x="",y="",fill="Lower s", title = paste(predator,"eating",this_prey,"in season",this_season))+
        facet_grid(prey_stage~pred_stage)
    }
  }
  
  s_length <- length(p)
  
  s_plot <- wrap_plots(p, ncol = 4) # 4 columns arrange the plots by seasons
  
  s_file <- paste("diets_from_init", predator, "s.png", sep="/")
  
  ##########################################################################################
  # VERTICAL OVERLAP
  # pull vertical distributions
  # of the predator
  v_pred <- v %>% filter(species == predator)
  v_overlap <- data.frame()
  for(i in 1:length(pred_stages)){
    
    fav_prey_df <- this_pprey %>%
      filter(PredatorAgeClass == pred_stages[i]) %>%
      group_by(Predator, PredatorAgeClass, PreyAgeClass) %>%
      arrange(desc(pprey)) %>%
      filter(pprey > 0) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    fav_prey <- fav_prey_df %>% pull(Prey) %>% unique()
    
    if(length(fav_prey)==0) {next}
    
    # pull vertical dists of the favorite prey for this life stage of the predator
    v_prey <- v %>% filter(species %in% fav_prey)
    for(j in 1:length(fav_prey)){ # prey items eaten by this life stage of this predator
      
      prey_stages <- fav_prey_df %>% filter(Prey == fav_prey[j]) %>% pull(PreyAgeClass) %>% unique()
      
      for(k in 1:length(prey_stages)){ # life stages of this prey for the life stage of the predator
        
        v_pred_tmp <- v_pred %>% 
          filter(stage == pred_stages[i]) %>%
          rename(pred = species, pred_stage = stage, v_pred = v)
        
        v_prey_tmp <- v_prey %>% filter(species == fav_prey[j], stage == prey_stages[k]) %>%
          rename(prey = species, prey_stage = stage, v_prey = v)
        
        # join
        v_overlap_tmp <- v_pred_tmp %>%
          left_join(v_prey_tmp, by = c("time","z")) %>%
          rowwise() %>%
          mutate(lo = min(c(v_pred, v_prey))) %>%
          group_by(time) %>%
          mutate(overlap = sum(lo)) %>% # a value of 0 means no overlap, 1 means complete overlap (identical dists, expected for cannibalism)
          ungroup() %>%
          dplyr::select(time, z, pred, pred_stage, prey, prey_stage, lo, overlap)
        
        
        # now bind into one dataframe for plotting
        v_overlap <- rbind(v_overlap, v_overlap_tmp)
        
      }
    }
    
  }
  
  # drop na's and adjust the stages for facetting
  v_overlap <- v_overlap %>%
    drop_na(prey) %>%
    mutate(pred_stage = paste("predator", pred_stage, sep = ":"),
           prey_stage = paste("prey", prey_stage, sep = ":"))
  
  p <- list()
  for(i in 1:length(unique(v_overlap$prey))){
    
    this_prey <- unique(v_overlap$prey)[i]
    
    vp <- v_overlap %>%
      filter(prey == unique(v_overlap$prey)[i]) 
    
    p[[i]] <- vp %>%
      ggplot()+
      geom_tile(aes(x = time, y = factor(-z), fill = lo)) +
      scale_fill_viridis() +
      theme_bw() +
      geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=round(overlap,2)), color = "red")+
      labs(x="",y="z (1 = surface)",fill="Lower VERT", title = paste(predator,"eating",this_prey))+
      facet_grid(prey_stage~pred_stage)
    
  }
  
  v_plot <- wrap_plots(p, ncol = 2) # 4 columns arrange the plots by seasons
  
  v_file <- paste("diets_from_init", predator, "v.png", sep="/")
  
  ##########################################################################################
  # GAPE SIZE OVERLAP
  # pull gape size information of the predator
  g_pred <- g %>% filter(species == predator)
  
  # pull weight-at-age information from initial conditions
  # of the predator
  fg_atts <- grps %>% filter(Code==predator)
  
  if(fg_atts$BiomassType!="vertebrate") stop("weight at age only for vertebrates.")
  
  #Extract from the init.nc file the appropriate reserve N time series variables
  resN_vars <- hyper_vars(init) %>% # all variables in the .nc file active grid
    filter(grepl("_ResN",name)) %>% # filter for reserve N
    filter(grepl(predator_name,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate structural N time series variables
  strucN_vars <- hyper_vars(init) %>% # all variables in the .nc file active grid
    filter(grepl("_StructN",name)) %>% # filter for structural N
    filter(grepl(predator_name,name)) # filter for specific functional group
  
  # function that pulls the fillvalue
  get_nc_fillval <- function(varid) {
    fillval <- ncdf4::ncatt_get(nc = init_nc, varid = varid, attname = "fill.value")
    fillval_df <- data.frame(varid, fillval[2])
    fillval_df$species <- gsub("[0-9]", "",gsub("_StructN","",gsub("_ResN","",fillval_df$varid)))
    fillval_df$age <- as.numeric(gsub("[^0-9]", "",fillval_df$varid))
    #fillval_df$var <- sub(".*[0-9]_(.*)", "\\1", fillval_df$varid)
    fillval_df <- fillval_df %>% dplyr::select(species, age, value)
    return(fillval_df)
  }
  
  # Use pmap with only the variable names
  resN <- purrr::map(resN_vars$name, get_nc_fillval) %>% bind_rows()
  strucN <- purrr::map(strucN_vars$name, get_nc_fillval) %>% bind_rows()
  
  # combine, get total mgN, and convert to g 
  weight <- resN %>% 
    left_join(strucN, by = c("species","age")) %>% 
    mutate(totN = value.x+value.y,
           weight_g = totN * 20 * 5.7 / 1000) %>%
    dplyr::select(species,age,weight_g)
  
  # add upper and lower limit from gape size
  weight <- weight %>%
    mutate(low = weight_g * (g_pred %>% filter(limit == "low") %>% pull(g)),
           high = weight_g * (g_pred %>% filter(limit == "high") %>% pull(g)))
  
  # add stage and species code
  pred_weight <- weight %>%
    left_join(grps %>% dplyr::select(Code, Name), by = c("species"="Name")) %>%
    left_join(agemat, by = c("Code"="species")) %>%
    mutate(age = age -1,
           stage = ifelse(age < agemat, "juvenile", "adult")) %>%
    dplyr::select(species, stage, age, low, high)
  
  # now for the prey species, but by predator stage?
  eaten <- data.frame()
  for(i in 1:length(pred_stages)){
    
    fav_prey_df <- this_pprey %>%
      filter(PredatorAgeClass == pred_stages[i]) %>%
      group_by(Predator, PredatorAgeClass, PreyAgeClass) %>%
      arrange(desc(pprey)) %>%
      filter(pprey > 0) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    fav_prey <- fav_prey_df %>% dplyr::select(Prey) %>% distinct()
    
    if(length(fav_prey)==0) {next}
    
    prey_names <- fav_prey %>% 
      left_join(grps %>% dplyr::select(Code, Name), by = c("Prey"="Code")) %>% 
      pull(Name)
    
    pattern <- paste(prey_names, collapse = "|") # Create a single pattern string
    
    resN_vars <- hyper_vars(init) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(pattern,name)) # filter for specific functional group
    
    strucN_vars <- hyper_vars(init) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(pattern,name)) # filter for specific functional group
    
    # Use pmap with only the variable names
    resN <- purrr::map(resN_vars$name, get_nc_fillval) %>% bind_rows()
    strucN <- purrr::map(strucN_vars$name, get_nc_fillval) %>% bind_rows()
    
    # combine, get total mgN, and convert to g 
    weight <- resN %>% 
      left_join(strucN, by = c("species","age")) %>% 
      mutate(totN = value.x+value.y,
             weight_g = totN * 20 * 5.7 / 1000) %>%
      dplyr::select(species,age,weight_g)
    
    # add stage and species code
    weight <- weight %>%
      left_join(grps %>% dplyr::select(Code, Name), by = c("species"="Name")) %>%
      left_join(agemat, by = c("Code"="species")) %>%
      mutate(age = age -1,
             stage = ifelse(age < agemat, "juvenile", "adult")) %>%
      dplyr::select(species, stage, age, weight_g)
    
    # get age classes for this stage of the predator
    pred_weight_tmp <- pred_weight %>%
      filter(stage == pred_stages[i])
    
    pred_ages <- pred_weight_tmp %>% pull(age)
    
    # we need to do this by prey life stage or else they'll show in the table
    # do it by prey
    # and by life stage of the prey as informed in fav_prey_df
    for(j in 1:nrow(fav_prey)){
      
      this_prey <- fav_prey[j,1] %>% pull()
      
      # identify stages of the prey for this predator stage
      prey_stages <- fav_prey_df %>% filter(Prey == this_prey) %>% pull(PreyAgeClass)
      
      for(k in 1:length(prey_stages)){
        
        prey_weight_tmp <- weight %>%
          filter(species == prey_names[j], stage == prey_stages[k])
        
        # loop over predator's age classes
        
        for(l in 1:length(pred_ages)){
          
          lo <- pred_weight_tmp[l,]$low
          hi <- pred_weight_tmp[l,]$high
          
          ages_eaten <- prey_weight_tmp %>%
            mutate(eaten = ifelse(between(weight_g, lo, hi), 1, 0)) %>%
            filter(eaten > 0) %>%
            pull(age)
          

          eaten_tmp <- tibble(predator,
                              "stage" = pred_stages[i],
                              "age" = pred_ages[l],
                              "prey" = fav_prey[j,1] %>% pull(Prey),
                              "eaten_prey_ages" = list(ages_eaten))
          
          eaten <- rbind(eaten, eaten_tmp)
          
        }
        
        
      }
      
      
    }
    
  }
  
  eaten <- eaten %>% group_by(predator,stage,age,prey) %>%
    summarize(eaten_prey_ages = list(sort(unlist(eaten_prey_ages))))
  
  # rearrange and write out a table
  eaten_wide <- eaten %>%
    pivot_wider(id_cols = c(predator,stage,age), names_from = prey, values_from = eaten_prey_ages)
  
  # Function to replace numeric(0) with NA
  replace_numeric0_with_NA <- function(x) {
    if(is.list(x)) {
      lapply(x, function(y) if(length(y) == 0) "" else y)
    } else {
      x
    }
  }
  
  # Apply this function to each column of the data frame
  eaten_wide[] <- lapply(eaten_wide, replace_numeric0_with_NA)
  
  # sort
  eaten_wide <- eaten_wide %>% arrange(age)
  
  # Create a kable table
  gape_kable <- eaten_wide %>%
    knitr::kable(format = "markdown", booktabs = TRUE) 
  
  # Save the table as PDF
  # Using R Markdown
  # For this, you'll need to put your table creation code in an R Markdown document
  # and then knit the document to PDF.
  
  yaml_header <- paste0(
    "---\n",
    "title: '", predator, "s favorite prey items (with age class)'\n",
    "output: \n",
    "  pdf_document:\n",
    "    latex_engine: xelatex\n",
    "header-includes:\n",
    "  - '\\usepackage{graphicx}'\n",
    "  - '\\usepackage{grffile}'\n",
    "  - '\\tiny' # Set the font size\n",
    "geometry: landscape\n",
    "---\n"
  )
  
  g_file <- paste("diets_from_init", predator, "g.Rmd", sep="/")
  
  ##########################################################################################
  # write everything out
  # make a folder
  dir.create(paste("diets_from_init", predator, sep="/"), recursive = T)
  
  # diet preferences
  ggsave(d_file, d_plot, width = 8, height = 7)
  
  # horizontal distributions
  ggsave(s_file, s_plot, width = 25, height = s_length)
  
  # vertical distributions
  ggsave(v_file, v_plot, width = 10, height = 12)
  
  # gape-based interactions
  cat(yaml_header, gape_kable, sep="\n", file=g_file)
  rmarkdown::render(g_file, output_format = "pdf_document")
  
}

# apply function
purrr::map(verts, possibly(diet_from_init,NA))
