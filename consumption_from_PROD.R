# Read PROD.nc

library(tidyverse)
library(rbgm)
library(ggsidekick)
library(maps)
library(mapdata)
library(tidync)
library(sf)
library(viridis)
library(ggh4x)
library(RColorBrewer)

select <- dplyr::select

# make new folder to save plots in
now <- gsub(' ', '_', gsub(':', '.', Sys.time()))

run_base <- 1328 # this is the control run we compare results to
run_new <- 1329 # just warm (temp) forcings

# set paths to directories
dir_base <- paste0('../../output_files/data/out_', run_base, '/')
dir_new <- paste0('../../output_files/data/out_', run_new, '/')

# Read files --------------------------------------------------------------
# Geography: read bgm
fl <- '../data/GOA_WGS84_V4_final.bgm'
bgm <- read_bgm(fl)
goa_sf <- box_sf(bgm)
st_crs(goa_sf) <- st_crs(attr(goa_sf$geometry, "crs")$proj)
boundary_boxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes

# Biology: read Groups.csv
grps <- read.csv('../data/GOA_Groups.csv', header = T)

# set up a functional group types table
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

# producers and detritus
prods_gtype <- c('PHYTOBEN', 'LG_PHY', 'SM_PHY', 'SED_BACT', 'PL_BACT', 'CARRION', 'LAB_DET', 'REF_DET')
# predators:
preds_gtype <- setdiff((grps %>% pull(GroupType) %>% unique()), prods_gtype)

# and corresponding codes
pred_codes <- grps %>% filter(GroupType %in% preds_gtype) %>% pull(Code)

# age class at maturity to split age classes between juveniles and adults
agemat <- read.csv('../data/agemat.csv', header = T)

# read in biomass CSV
biomage_file <- list.files(dir_base, pattern = "AgeBiom", full.names = T)[1]
biomage <- read.table(biomage_file, sep = " ", header = T)

# Model output: read NetCDF files
# We need both the PROD and the regular nc output files
# base run
# regular
out_fl_base <- paste0(dir_base, 'outputGOA0', run_base, '_test.nc')
out_base <- tidync(out_fl_base)
this_nc_base <- ncdf4::nc_open(out_fl_base)
# PROD.nc
out_fl_base_prod <- paste0(dir_base, 'outputGOA0', run_base, '_testPROD.nc')
out_base_prod <- tidync(out_fl_base_prod)
this_nc_base_prod <- ncdf4::nc_open(out_fl_base_prod)
# # new run
# out_fl_new <- paste0(dir_new, 'outputGOA0', run_new, '_testPROD.nc')
# out_new <- tidync(out_fl_new)
# this_nc_new <- ncdf4::nc_open(out_fl_new)

# derived values for output
depths <- out_base %>% hyper_filter(t=t==0) %>% hyper_tibble(select_var="dz") %>% dplyr::select(-t)
# glimpse(depths)

# # volumes of each layer
# # extract from base run
volumes <- out_base %>% hyper_filter(t=t==0) %>% hyper_tibble(select_var="volume") %>% dplyr::select(-t)
# careful with indexing of volumes...
box_volumes <- volumes %>% mutate(b = b-1) %>% group_by(b) %>% summarize(volume = sum(volume))

#
# # time dimension
# # extract from base run
ts <- ncdf4::ncvar_get(this_nc_base,varid = "t") %>% as.numeric
tyrs <- ts/(60*60*24*365)

areas <- volumes %>% filter(z==max(z)) %>% dplyr::select(b,volume) %>% rename(area=volume)

# read consumption time series by age class
# this also needs to pull in numbers, RN, and SN from the out.nc file, to apportion consumption to box layers
plot_eat_timeseries <- function(fg, out, outprod, this.nc, this.nc.prod, run){
  
  # These are used for testing / debugging
  # fg = "Pollock"
  # out = out_base
  # outprod = out_base_prod 
  # this.nc = this_nc_base
  # this.nc.prod = this_nc_base_prod 

  print(paste("Doing", fg))
  
  # get the attributes associated with each functional group
  fg_atts <- grps %>% filter(Name==fg)
  
  if(fg_atts$BiomassType!="vertebrate") stop("weight at age only for vertebrates.")
  
  #Extract from the output PROD.nc file the Eat variables per age class
  eat_vars <- outprod %>%
    activate("D1,D0") %>%
    hyper_vars() %>% # all variables in the .nc file active grid
    filter(grepl("_Eat",name)) %>% # filter for Eat
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate reserve N time series variables
  resN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_ResN",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate structural N time series variables
  strucN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_StructN",name)) %>% # filter for structural N
    filter(grepl(fg,name)) # filter for specific functional group
  
  # Get numbers by box
  abun_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
    filter(grepl("_Nums",name)) %>% # filter for abundance variables
    filter(grepl(fg,name)) # filter for specific functional group
  
  if(nrow(eat_vars)==0) {return("no data.")}
  else {
    # # Actually pull the data from the .nc files
    eat <- purrr::map(eat_vars$name,ncdf4::ncvar_get,nc=this.nc.prod) 
    resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this.nc) 
    strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this.nc)
    nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this.nc) #numbers by age group,box,layer,time
    
    # # add the two weight matrices to get total nitrogen weight
    rnsn <- purrr::map2(resN,strucN,`+`)
    
    # check the sanity of this step
    # check1 <- resN[[1]][4,,] # non-0 values in boundary boxes because of fillvals in init.nc, but accounted for in nums which are 0 in BB's
    # check2 <- strucN[[1]][4,,]
    # check3 <- rnsn[[1]][4,,] # seems all good

    # # multiply by numbers to get total biomass per cell (unit is irrelevant here because we will use this to calculate proportions)
    biomass <- purrr::map2(rnsn,nums,`*`)
    
    # check this step
    # throughout, remember that layer 7 is sediment, layer 6 is the surface, and so on...
    # check4 <- nums[[1]][4,,]
    # check5 <- biomass[[1]][4,,] # this is sane
    
    # now we have a list for eat and a list for biomass
    # could collapse the dimensions here and work with data frames instead, but it will be slow.
    # GOAL: take Q by ts by box, and apportion it to layers
    # Initialize a list or array to store the results
    distributed_eat <- list()
    
    # Loop through each age class
    for (i in 1:length(eat)) {

      # Calculate sum of biomass across depth layers (7, 109, 251 -> 109, 251)
      sum_biomass <- apply(biomass[[i]], c(2, 3), sum)
      
      # # check
      # check6 <- biomass[[1]][1,,] +
      #   biomass[[1]][2,,] +
      #   biomass[[1]][3,,] +
      #   biomass[[1]][4,,] +
      #   biomass[[1]][5,,] +
      #   biomass[[1]][6,,] +
      #   biomass[[1]][7,,] -
      #   sum_biomass 
      # max(check6) # this works
        
      # Expand sum_biomass to a 3D array to match dimensions of biomass[[i]]
      expanded_sum_biomass <- array(rep(sum_biomass, each = 7), dim = dim(biomass[[i]]))
      
      # check7 <- expanded_sum_biomass[1,,]
      # check8 <- expanded_sum_biomass[4,,] # this works
      
      # Calculate proportions
      proportions <- biomass[[i]] / expanded_sum_biomass
      
      # check, these need to add up to 1 over the water column
      # check9 <- proportions[1,,] + 
      #   proportions[2,,] + 
      #   proportions[3,,] + 
      #   proportions[4,,] +
      #   proportions[5,,] +
      #   proportions[6,,] +
      #   proportions[7,,]
      # min(check9, na.rm = T) # this works, all proportions add up to 1 across the water column
      
      
      # Initialize an array to store the distributed values for this age class
      distributed_values <- array(dim = dim(biomass[[i]]))
      
      # Distribute `eat` values across depth layers
      for (j in 1:7) {
        distributed_values[j,,] <- eat[[i]] * proportions[j,,]
      }
      
      # Store the distributed values in the result structure
      distributed_eat[[i]] <- distributed_values
    }
    
    # check
    # check10 <- distributed_eat[[5]][1,,] +
    #   distributed_eat[[5]][2,,] +
    #   distributed_eat[[5]][3,,] +
    #   distributed_eat[[5]][4,,] +
    #   distributed_eat[[5]][5,,] +
    #   distributed_eat[[5]][6,,] +
    #   distributed_eat[[5]][7,,] -
    #   eat[[5]]
    # min(check10, na.rm = T) # this works, they all sum up to the original eat values
    
    # this contains values in mg N m-3 d-1, but disaggregated by depth layer based on biomass per depth layer
    
    # mind that layers are indexed like 7 = sediment, 6 = surface, etc, so make sure that volumes match this
    # Convert the volume data frame from long to wide format
    # remember, layer 7 is the sediment
    volume_matrix <- matrix(nrow = 7, ncol = 109)
    for (i in 1:nrow(volumes)) {
      volume_matrix[volumes$z[i], volumes$b[i]] <- volumes$volume[i]
    }
    
    # check
    # check11 <- distributed_eat[[1]][,,4] # volumes matrix and distributed_eat slices have the correct shape
    
    # Iterate through each element in the distributed_eat list
    for (i in 1:length(distributed_eat)) {
      # Multiply each value by the corresponding volume
      for (j in 1:dim(distributed_eat[[i]])[3]) {  # iterating over time steps
        distributed_eat[[i]][,,j] <- distributed_eat[[i]][,,j] * volume_matrix
      }
    }
    
    # Now distributed_eat contains values in mg N d-1
    # let's go to mt d-1
    mgn_to_mt <- 20*5.7/1000000000
    
    # Iterate through each element in the distributed_eat list
    for (i in 1:length(distributed_eat)) {
      # Multiply the entire array by the scalar
      distributed_eat[[i]] <- distributed_eat[[i]] * mgn_to_mt
    }
    
    # now we have values in mt d-1, by cell and by age class.
    # let's turn this into a data frame to better work with it
    eat_df <- data.frame()
    
    # Iterate through each age class in distributed_eat
    for (age in 1:length(distributed_eat)) {
      # Get the array for this age class
      eat_array <- distributed_eat[[age]]
      
      # Flatten the array into a matrix
      # stack time steps on top of each other
      flattened_matrix <- apply(eat_array,2,c) # rows are rep(1:7, 251)
      
      # check
      # check12 <- eat_array[,,1]
      # check13 <- eat_array[,,2]
      # check14 <- eat_array[,,3]
      # check15 <- head(flattened_matrix, 3*7) # concatenation works ok
      
      # Convert the matrix to a data frame
      temp_df <- as.data.frame(flattened_matrix)
      
      # Add columns for age, depth, location, and time step
      temp_df <- temp_df %>%
        `colnames<-`(0:(ncol(.)-1)) %>% # name columns after box id (numbered from 0)
        mutate(age = age,
               depth = rep(1:dim(eat_array)[1], dim(eat_array)[3]),
               ts = rep(1:dim(eat_array)[3], each = dim(eat_array)[1])) %>%
        pivot_longer(-c(age,depth,ts), names_to = "b", values_to = "daily_eat_mt") %>%
        mutate(b = as.numeric(b),
               daily_eat_mt = replace(daily_eat_mt, is.nan(daily_eat_mt), 0)) # turn NaN to 0
      
      # now sum across all cells (depth layers and boxes) because we are dealing with tons and not values per m3
      temp_df_box <- temp_df %>%
        group_by(age,ts) %>%
        summarise(daily_eat_mt = sum(daily_eat_mt))
      
      # now bind the age dataframes
      eat_df <- rbind(eat_df, temp_df_box)
    }
    
    # now interpolate to get daily values from 73-day output. We need daily values if we want to estimate annual consumption
    # Eat may not change linearly between time steps so this has error aound it
    # make complete df
    eat_df_complete <- eat_df %>%
      mutate(ts = ts-1) %>% # number from 0
      filter(ts > 0) %>% # drop first time step (values are 0)
      mutate(ts = ts * 73) %>% # effectively start from day 73
      group_by(age) %>%
      complete(ts = full_seq(ts, 1), fill = list(daily_eat_mt = NA)) %>%
      ungroup()
    
    # interpolate
    eat_df_daily <- eat_df_complete %>%
      group_by(age) %>%
      arrange(age) %>%
      mutate(daily_eat_mt = approx(ts, daily_eat_mt, xout = ts)$y) %>%
      ungroup()
    
    # add species name and turn back to years
    eat_df_daily <- eat_df_daily %>%
      mutate(Name = fg_atts$Name,
             LongName = fg_atts$LongName) %>%
      mutate(ts = ts / 365)
    
    # drop burn-in period
    eat_df_daily <- eat_df_daily %>% filter(ts >= 30)
    
    # final_df now contains all the data from distributed_eat in a long format
    
    return(eat_df_daily)
  }
}

# make plots with time series of lines by age class, and one with the sum of all age classes
# plot eat time series
fg_to_plot <- grps %>% filter(GroupType %in% c("MAMMAL","BIRD","SHARK","FISH")) %>% pull(Name)
codes_to_plot <- grps %>% filter(GroupType %in% c("MAMMAL","BIRD","SHARK","FISH")) %>% pull(Code)
eat_base <- bind_rows(purrr::map(fg_to_plot, plot_eat_timeseries, 
                                 out = out_base,
                                 outprod = out_base_prod, 
                                 this.nc = this_nc_base,
                                 this.nc.prod = this_nc_base_prod, 
                                 run = 'base'))

# plot
# this out the series to plot it though
p <- eat_base %>%
  filter(ts %in% 30:50) %>%
  ggplot(aes(x = ts, y = daily_eat_mt, color = factor(age)))+
  geom_line()+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = "Year", y = "Daily consumption (mt)", color = "Age")+
  facet_wrap(~LongName, scales = "free", nrow = 10)
p

ggsave("consumption_from_PROD.png", p, width = 8, height = 16)

# now sum all age classes and get annual consumption
annual_consumption <- eat_base %>%
  mutate(year = ceiling(ts)) %>%
  group_by(LongName, year) %>%
  summarise(Q = sum(daily_eat_mt))

# get biomass
biomage1 <- biomage %>% 
  pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
  separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
  filter(Code %in% codes_to_plot) %>%
  group_by(Time, Code) %>%
  summarise(biomass_mt = sum(biomass_mt)) %>%
  ungroup() %>%
  left_join(grps %>% dplyr::select(Code, LongName), by = "Code") %>%
  mutate(year = Time / 365)

# join and get QB
# look at Tier 3 species
qb_t3 <- annual_consumption %>%
  left_join(biomage1, by = c("LongName", "year")) %>%
  filter(Code %in% c("POL","COD","ATF","POP","HAL","FFS","FFD","FHS","REX","SBF","RFS","RFP")) %>%
  mutate(QB = Q/biomass_mt) 

qb_t3 %>%
  ggplot(aes(x = year, y = QB, color = LongName))+
  geom_line()+
  theme_bw()

# The method that apportions the "Eat" tracer to depth layers by fish biomass returns low QB values (0.2-0.7 for Tier 3 stocks)
