
get_lmat <- function(this_genus,this_species){ # returns mean lmat for the species according to FishLife
  this_taxon <- fl_names[grepl(paste(this_genus, this_species,sep='.*'),fl_names,ignore.case = T)]
  lmat_this_species <- exp(fl_params[this_taxon,7]+0.5*fl_covmat[this_taxon,7,7]) # see code in tutorial.Rmd at https://github.com/James-Thorson-NOAA/FishLife/tree/main/vignettes
  return(lmat_this_species)
}

apply_Froese <- function(this_genus, this_species){
  
  Sys.setenv(JAGS_HOME='C:/Program Files/JAGS/JAGS-4.3.0')
  library(R2jags)  # Interface with JAGS
  #### R and JAGS code for estimating LWR-parameters from previous studies for a given species
  #### Developed by James Thorson in October 2012, modified by Rainer Froese in November/December 2012 
  #### Read in data
  #setwd('catch_to_CPUE_AKFIN/')  # Use current R working directory
  File = paste0(getwd(), '/')
  
  runif(1)         # sets random seed  
  
  DataFile = "BodyShape_3.csv" # The file BodyShape_3.csv contains over 5000 LWR studies
  Genus = this_genus  # generic name, e.g. "Engraulis"
  Species = this_species # specific name, e.g. "encrasicolus"
  Data = read.csv(DataFile, header=TRUE)
  DataG = Data[Data$Genus == Genus,] # Select data for Genus
  DataGS = DataG[DataG$Species == Species,] # select data for Species within Genus
  
  # Define data
  Keep = which(DataGS$Score>0) # exclude studies with Score = zero
  wts = DataGS$Score[Keep]  # Un-normalized weights (so that Cov is comparable among analyses)
  a = DataGS$a[Keep] # vector with estimates of parameter 'a' from selected studies
  b = DataGS$b[Keep] # vector with estimates of parameter 'b' from selected studies
  Bshape = DataGS$Bshape[Keep][1] # one of: "eel-like", "elongated", "fusiform", "short & deep"
  #Bshape = levels(Bshape)[as.numeric(Bshape)] # turns Bshape from type factor into text
  
  # if no body shape is indicated, use pre-estimated priors across all shapes 
  prior_mean_log10a = -2.0
  prior_sd_log10a = 0.313
  prior_tau_log10a = 1/prior_sd_log10a^2
  prior_mean_b = 3.04 
  prior_sd_b = 0.119
  prior_tau_b = 1/prior_sd_b^2
  if (Bshape == "eel-like") { # eel-like prior for log(a) and b
    prior_mean_log10a = -2.99
    prior_sd_log10a = 0.175
    prior_tau_log10a = 1/prior_sd_log10a^2
    prior_mean_b = 3.06 
    prior_sd_b = 0.0896
    prior_tau_b = 1/prior_sd_b^2}
  if (Bshape == "elongated") { # elongate prior for log(a) and b
    prior_mean_log10a = -2.41
    prior_sd_log10a = 0.171
    prior_tau_log10a = 1/prior_sd_log10a^2
    prior_mean_b = 3.12 
    prior_sd_b = 0.09
    prior_tau_b = 1/prior_sd_b^2}
  if (Bshape == "fusiform") { # fusiform prior for log(a) and b
    prior_mean_log10a = -1.95
    prior_sd_log10a = 0.173
    prior_tau_log10a = 1/prior_sd_log10a^2
    prior_mean_b = 3.04 
    prior_sd_b = 0.0857
    prior_tau_b = 1/prior_sd_b^2}
  if (Bshape == "short & deep") { # short & deep prior for log(a) and b
    prior_mean_log10a = -1.7
    prior_sd_log10a = 0.175
    prior_tau_log10a = 1/prior_sd_log10a^2
    prior_mean_b = 3.01 
    prior_sd_b = 0.0905
    prior_tau_b = 1/prior_sd_b^2}
  
  # Priors for measurement error (= sigma)
  # here given as mu and r for gamma distribution
  SD_rObs_log10a = 6520
  SD_muObs_log10a = 25076 
  SD_rObs_b = 6808
  SD_muObs_b = 37001
  
  # Define JAGS model 
  Model = "
model {               
  sigma1 ~ dgamma( SD_rObs_log10a, SD_muObs_log10a) # given the data and the prior, establish posterior distribution for measurement error in log10a  
  sigma2 ~ dgamma( SD_rObs_b, SD_muObs_b) # given the data and the prior, establish posterior distribution for measurement error in log10a   
  
  ro ~ dunif(-0.99,0)     # uniform prior for negative correlation between log10a and b
  abTrue[1] ~ dnorm(prior_mean_log10a,prior_tau_log10a) # given the data and the prior, get normal posterior distribution for log10a
  abTrue[2] ~ dnorm(prior_mean_b,prior_tau_b) # given the data and the prior, get normal posterior distribution for b
  CovObs[1,1] <- pow(sigma1,2)  
  CovObs[2,2] <- pow(sigma2,2) 
  CovObs[1,2] <- ro * sigma1 * sigma2 
  CovObs[2,1] <- CovObs[1,2]
  TauObs[1:2,1:2] <- inverse(CovObs[1:2,1:2]) # create inverse covariance matrix
  for(i in 1:N){
    TauObsI[i,1:2,1:2] <- TauObs[1:2,1:2] * pow(Weights[i],2)   # converts prior SD into prior weighted precision
    ab[i,1:2] ~ dmnorm(abTrue[1:2],TauObsI[i,1:2,1:2]) # given the data, the priors and the covariance, create multivariate normal posteriors for log(a) and b 
   }
}
"

# Write JAGS model 
cat(Model, file=paste(File,"dmnorm.bug",sep=""))
# JAGS settings
Nchains = 3	# number of MCMC chains to be used in JAGS
Nburnin = 1e4 # number of burn-in runs, to be discarded; 10000 iterations for burn-in
Niter = 3e4 # number of iterations after burn-in; 2e4 = 20000 iterations
Nthin = 1e1 # subset of iterations to be used for analysis; 1e1 = every 10th iteration 
# Run JAGS: define data to be passed on in DataJags; determine parameters to be returned in Param2Save; call JAGS with function Jags()
DataJags = list(ab=cbind(log10(a),b), N=length(a), Weights=wts, prior_mean_b=prior_mean_b, 
                prior_tau_b=prior_tau_b, prior_mean_log10a=prior_mean_log10a, prior_tau_log10a=prior_tau_log10a, 
                SD_rObs_log10a=SD_rObs_log10a, SD_muObs_log10a=SD_muObs_log10a,
                SD_rObs_b=SD_rObs_b, SD_muObs_b=SD_muObs_b)
Params2Save = c("abTrue","sigma1","sigma2","ro")
Jags <- jags(model.file=paste(File,"dmnorm.bug",sep=""), working.directory=NULL, data=DataJags, parameters.to.save=Params2Save, n.chains=Nchains, n.thin=Nthin, n.iter=Niter, n.burnin=Nburnin)
Jags$BUGSoutput # contains the results from the JAGS run
# Analyze output
abTrue <- Jags$BUGSoutput$sims.list$abTrue
mean_log10a  <- mean(abTrue[,1]) # true mean of log10(a)
sd_log10a    <- sd(abTrue[,1])   # true SE of log10(a)
mean_b       <- mean(abTrue[,2])         # true mean of b
sd_b         <- sd(abTrue[,2])           # true SE of b
mean_sigma_log10a <- mean(Jags$BUGSoutput$sims.list$sigma1) # measurement error of log10(a)
sd_sigma_log10a <- apply(Jags$BUGSoutput$sims.list$sigma1, 2, sd)
mean_sigma_b    <- mean(Jags$BUGSoutput$sims.list$sigma2) # measurement error of b
sd_sigma_b		<- apply(Jags$BUGSoutput$sims.list$sigma2, 2, sd)
ro_ab        <- mean(Jags$BUGSoutput$sims.list$ro) # measurement correlation of log10(a),b

return(c(mean_log10a, mean_b))
}

