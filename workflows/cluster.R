###############################################################################
###                              kmeans                                     ###
###############################################################################
# A workflow to perform k-means clustering on sites


###############################################################################
### Generic Workflow Setup ###

# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")

### End Setup ###
###############################################################################

# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# SET USER INPUTS FOR WORKFLOW

sitecsv = "longsites.csv"
seq = c(729)
cluster_include_met = T
cluster_include_lags = F

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load libraries
library(ClusterR)

# Load sites
sites = read_csv(paste0("inputs/",
                        sitecsv),
                 col_names=FALSE,
                 show_col_types = F) %>% 
        unlist()

# Create inputs for parallel clustering
g <- expand.grid(sites=sites,
                 k=seq)

# Define the parallel worker function
worker = function(Site,k){
  
  # Load the input file and extract required data
  load(paste0("inputs/FLUXNET_processed/",
              Site,
              ".Rdata"))
  
  # Define climate variables
  metvars = c("LAI_MOD",
              "SWdown",
              "Tair",
              "Wind",
              "Precip",
              "VPD")
  
  lagvars = c("P1to30",
              "P31to90",
              "P91to180",
              "P181to365",
              "P366to730",
              "P731to1095",
              "P1096to1460",
              "T1",
              "T2to7",
              "T8to14",
              "T15to30")
  # Combine into clustering inputs based on requirements
  climvars = c(metvars[cluster_include_met],
               lagvars[cluster_include_lags])
  
  # Select and scale the requested data
  cluster_data = Data %>% 
    na.omit() %>% # Remove na rows from start
    select(all_of(climvars)) %>%
    scale() # Scale the data
  
  # Perform k-means clustering using the kmeans++ initialisation
  km = KMeans_rcpp(data = cluster_data,
                  clusters = k,
                  max_iters = 1000,
                  num_init = 50,
                  tol = 1e-05)
  
  # Attach the clustering inputs to the clustering output
  km$cluster_data = cluster_data
  
  # Save the file
  save(km,
       file=paste0("outputs/kmeans/clusterings/",
                    Site,
                   "_",
                   c("met")[cluster_include_met],
                   c("lag")[cluster_include_lags],
                   "_",
                   k,
                   "c.Rdata"))
  
  return(km)
}

# Perform k-means clustering in parallel
foo = mcmapply(FUN = function(Site,k) worker(Site,k),
               g$sites,
               g$k,
               mc.cores = 16)

###############################################################################
### Generic Workflow Cleanup ###

# Message end and elapsed time
End = Sys.time()
Duration = difftime(End,Start,unit = "mins")

message("Ending workflow at ",End)
message("Elapsed Time: ",round(Duration,2)," mins")

# Tidy up
rm(list=ls())

### End Cleanup ###
###############################################################################