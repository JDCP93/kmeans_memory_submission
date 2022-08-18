###############################################################################
###                         FluxnetProcess                                  ###
###############################################################################
# A workflow to extract the data from PLUMBER2 netcdfs and turn it into an Rdata
# file of a dataframe


###############################################################################
### Generic Workflow Setup ###

# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")

# Source functions
source("functions/function_fluxnetprocess.R")

### End Setup ###
###############################################################################

# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# SET USER INPUTS FOR WORKFLOW

sitecsv = "longsites.csv"
check_if_file_exists = F

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

### WORKFLOW 

# Load sites
sites = read_csv(paste0("inputs/",sitecsv),col_names=FALSE) %>%
        as.matrix()

# Initialise input dataframe
df = data.frame()

# Process the raw data
for (site in sites){
 if (!file.exists(paste0("inputs/FLUXNET_processed/",site,".Rdata")) &
     check_if_file_exists == T){
    FluxnetProcess(site)
  } else if (file.exists(paste0("inputs/FLUXNET_processed/",site,".Rdata")) &
                check_if_file_exists == T){
    message(paste0("File exists for ",site))
  } else {
    FluxnetProcess(site)
  }
  # Reload the data
  load(paste0("inputs/FLUXNET_processed/",site,".Rdata"))
  # Place into combined dataframe
  Data$Site = site
  df = rbind(df,Data)
}

# Save combined dataframe
save(df,file = "inputs/FLUXNET_processed/allsites.Rdata")

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