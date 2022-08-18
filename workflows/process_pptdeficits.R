### Workflow to calculate long-term site metrics such as MAP and PET

### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")

### End Setup ###
###############################################################################

# Load sites
sites = read_csv(paste0("inputs/longsites.csv"),
                 col_names=FALSE,
                 show_col_types = F) %>% 
  unlist()

PPTdef = data.frame()

# Calculate the metrics for each site
for (Site in sites){

  # Load the FLUXNET data
  load(paste0("inputs/FLUXNET_processed/",Site,".Rdata"))
  
  # Calculate the mean rainfall for each month
  MMP =  Data %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Precip = sum(Precip,na.rm=T)) %>%
    ungroup() %>%
    mutate(Month = substr(YearMon,6,7)) %>%
    group_by(Month) %>%
    summarise(meanMP = mean(Precip),
              medianMP = median(Precip)) %>% 
    ungroup() %>%
    mutate(CumMeanP = cumsum(meanMP),
           CumMedianP = cumsum(medianMP))
  
  # Calculate actual monthly rainfall
  MP =  Data %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Precip = sum(Precip,na.rm=T)) %>%
    ungroup() %>%
    mutate(Year = substr(YearMon,1,4),
      Month = substr(YearMon,6,7)) %>%
    group_by(Year) %>%
    mutate(YearCum = cumsum(Precip))
  
  df = merge(MP,MMP) %>% select(-Month,-Year) %>%
    arrange(YearMon) %>%
    mutate(cumdeficitmean = YearCum-CumMeanP,
           cumdeficitmedian = YearCum-CumMedianP,
           deficitmean = Precip-meanMP,
           deficitmedian = Precip-medianMP)
  
  df$site = Site
  
  PPTdef = rbind(PPTdef,df)
}


# save(INFO,file = "inputs/siteinfo.Rdata")
save(PPTdef,file = "inputs/sitePPTdeficits.Rdata")
