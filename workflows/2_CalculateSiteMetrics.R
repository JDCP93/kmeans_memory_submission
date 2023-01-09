### Workflow to calculate long-term site metrics such as MAP and PET

### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Source library
library(SPEI)

# Make sure directories and libraries are set up
source("workflows/setup.R")

### End Setup ###
###############################################################################

# Load sites
sites = read_csv(paste0("inputs/longsites.csv"),
                 col_names=FALSE,
                 show_col_types = F) %>% 
  unlist()

# Define dataframe
INFO = data.frame("Site" = sites,
                  MAP=0,
                  MAT=0,
                  MAPET=0,
                  Mths=0,
                  Yrs=0)

PPTSeries = data.frame()

# Calculate the metrics for each site
for (Site in sites){
  # Load the FLUXNET data
  load(paste0("inputs/FLUXNET_processed/",Site,".Rdata"))
  
  print(Site)
  # Turn the daily data into monthly data
  MonthlyData =  Data %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Tmax = max(Tmax)-273.15, # Need in Celsius
              Tave = mean(Tair,na.rm=T)-273.15,
              Tmin = min(Tmin)-273.15,
              U2 = mean(Wind,na.rm=T),
              Rs = mean(SWdown,na.rm=T)/100000, # Need in MJ/m2
              P = mean(Psurf,na.rm=T),
              RH = mean(RH,na.rm=T),
              Precip = sum(Precip,na.rm=T),
              ET = (sum(Qle,na.rm=T)/2450)/1000) # Turn from J/m2 to g/m2 to kg/m2=mm
  
  # Calculate the potential evapotranspiration                  
  MonthlyData$PET = as.numeric(thornthwaite(MonthlyData$Tave,lat=Data$Lat[1]))
  
  # Calculate the yearly data
  YearlyData = MonthlyData %>%  mutate(Year = substr(YearMon,1,4)) %>%
    group_by(Year) %>%
    summarise(AP = sum(Precip,na.rm=T),
              AT = mean(Tave,na.rm=T),
              APET = sum(PET,na.rm=T),
              CVP = sd(Precip,na.rm=T)/mean(Precip,na.rm=T),
              AET = sum(ET,na.rm=T))
  
  MonthlyData$site = Site
  
  # Calculate the yearly statistics
  INFO$MAP[INFO$Site == Site] = colMeans(YearlyData[,2],na.rm=T)
  INFO$MAT[INFO$Site == Site] = colMeans(YearlyData[,3],na.rm=T)
  INFO$MAPET[INFO$Site == Site] = colMeans(YearlyData[,4],na.rm=T)
  INFO$CVP[INFO$Site == Site] = colMeans(YearlyData[,5],na.rm=T)
  INFO$MAET[INFO$Site == Site] = colMeans(YearlyData[,6],na.rm=T)
  INFO$Mths[INFO$Site == Site] = nrow(MonthlyData)
  INFO$Yrs[INFO$Site == Site] = nrow(YearlyData)
  
  PPTSeries = rbind(PPTSeries,MonthlyData)
}

# Calculate the metrics
INFO$AI = INFO$MAPET/INFO$MAP
INFO$EI = INFO$MAET/INFO$MAP

# Save the outputs
save(INFO,file = "inputs/siteinfo.Rdata")
save(PPTSeries,file = "inputs/sitePPTseries.Rdata")
