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
  # Calculate P-T PET
  source("functions/functions_PET.R")
  #Convert SWdown from W m-2 to MJ/m2/day
  conv_mj_m2_d <- 86400 * 0.001 * 0.001
  swrad = Data$SWdown/1800/24
  lat = Data$Lat[1]
  temp = Data$Tair
  #Calculate fractional sunshine hours from SWdown
  sun_frac <- sw2sunshine(day = rep(1:365, length.out=length(swrad)), lat=lat, swrad = swrad * conv_mj_m2_d)
  
  #Calculate daily PET 
  PET_pt <- daily_pet(T=temp-273.15, clt=1-sun_frac, day=rep(1:365, length.out=length(swrad)), lat=lat)
  
  Data$PET_pt = PET_pt
  
  # Turn the daily data into monthly data
  foo =  Data %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Tmax = max(Tmax)-273.15, # Need in Celsius
              Tave = mean(Tair,na.rm=T)-273.15,
              Tmin = min(Tmin)-273.15,
              U2 = mean(Wind,na.rm=T),
              Rs = mean(SWdown,na.rm=T)/100000, # Need in MJ/m2
              P = mean(Psurf,na.rm=T),
              RH = mean(RH,na.rm=T),
              Precip = sum(Precip,na.rm=T),
              ET = (sum(Qle,na.rm=T)/2450)/1000, # Turn from J/m2 to g/m2 to kg/m2=mm
              PET_pt = sum(PET_pt)) 
  
  # Calculate the potential evapotranspiration                  
  foo$PET_th = as.numeric(thornthwaite(foo$Tave,lat=Data$Lat[1]))
  foo$PET_pm = as.numeric(penman(foo$Tmin,foo$Tmax,foo$U2,lat = Data$Lat[1],Rs = foo$Rs, z = Data$Elevation[1], P = foo$P, RH = foo$RH))
  
  # Calculate the yearly data
  bar = foo %>%  mutate(Year = substr(YearMon,1,4)) %>%
    group_by(Year) %>%
    summarise(AP = sum(Precip,na.rm=T),
              AT = mean(Tave,na.rm=T),
              APET_th = sum(PET_th,na.rm=T),
              APET_pm = sum(PET_pm,na.rm=T),
              APET_pt = sum(PET_pt,na.rm=T),
              CVP = sd(Precip,na.rm=T)/mean(Precip,na.rm=T),
              AET = sum(ET,na.rm=T))
  
  foo$site = Site
  
  # Calculate the yearly statistics
  INFO$MAP[INFO$Site == Site] = colMeans(bar[,2],na.rm=T)
  INFO$MAT[INFO$Site == Site] = colMeans(bar[,3],na.rm=T)
  INFO$MAPET_th[INFO$Site == Site] = colMeans(bar[,4],na.rm=T)
  INFO$MAPET_pm[INFO$Site ==Site] = colMeans(bar[,5],na.rm=T)
  INFO$MAPET_pt[INFO$Site ==Site] = colMeans(bar[,6],na.rm=T)
  INFO$CVP[INFO$Site == Site] = colMeans(bar[,7],na.rm=T)
  INFO$MAET[INFO$Site == Site] = colMeans(bar[,8],na.rm=T)
  INFO$Mths[INFO$Site == Site] = nrow(foo)
  INFO$Yrs[INFO$Site == Site] = nrow(bar)
  
  PPTSeries = rbind(PPTSeries,foo)
}

INFO$AI = INFO$MAPET_th/INFO$MAP
INFO$AI_pm = INFO$MAPET_pm/INFO$MAP
INFO$AI_pt = INFO$MAPET_pt/INFO$MAP
INFO$EI = INFO$MAET/INFO$MAP

save(INFO,file = "inputs/siteinfo.Rdata")
save(PPTSeries,file = "inputs/sitePPTseries.Rdata")
