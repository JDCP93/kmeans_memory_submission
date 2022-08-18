FluxnetProcess = function(Site){

  # A function to extract the necessary data from the Fluxnet netcdf for a 
  # given site. This is designed to work with data that has been cleaned by Anna Ukkola. 
  # 
  # ############################################################################
  # Function inputs and outputs
  # ############################################################################
  #
  # INPUTS:
  #  - Site: A character vector with the FluxNet siteID
  #  
  #  OUTPUTS:
  #  - Site.Rdata: a .Rdata file containing "Data", a dataframe of daily met
  #                     and flux data, as well as antecedent rainfall.
  
  library(ncdf4)
  # ##################
  # Extract raw data from files
  # ##################
  # 
  # Load in the Flux data we need:
  Files = list.files("inputs/FLUXNET_raw/",pattern = Site)
  # Read the data into R 
  FluxNC = nc_open(paste0("inputs/FLUXNET_raw/",Files[1]))
  
  # Convert time into datetimes
  # Find the origin of the times
  time_from = substr(ncatt_get(FluxNC, "time")$units, 15, 33)
  # times are in seconds since origin
  Time = ncvar_get(FluxNC,"time") %>%
          as.POSIXct(origin=time_from, tz = "GMT")

  # List the variables we want to extract
  Variables = c("NEE",
                "GPP",
                "Qle",
                "Qh")
  # Extract the variables we require
  Data = data.frame(Time)
  for (Var in Variables){
    Data[Var] = ncvar_get(FluxNC,Var)
  }
  
  # Met data
  MetNC = nc_open(paste0("inputs/FLUXNET_raw/",Files[2]))
  # List the variables we want to extract 
  Variables = c("Tair",
                "SWdown",
                "VPD",
                "Precip",
                "LAI",
                "LAI_alternative",
                "Wind",
                "Psurf",
                "RH",
                "latitude",
                "elevation")
  # Extract the variables we require
  for (Var in Variables){
    Data[Var] = ncvar_get(MetNC,Var)
  }
  
  # We assume all data is good-quality since Anna Ukkola gapfilled this
  # However, still omit any NA rows from the beginning
  if(rle(complete.cases(Data))$values[1] == F){
    Data = Data[-(1:rle(complete.cases(Data))$lengths[1]),]
  }
  
  # Retime data to daily 
  # multiply rates by 1800 to go from 1/s to 1/30mins before summing to 1/day
  Data = Data %>%
          mutate(Date = as.Date(Time)) %>%
          group_by(Date) %>%
          summarise(NEE=sum(NEE)*1800, # umol/m2/s to umol/m2day 
                    GPP=sum(GPP)*1800, # umol/m2/s to umol/m2day 
                    Qle=sum(Qle)*1800, # J/m2s to J/m2day
                    Qh=sum(Qh)*1800,  # J/m2s to J/m2day
                    Tmax=max(Tair), # K
                    Tmin=min(Tair), # K
                    Tair=mean(Tair), # K
                    SWdown=sum(SWdown)*1800, # J/m2s to J/m2day
                    VPD=mean(VPD)/10, # hPa to kPa
                    Precip=sum(Precip)*1800, # mm/s to mm/day
                    LAI_MOD=mean(LAI), # m/m
                    LAI_CLS=mean(LAI_alternative), # m/m
                    Wind=mean(Wind),# m/s
                    Psurf = mean(Psurf), # Pa
                    RH = mean(RH), # %
                    Lat = mean(latitude),
                    Elevation = mean(elevation)) 
  
  # PPT lags
  # for each day in the main data series, calculate the antecedent rainfall
  for (i in 15:nrow(Data)){
    Data[i,"P1to14"] = sum(Data$Precip[(i-1):(i-14)])
    if(i>30){Data[i,"P15to30"] = sum(Data$Precip[(i-15):(i-30)])}
    if(i>60){Data[i,"P31to60"] = sum(Data$Precip[(i-31):(i-60)])}
    if(i>120){Data[i,"P61to120"] = sum(Data$Precip[(i-61):(i-120)])}
    if(i>180){Data[i,"P121to180"] = sum(Data$Precip[(i-121):(i-180)])}
    if(i>270){Data[i,"P181to270"] = sum(Data$Precip[(i-181):(i-270)])}
    if(i>365){Data[i,"P271to365"] = sum(Data$Precip[(i-271):(i-365)])}
    if(i>30){Data[i,"P1to30"] = sum(Data$Precip[(i-1):(i-30)])}
    if(i>90){Data[i,"P31to90"] = sum(Data$Precip[(i-31):(i-90)])}
    if(i>180){Data[i,"P91to180"] = sum(Data$Precip[(i-91):(i-180)])}
    if(i>365){Data[i,"P181to365"] = sum(Data$Precip[(i-181):(i-365)])}
    if(i>730){Data[i,"P366to730"] = sum(Data$Precip[(i-366):(i-730)])}
    if(i>1095){Data[i,"P731to1095"] = sum(Data$Precip[(i-731):(i-1095)])}
    if(i>1460){Data[i,"P1096to1460"] = sum(Data$Precip[(i-1096):(i-1460)])}
  }

  # Tair lags
  # for each day in the main data series, calculate the antecedent temperature
  for (i in 2:nrow(Data)){
    Data[i,"T1"] = mean(Data$Tair[(i-1)])
    if(i>7){Data[i,"T2to7"] = sum(Data$Tair[(i-2):(i-7)])}
    if(i>14){Data[i,"T8to14"] = sum(Data$Tair[(i-8):(i-14)])}
    if(i>21){Data[i,"T15to21"] = sum(Data$Tair[(i-15):(i-21)])}
    if(i>30){Data[i,"T22to30"] = sum(Data$Tair[(i-22):(i-30)])}
    if(i>7){Data[i,"T1to7"] = sum(Data$Tair[(i-1):(i-7)])}
    if(i>30){Data[i,"T15to30"] = sum(Data$Tair[(i-15):(i-30)])}
    if(i>60){Data[i,"T31to60"] = sum(Data$Tair[(i-31):(i-60)])}
    if(i>180){Data[i,"T61to180"] = sum(Data$Tair[(i-61):(i-180)])}
    if(i>365){Data[i,"T181to365"] = sum(Data$Tair[(i-181):(i-365)])}
    if(i>730){Data[i,"T366to730"] = sum(Data$Tair[(i-366):(i-730)])}
  }
  
  # Calculate Ecosystem Limitation Indices (ELIs)
  calculateELI(Site, Data)
  
  # Plot and save the Inputs
  InputPlot(Site,Data)
  
  # Create output
  save(Data,file=paste0("inputs/FLUXNET_processed/",Site,".Rdata"))
  
  message("Input processed for ",Site)
}
