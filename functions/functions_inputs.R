################################################################################

### Functions for stuff

################################################################################

calculateELI = function(Site, Data){
  # Function to calculate the ELI using various variables for energy, water and 
  # vegetation functioning
  # Inputs: 
  #       Site = site name
  #       Data = dataframe of inputs from FluxnetProcess
  len = 1:nrow(na.omit(Data))
  # Detrend data
  detrend = Data %>% na.omit() %>% 
    mutate(across(NEE:P1096to1460,
                  ~ .x - fitted(lm(.x ~ len, na.action = na.exclude))))
  # Calculate mean yearly cycle
  season_cycle = detrend %>% group_by(Month=month(Date),Day=day(Date)) %>%
    mutate(across(NEE:P1096to1460,mean)) %>%
    ungroup() %>%
    select(-Day,-Month)
  
  # Calculate anomalies
  anom = cbind("Date"=na.omit(Data)$Date,
               detrend[,2:ncol(Data)] - season_cycle[,2:ncol(Data)])
  
  # Split anomalies into Energy, Water and Veg variables
  Eanom = cbind("Tair" = anom$Tair,
                "VPD"=anom$VPD,
                "SWdown"=anom$SWdown)
  
  Wanom = cbind("PPT" = anom$Precip,
                "P1to30" = anom$P1to30)
  
  Vanom = cbind("NEE" = anom$NEE,
                "Qle" = anom$Qle,
                "GPP" = anom$GPP,
                "LAI" = anom$LAI_MOD)
  
  # Calculate the correlations
  EVcor = cor(Eanom,Vanom) 
  WVcor = cor(Wanom,Vanom)
  
  # Subtract and assign names
  ELI = lapply(1:3, function(x) EVcor[x,]-WVcor)
  names(ELI) = c("Tair","VPD","SWdown")
  
  ELI = ELI %>% as.data.frame() %>% 
    mutate("W" = rownames(.)) %>% 
    pivot_longer(!W, names_to = c("E","V"),names_sep = "\\.") %>% 
    mutate("Site" = Site)
  
  # Save
  save(ELI,file=paste0("inputs/ELI/",Site,"_ELI.Rdata"))
}