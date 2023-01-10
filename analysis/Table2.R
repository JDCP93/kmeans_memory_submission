### Workflow to plot timeseries and obs vs pred for sites

### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")
library(zoo)


### End Setup ###
###############################################################################
###############################################################################

k = 729

###############################################################################
###############################################################################


# Load sites
sites = read_csv(paste0("inputs/longsites.csv"),
                 col_names=FALSE,
                 show_col_types = F) %>% 
  unlist()

All.All = data.frame()
All.Monthly = data.frame()
All.Daily = data.frame()
Metrics = data.frame()

for (Site in sites){
  
  # Load the instantaneous model output
  load(paste0("outputs/kmeans/NEE/",
              Site,
              "_metcluster_metregress_",
              k,
              "c.Rdata"))
  
  instantaneous = output$Flux_df %>% select(Date,Flux_pred)
  names(instantaneous) = c("Date","Instantaneous")
  
  # Load the historical model output
  load(paste0("outputs/kmeans/NEE/",          
              Site,
              "_metlagcluster_metlagregress_",
              k,
              "c.Rdata"))
  
  historical = output$Flux_df %>% select(Date,Flux,Flux_pred)
  names(historical) = c("Date","Observed","Historical")
  
  # Extract the scaling data from before modelling occurred
  NEE.df.center = attr(historical$Observed,c("scaled:center"))
  NEE.df.scale = attr(historical$Observed,c("scaled:scale"))
  
  NEE.df = merge(instantaneous,historical,by.x = "Date", by.y = "Date")
  NEE.df$Flux = "NEE"
  # Take negative to turn NEE into NEP
  NEE.df[,c("Observed","Instantaneous","Historical")] = -NEE.df[,c("Observed","Instantaneous","Historical")]
  
  # Load the instantaneous model output
  load(paste0("outputs/kmeans/Qle/",
              Site,
              "_metcluster_metregress_",
              k,
              "c.Rdata"))
  
  instantaneous = output$Flux_df %>% select(Date,Flux_pred)
  names(instantaneous) = c("Date","Instantaneous")
  
  # Load the historical model output
  load(paste0("outputs/kmeans/Qle/",          
              Site,
              "_metlagcluster_metlagregress_",
              k,
              "c.Rdata"))
  
  historical = output$Flux_df %>% select(Date,Flux,Flux_pred)
  names(historical) = c("Date","Observed","Historical")
  
  Qle.df.center = attr(historical$Observed,c("scaled:center"))
  Qle.df.scale = attr(historical$Observed,c("scaled:scale"))
  
  Qle.df = merge(instantaneous,historical,by.x = "Date", by.y = "Date")
  Qle.df$Flux = "Qle"
  
  df = rbind(NEE.df,Qle.df)
  
  # Calculate the monthly statistics for NEE
  MonthAvg.NEE = NEE.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Observed_month_mean = mean(Observed*Qle.df.scale+Qle.df.center),
              Instantaneous_month_mean = mean(Instantaneous*Qle.df.scale+Qle.df.center),
              Historical_month_mean = mean(Historical*Qle.df.scale+Qle.df.center)) %>%
    mutate(Mon = substr(YearMon,6,7)) %>%
    group_by(Mon) %>%
    summarise(Observed_mean = mean(Observed_month_mean),
              Observed_sd = sd(Observed_month_mean),
              Instantaneous_mean = mean(Instantaneous_month_mean),
              Instantaneous_sd = sd(Instantaneous_month_mean),
              Historical_mean = mean(Historical_month_mean),
              Historical_sd = sd(Historical_month_mean)) %>%
    mutate(Date = as.Date(paste0("2000-",Mon,"-01",format = "%Y-%m-%d")),
           "Flux" = "NEE",
           Observed_sd = ifelse(is.na(Observed_sd), 0.1*Observed_mean, Observed_sd),
           Instantaneous_sd = ifelse(is.na(Instantaneous_sd), 0.1*Instantaneous_mean, Instantaneous_sd),
           Historical_sd = ifelse(is.na(Historical_sd), 0.1*Historical_mean, Historical_sd),
           InstantPoorMean = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean-Observed_sd 
                                      | Instantaneous_mean-Instantaneous_sd > Observed_mean+Observed_sd),T,F),
           HistoricalPoorMean = if_else((Historical_mean+Historical_sd < Observed_mean-Observed_sd 
                                         | Historical_mean-Historical_sd > Observed_mean+Observed_sd),T,F),
           InstantPoorMean_percent = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean*0.9 
                                              | Instantaneous_mean-Instantaneous_sd > Observed_mean*1.1),T,F),
           HistoricalPoorMean_percent = if_else((Historical_mean < Observed_mean*0.9
                                                 | Historical_mean-Historical_sd > Observed_mean*1.1),T,F),
           InstantPoorVAR = if_else(Instantaneous_sd<0.5*Observed_sd | 0.5*Instantaneous_sd>Observed_sd,T,F),
           HistoricalPoorVAR = if_else(Historical_sd<0.5*Observed_sd | 0.5*Historical_sd>Observed_sd,T,F))
  
  # Calculate the daily statistics for NEE
  DailyAvg.NEE = NEE.df %>% mutate(MonDay = format(as.Date(Date), "%m-%d")) %>%
    group_by(MonDay) %>%
    summarise(Observed_mean = mean(Observed*NEE.df.scale+NEE.df.center,na.rm=T),
              Observed_sd = sd(Observed*NEE.df.scale+NEE.df.center,na.rm=T),
              Instantaneous_mean = mean(Instantaneous*NEE.df.scale+NEE.df.center,na.rm=T),
              Instantaneous_sd = sd(Instantaneous*NEE.df.scale+NEE.df.center,na.rm=T),
              Historical_mean = mean(Historical*NEE.df.scale+NEE.df.center,na.rm=T),
              Historical_sd = sd(Historical*NEE.df.scale+NEE.df.center,na.rm=T)) %>%
    mutate(Date = as.Date(paste0("2000-",MonDay),format = "%Y-%m-%d"),
           "Flux" = "NEE",
           Observed_sd = ifelse(is.na(Observed_sd), 0.1*Observed_mean, Observed_sd),
           Instantaneous_sd = ifelse(is.na(Instantaneous_sd), 0.1*Instantaneous_mean, Instantaneous_sd),
           Historical_sd = ifelse(is.na(Historical_sd), 0.1*Historical_mean, Historical_sd),
           InstantPoorMean = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean-Observed_sd 
                                       | Instantaneous_mean-Instantaneous_sd > Observed_mean+Observed_sd),T,F),
           HistoricalPoorMean = if_else((Historical_mean+Historical_sd < Observed_mean-Observed_sd 
                                      | Historical_mean-Historical_sd > Observed_mean+Observed_sd),T,F),
           InstantPoorMean_percent = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean*0.9 
                                      | Instantaneous_mean-Instantaneous_sd > Observed_mean*1.1),T,F),
           HistoricalPoorMean_percent = if_else((Historical_mean < Observed_mean*0.9
                                         | Historical_mean-Historical_sd > Observed_mean*1.1),T,F),
           InstantPoorVAR = if_else(Instantaneous_sd<0.5*Observed_sd | 0.5*Instantaneous_sd>Observed_sd,T,F),
           HistoricalPoorVAR = if_else(Historical_sd<0.5*Observed_sd | 0.5*Historical_sd>Observed_sd,T,F))
  
  
  # Calculate the monthly statistics for LE
  MonthAvg.Qle = Qle.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Observed_month_mean = mean(Observed*Qle.df.scale+Qle.df.center),
              Instantaneous_month_mean = mean(Instantaneous*Qle.df.scale+Qle.df.center),
              Historical_month_mean = mean(Historical*Qle.df.scale+Qle.df.center)) %>%
    mutate(Mon = substr(YearMon,6,7)) %>%
    group_by(Mon) %>%
    summarise(Observed_mean = mean(Observed_month_mean),
              Observed_sd = sd(Observed_month_mean),
              Instantaneous_mean = mean(Instantaneous_month_mean),
              Instantaneous_sd = sd(Instantaneous_month_mean),
              Historical_mean = mean(Historical_month_mean),
              Historical_sd = sd(Historical_month_mean)) %>%
    mutate(Date = as.Date(paste0("2000-",Mon,"-01",format = "%Y-%m-%d")),
           "Flux" = "Qle",
           Observed_sd = ifelse(is.na(Observed_sd), 0.1*Observed_mean, Observed_sd),
           Instantaneous_sd = ifelse(is.na(Instantaneous_sd), 0.1*Instantaneous_mean, Instantaneous_sd),
           Historical_sd = ifelse(is.na(Historical_sd), 0.1*Historical_mean, Historical_sd),
           InstantPoorMean = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean-Observed_sd 
                                      | Instantaneous_mean-Instantaneous_sd > Observed_mean+Observed_sd),T,F),
           HistoricalPoorMean = if_else((Historical_mean+Historical_sd < Observed_mean-Observed_sd 
                                         | Historical_mean-Historical_sd > Observed_mean+Observed_sd),T,F),
           InstantPoorMean_percent = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean*0.9 
                                              | Instantaneous_mean-Instantaneous_sd > Observed_mean*1.1),T,F),
           HistoricalPoorMean_percent = if_else((Historical_mean < Observed_mean*0.9
                                                 | Historical_mean-Historical_sd > Observed_mean*1.1),T,F),
           InstantPoorVAR = if_else(Instantaneous_sd<0.5*Observed_sd | 0.5*Instantaneous_sd>Observed_sd,T,F),
           HistoricalPoorVAR = if_else(Historical_sd<0.5*Observed_sd | 0.5*Historical_sd>Observed_sd,T,F))
  
  # Calculate the daily statistics for LE
  DailyAvg.Qle = Qle.df %>% mutate(MonDay = format(as.Date(Date), "%m-%d")) %>%
    group_by(MonDay) %>%
    summarise(Observed_mean = mean(Observed*Qle.df.scale+Qle.df.center,na.rm=T),
              Observed_sd = sd(Observed*Qle.df.scale+Qle.df.center,na.rm=T),
              Instantaneous_mean = mean(Instantaneous*Qle.df.scale+Qle.df.center,na.rm=T),
              Instantaneous_sd = sd(Instantaneous*Qle.df.scale+Qle.df.center,na.rm=T),
              Historical_mean = mean(Historical*Qle.df.scale+Qle.df.center,na.rm=T),
              Historical_sd = sd(Historical*Qle.df.scale+Qle.df.center,na.rm=T)) %>%
    mutate(Date = as.Date(paste0("2000-",MonDay),format = "%Y-%m-%d"),
           "Flux" = "Qle",
           Observed_sd = ifelse(is.na(Observed_sd), 0.1*Observed_mean, Observed_sd),
           Instantaneous_sd = ifelse(is.na(Instantaneous_sd), 0.1*Instantaneous_mean, Instantaneous_sd),
           Historical_sd = ifelse(is.na(Historical_sd), 0.1*Historical_mean, Historical_sd),
           InstantPoorMean = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean-Observed_sd 
                                      | Instantaneous_mean-Instantaneous_sd > Observed_mean+Observed_sd),T,F),
           HistoricalPoorMean = if_else((Historical_mean+Historical_sd < Observed_mean-Observed_sd 
                                         | Historical_mean-Historical_sd > Observed_mean+Observed_sd),T,F),
           InstantPoorMean_percent = if_else((Instantaneous_mean+Instantaneous_sd < Observed_mean*0.9 
                                              | Instantaneous_mean-Instantaneous_sd > Observed_mean*1.1),T,F),
           HistoricalPoorMean_percent = if_else((Historical_mean < Observed_mean*0.9
                                                 | Historical_mean-Historical_sd > Observed_mean*1.1),T,F),
           InstantPoorVAR = if_else(Instantaneous_sd<0.5*Observed_sd | 0.5*Instantaneous_sd>Observed_sd,T,F),
           HistoricalPoorVAR = if_else(Historical_sd<0.5*Observed_sd | 0.5*Historical_sd>Observed_sd,T,F))

  Monthly.df = rbind(MonthAvg.Qle,MonthAvg.NEE)
  Daily.df = rbind(DailyAvg.NEE,DailyAvg.Qle)
  
  df$Site = Site
  Monthly.df$Site = Site
  Daily.df$Site = Site
  
  All.All = rbind(All.All,df)
  All.Monthly = rbind(All.Monthly,Monthly.df)
  All.Daily = rbind(All.Daily,Daily.df)
  
  
 # Calculate timeseries metrics
 # Monthly NEE
 # Calculate monthly timeseries
  MonthFit.NEE = NEE.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Observed_month_mean = mean(Observed*Qle.df.scale+Qle.df.center),
              Instantaneous_month_mean = mean(Instantaneous*Qle.df.scale+Qle.df.center),
              Historical_month_mean = mean(Historical*Qle.df.scale+Qle.df.center))
  # Calculate instantaneous monthly metrics
  NME = sum(abs(MonthFit.NEE$Instantaneous_month_mean-MonthFit.NEE$Observed_month_mean))/sum(abs(mean(MonthFit.NEE$Observed_month_mean)-MonthFit.NEE$Observed_month_mean))
  MBE = sum(MonthFit.NEE$Instantaneous_month_mean-MonthFit.NEE$Observed_month_mean,na.rm=TRUE)/length(MonthFit.NEE$Instantaneous_month_mean)
  CCO = cor(MonthFit.NEE$Instantaneous_month_mean,MonthFit.NEE$Observed_month_mean,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(MonthFit.NEE$Instantaneous_month_mean,na.rm=TRUE)/sd(MonthFit.NEE$Observed_month_mean,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "NEE", "Time" = "Monthly", "Model" = "Instantaneous", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate daily instantaneous metrics
  NME = sum(abs(NEE.df$Instantaneous-NEE.df$Observed))/sum(abs(mean(NEE.df$Observed)-NEE.df$Observed))
  MBE = sum(NEE.df$Instantaneous-NEE.df$Observed,na.rm=TRUE)/length(NEE.df$Instantaneous)
  CCO = cor(NEE.df$Instantaneous,NEE.df$Observed,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(NEE.df$Instantaneous,na.rm=TRUE)/sd(NEE.df$Observed,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "NEE", "Time" = "Daily", "Model" = "Instantaneous", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate monthly historical metrics
  NME = sum(abs(MonthFit.NEE$Historical_month_mean-MonthFit.NEE$Observed_month_mean))/sum(abs(mean(MonthFit.NEE$Observed_month_mean)-MonthFit.NEE$Observed_month_mean))
  MBE = sum(MonthFit.NEE$Historical_month_mean-MonthFit.NEE$Observed_month_mean,na.rm=TRUE)/length(MonthFit.NEE$Historical_month_mean)
  CCO = cor(MonthFit.NEE$Historical_month_mean,MonthFit.NEE$Observed_month_mean,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(MonthFit.NEE$Historical_month_mean,na.rm=TRUE)/sd(MonthFit.NEE$Observed_month_mean,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "NEE", "Time" = "Monthly", "Model" = "Historical", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate daily historical metrics
  NME = sum(abs(NEE.df$Historical-NEE.df$Observed))/sum(abs(mean(NEE.df$Observed)-NEE.df$Observed))
  MBE = sum(NEE.df$Historical-NEE.df$Observed,na.rm=TRUE)/length(NEE.df$Historical)
  CCO = cor(NEE.df$Historical,NEE.df$Observed,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(NEE.df$Historical,na.rm=TRUE)/sd(NEE.df$Observed,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "NEE", "Time" = "Daily", "Model" = "Historical", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Monthly Qle
  MonthFit.Qle = Qle.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m")) %>%
    group_by(YearMon) %>%
    summarise(Observed_month_mean = mean(Observed*Qle.df.scale+Qle.df.center),
              Instantaneous_month_mean = mean(Instantaneous*Qle.df.scale+Qle.df.center),
              Historical_month_mean = mean(Historical*Qle.df.scale+Qle.df.center)) 
  # Calculate monthly instantaneous metrics
  NME = sum(abs(MonthFit.Qle$Instantaneous_month_mean-MonthFit.Qle$Observed_month_mean))/sum(abs(mean(MonthFit.Qle$Observed_month_mean)-MonthFit.Qle$Observed_month_mean))
  MBE = sum(MonthFit.Qle$Instantaneous_month_mean-MonthFit.Qle$Observed_month_mean,na.rm=TRUE)/length(MonthFit.Qle$Instantaneous_month_mean)
  CCO = cor(MonthFit.Qle$Instantaneous_month_mean,MonthFit.Qle$Observed_month_mean,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(MonthFit.Qle$Instantaneous_month_mean,na.rm=TRUE)/sd(MonthFit.Qle$Observed_month_mean,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "Qle", "Time" = "Monthly", "Model" = "Instantaneous", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate instantaneous daily metrics
  NME = sum(abs(Qle.df$Instantaneous-Qle.df$Observed))/sum(abs(mean(Qle.df$Observed)-Qle.df$Observed))
  MBE = sum(Qle.df$Instantaneous-Qle.df$Observed,na.rm=TRUE)/length(Qle.df$Instantaneous)
  CCO = cor(Qle.df$Instantaneous,Qle.df$Observed,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(Qle.df$Instantaneous,na.rm=TRUE)/sd(Qle.df$Observed,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "Qle", "Time" = "Daily", "Model" = "Instantaneous", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate monthly historical metrics
  NME = sum(abs(MonthFit.Qle$Historical_month_mean-MonthFit.Qle$Observed_month_mean))/sum(abs(mean(MonthFit.Qle$Observed_month_mean)-MonthFit.Qle$Observed_month_mean))
  MBE = sum(MonthFit.Qle$Historical_month_mean-MonthFit.Qle$Observed_month_mean,na.rm=TRUE)/length(MonthFit.Qle$Historical_month_mean)
  CCO = cor(MonthFit.Qle$Historical_month_mean,MonthFit.Qle$Observed_month_mean,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(MonthFit.Qle$Historical_month_mean,na.rm=TRUE)/sd(MonthFit.Qle$Observed_month_mean,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "Qle", "Time" = "Monthly", "Model" = "Historical", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
  # Calculate daily historical metrics
  NME = sum(abs(Qle.df$Historical-Qle.df$Observed))/sum(abs(mean(Qle.df$Observed)-Qle.df$Observed))
  MBE = sum(Qle.df$Historical-Qle.df$Observed,na.rm=TRUE)/length(Qle.df$Historical)
  CCO = cor(Qle.df$Historical,Qle.df$Observed,use = "complete.obs", method = "pearson")
  SDD = abs(1-sd(Qle.df$Historical,na.rm=TRUE)/sd(Qle.df$Observed,na.rm=TRUE))
  # Add to df
  metricrow = data.frame("Site" = Site, "Flux" = "Qle", "Time" = "Daily", "Model" = "Historical", "NME" = NME, "MBE" = MBE, "CCO" = CCO, "SDD" = SDD)
  Metrics = rbind(Metrics,metricrow)
  
}



MonthlySummary = All.Monthly %>% group_by(Site,Flux) %>%
  summarise(InstantPoorMean = if_else(sum(InstantPoorMean)>0,T,F),
            HistoricalPoorMean = if_else(sum(HistoricalPoorMean)>0,T,F),
            InstantPoorVAR = if_else(sum(InstantPoorVAR)>=3,T,F),
            HistoricalPoorVAR = if_else(sum(HistoricalPoorVAR)>=3,T,F))

DailySummary = All.Daily %>% group_by(Site,Flux) %>%
  summarise(InstantPoorMean = if_else(sum(InstantPoorMean)>0,T,F),
            HistoricalPoorMean = if_else(sum(HistoricalPoorMean)>0,T,F),
            InstantPoorMean_percent = if_else(sum(InstantPoorMean_percent)>30,T,F),
            HistoricalPoorMean_percent = if_else(sum(HistoricalPoorMean_percent)>30,T,F),
            InstantPoorVAR = if_else(sum(InstantPoorVAR)>=91,T,F),
            HistoricalPoorVAR = if_else(sum(HistoricalPoorVAR)>=91,T,F))


# Find the number of sites for each of the categories
sum(DailySummary$InstantPoorMean[DailySummary$Flux=="NEE"])/65
sum(DailySummary$HistoricalPoorMean[DailySummary$Flux=="NEE"])/65
sum(DailySummary$InstantPoorVAR[DailySummary$Flux=="NEE" & DailySummary$InstantPoorMean == F])/65
sum(DailySummary$HistoricalPoorVAR[DailySummary$Flux=="NEE" & DailySummary$HistoricalPoorMean == F])/65

sum(DailySummary$InstantPoorMean[DailySummary$Flux=="Qle"])/65
sum(DailySummary$HistoricalPoorMean[DailySummary$Flux=="Qle"])/65
sum(DailySummary$InstantPoorVAR[DailySummary$Flux=="Qle" & DailySummary$InstantPoorMean == F])/65
sum(DailySummary$HistoricalPoorVAR[DailySummary$Flux=="Qle" & DailySummary$HistoricalPoorMean == F])/65



sum(MonthlySummary$InstantPoorMean[MonthlySummary$Flux=="NEE"])/65
sum(MonthlySummary$HistoricalPoorMean[MonthlySummary$Flux=="NEE"])/65
sum(MonthlySummary$InstantPoorVAR[MonthlySummary$Flux=="NEE" & MonthlySummary$InstantPoorMean == F])/65
sum(MonthlySummary$HistoricalPoorVAR[MonthlySummary$Flux=="NEE" & MonthlySummary$HistoricalPoorMean == F])/65

sum(MonthlySummary$InstantPoorMean[MonthlySummary$Flux=="Qle"])/65
sum(MonthlySummary$HistoricalPoorMean[MonthlySummary$Flux=="Qle"])/65
sum(MonthlySummary$InstantPoorVAR[MonthlySummary$Flux=="Qle" & MonthlySummary$InstantPoorMean == F])/65
sum(MonthlySummary$HistoricalPoorVAR[MonthlySummary$Flux=="Qle" & MonthlySummary$HistoricalPoorMean == F])/65


Metrics %>% group_by(Time,Flux,Model) %>%
  summarise(NME_mean = mean(NME),
            NME_median = median(NME),
            MBE_mean = mean(MBE),
            MBE_median = median(MBE),
            CCO_mean = mean(CCO),
            CCO_median = median(CCO),
            SDD_mean = mean(SDD),
            SDD_median = median(SDD))



All.Daily %>% group_by(Site,Flux) %>%
  summarise(InstantPoorMean = sum(InstantPoorMean),
            HistoricalPoorMean = sum(HistoricalPoorMean),
            InstantPoorMean_percent = sum(InstantPoorMean_percent),
            HistoricalPoorMean_percent = sum(HistoricalPoorMean_percent),
            InstantPoorVAR = sum(InstantPoorVAR),
            HistoricalPoorVAR = sum(HistoricalPoorVAR))
