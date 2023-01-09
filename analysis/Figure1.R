### Workflow to plot timeseries for 4 specific sites

### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")
library(zoo)
library(grid)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(viridisLite)

### End Setup ###
###############################################################################
###############################################################################

# Load sites
sites = c("CH-Dav","DE-Kli","DK-ZaH","US-SRM")

###############################################################################
###############################################################################

# Define the models used
models = c("Observed",
           "Instantaneous",
           "Historical")
labels = c("Observed Flux",
           "Historical Model",
           "Instantaneous Model")

for (Site in sites){
  
  ### NEP
  # Load the instantaneous model output
  load(paste0("outputs/kmeans/NEE/",
              Site,
              "_metcluster_metregress_729c.Rdata"))
  
  instantaneous = output$Flux_df %>% select(Date,Flux_pred)
  names(instantaneous) = c("Date","Instantaneous")
  
  # Load the historical model output
  load(paste0("outputs/kmeans/NEE/",          
              Site,
              "_metlagcluster_metlagregress_729c.Rdata"))
  
  historical = output$Flux_df %>% select(Date,Flux,Flux_pred)
  names(historical) = c("Date","Observed","Historical")
  
  # Find how the flux was scaled prior to modelling
  NEE.df.center = attr(historical$Observed,c("scaled:center"))
  NEE.df.scale = attr(historical$Observed,c("scaled:scale"))
  
  NEE.df = merge(instantaneous,historical,by.x = "Date", by.y = "Date")
  NEE.df$Flux = "NEE"
  # Invert to make NEE into NEP
  NEE.df[,models] = -NEE.df[,models]
  
  ### Latent Heat
  # Load the instantaneous model output
  load(paste0("outputs/kmeans/Qle/",
              Site,
              "_metcluster_metregress_729c.Rdata"))
  
  instantaneous = output$Flux_df %>% select(Date,Flux_pred)
  names(instantaneous) = c("Date","Instantaneous")
  
  # Load the historical model output
  load(paste0("outputs/kmeans/Qle/",          
              Site,
              "_metlagcluster_metlagregress_729c.Rdata"))
  
  historical = output$Flux_df %>% select(Date,Flux,Flux_pred)
  names(historical) = c("Date","Observed","Historical")
  
  # Find how the flux was scaled prior to modelling
  Qle.df.center = attr(historical$Observed,c("scaled:center"))
  Qle.df.scale = attr(historical$Observed,c("scaled:scale"))
  
  Qle.df = merge(instantaneous,historical,by.x = "Date", by.y = "Date")
  Qle.df$Flux = "Qle"
  
  # Combine dataframes
  df = rbind(NEE.df,Qle.df)
  
  # Plot for NEP
  NEP.Plot = NEE.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m"),
                               Observed = Observed*NEE.df.scale+NEE.df.center,
                               Instantaneous = Instantaneous*NEE.df.scale+NEE.df.center,
                               Historical = Historical*NEE.df.scale+NEE.df.center) %>%
    group_by(YearMon) %>%
    summarise_at(c("Observed","Instantaneous","Historical"),
                 list(min = min, mean = mean,max = max,sd = sd)) %>%
    mutate(Date = as.Date(paste0(YearMon,"-01",format = "%Y-%m-%d"))) %>%
    pivot_longer(cols = 2:13, names_to = c("Model",".value"),names_sep = "_") %>%
    mutate(Model = factor(Model, levels = models)) %>%
    ggplot() +
    geom_ribbon(aes(x=Date, ymin = (mean-sd)/1000, ymax = (mean+sd)/1000,fill = Model),
                alpha = 0.2) +
    geom_line(aes(x = Date, y = mean/1000,color = Model),
              alpha = 0.9, 
              size = 1.5) +
    scale_color_viridis_d(option = "inferno", end = 0.9) +
    scale_fill_viridis_d(option = "inferno", end = 0.9) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(NEE.df)$Date[1],
                            tail(na.omit(NEE.df)$Date,n=1))) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(text = element_text(size=20),
          legend.position = "none") +
    ylab(expression("NEP"~"(mmol"~"/"~"m"^"2"*"day)")) +
    xlab("")
  
  name = paste0(Site,".NEP.Plot")
  assign(name,NEP.Plot)

  LE.Plot = Qle.df %>% mutate(YearMon = format(as.Date(Date), "%Y-%m"),
                              Observed = Observed*Qle.df.scale+Qle.df.center,
                              Instantaneous = Instantaneous*Qle.df.scale+Qle.df.center,
                              Historical = Historical*Qle.df.scale+Qle.df.center) %>%
    group_by(YearMon) %>%
    summarise_at(c("Observed","Instantaneous","Historical"),
                 list(min = min, mean = mean,max = max,sd = sd)) %>%
    mutate(Date = as.Date(paste0(YearMon,"-01",format = "%Y-%m-%d"))) %>%
    pivot_longer(cols = 2:13, names_to = c("Model",".value"),names_sep = "_") %>%
    mutate(Model = factor(Model,levels = models)) %>%
    ggplot() +
    geom_ribbon(aes(x=Date, ymin = (mean-sd)/86400, ymax = (mean+sd)/86400,fill = Model),
                alpha = 0.2) +
    geom_line(aes(x = Date, y = mean/86400,color = Model),
              size = 1.5, 
              alpha = 0.9) +
    scale_color_viridis_d(option = "inferno", end = 0.9) +
    scale_fill_viridis_d(option = "inferno", end = 0.9) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(Qle.df)$Date[1],
                            tail(na.omit(Qle.df)$Date,n=1))) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(text = element_text(size=20),
          legend.position = "none") +
    ylab(expression("LE"~"(W"~"/"~"m"^"2"*")")) +
    xlab("")

  name = paste0(Site,".LE.Plot")
  assign(name,LE.Plot)
}

# Create a plot to extract the legend
Legend.Plot = df %>% select(Date,Flux,Observed,Instantaneous,Historical) %>%
  pivot_longer(names_to = "Model",cols = 4:5, values_to = "Predicted") %>%
  ggplot() +
  geom_smooth(aes(x = Observed, y = Predicted, color = Model, fill = Model), 
              method = "lm", 
              se = T, 
              fullrange = T, 
              size = 3, 
              alpha = 0.2) +
  geom_point(aes(x = 1, y = 1, fill = "Bluh", color = "Bluh")) +
  scale_color_viridis_d("",
                        option = "inferno", 
                        end = 0.9, 
                        labels = labels) +
  scale_fill_viridis_d("",
                       option = "inferno",
                       end = 0.9, 
                       labels = labels) + 
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = "top",       
        legend.title.align=0.5,
        legend.key.size= unit(1, 'cm')) +
  guides(fill = guide_legend(nrow = 1,title.position = "top",reverse=F),
         color = guide_legend(nrow = 1,title.position = "top",reverse=F))

legend = get_legend(Legend.Plot)

# Save plot
png(filename = "images/Figure1.png",
    width = 16,
    height = 16,
    units = "in",
    res = 320)

  print(
        ggarrange(`CH-Dav.NEP.Plot`,
                  `DE-Kli.LE.Plot`,
                  `DK-ZaH.NEP.Plot`,
                  `US-SRM.LE.Plot`,
                  legend,
                  labels = c("(a)","(b)","(c)","(d)"), 
                  font.label = list(size = 20),
                  nrow = 5,
                  hjust = -0.1,
                  vjust = 1.1, 
                  heights = c(3,3,3,3,1)))
dev.off()



