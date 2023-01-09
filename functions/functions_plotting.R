################################################################################

### Plotting functions for use with k-means clustering

################################################################################

# Libraries
library(ggplot2)
library(ggridges)
library(scales)
library(ggnewscale)

################################################################################
InputPlot = function(Site,Data){
  # Function to plot the inputs (met and flux) from the processed FLUXNET files
  # Plot is faceted and allows a quick check that input timeseries are sensible
  # Inputs:
  #       Site = site name
  #       Data = dataframe from output of FluxnetProcess function
  
  
  df = Data %>%
    select("Date",
           "NEE",
           "Qle",
           "SWdown",
           "Tair",
           "Wind",
           "VPD",
           "Precip",# Select long-term rainfall
           "LAI_MOD") %>%
    pivot_longer(NEE:LAI_MOD)
  df$name = factor(df$name,levels = unique(df$name))
  
  Plot = df %>% ggplot() +
    geom_line(aes(x=Date,y=value)) +
    facet_wrap(.~name, scales="free") +
    theme_bw() +
    ggtitle(Site)
  
  ggsave(filename = paste0("inputs/images/",Site,".png"), 
         plot = Plot,
         device = "png",
         dpi = 320,
         width = 10, height = 5)
}

################################################################################
plotobsvpred = function(Site,k,df){
  # Function to plot the timeseries of observed NEE against modelled NEE
  # Inputs:
  #       Site = name of Site
  #       k = number of clusters
  #       df = dataframe of flux, with observations in column 2 and modelled
  #            values in column 3
  
  Plot = df %>% pivot_longer(cols=2:3) %>%
    ggplot() +
    geom_line(aes(x=Date,y=value,colour=name),size=1) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(df)$Date[1],tail(na.omit(df)$Date,n=1))) +
    
    theme_bw() +
    xlab("") +
    ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(Site," Model - ",k," clusters"))
}

################################################################################
plotcoeffs_Abs = function(Site,k,coeff_df,Data){
  # Function to plot timeseries of the absolute coefficient values from the 
  # regressions
  # Inputs:
  #       Site = name of site
  #       k = number of clusters
  #       coeff_df = pivotted dataframe of coefficient values with coefficient
  #                 names in column "name" and values (non-absolute) in column
  #                 "value"
  
  Plot = coeff_df %>%
    ggplot() +
    geom_linerange(aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       colour=abs(value)),
                   size=10) +
    binned_scale(aesthetics = "color",
                 scale_name = "stepsn",
                 name = "Coefficient\nValue",
                 palette = function(x) c("white",'#d1e5f0','#67a9cf','#2166ac'),
                 breaks = c(0.1, 1,5),
                 oob = scales::squish,
                 limits = c(0, 10),
                 show.limits = TRUE, 
                 guide = "colorsteps")+
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],
                            tail(na.omit(coeff_df)$Date,n=1)))  +
    theme_bw() +
    coord_flip() +
    xlab("Driver") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    ggtitle(paste0(Site," Abs Lag Weights - ",k," clusters - Filtered"))
}

################################################################################
plotcoeffs = function(Site,k,coeff_df,Data){
  # Function to plot timeseries of the coefficient values from the regressions
  # Inputs:
  #       Site = name of site
  #       k = number of clusters
  #       coeff_df = pivotted dataframe of coefficient values with coefficient
  #                 names in column "name" and values in column "value"
  
  Plot = coeff_df %>%
    ggplot() +
    geom_linerange(aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       colour=value),
                   size=5) +
    binned_scale(aesthetics = "color",
                 scale_name = "stepsn",
                 name = "Coefficient\nValue",
                 palette = function(x) c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#ffffff','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419'),
                 breaks = c(-1000,-100,-10,-1, -0.1, 0.1, 1,10,100,1000),
                 labels = function(x) as.character(x),
                 oob = scales::squish,
                 limits = c(-100000,100000),
                #show.limits = TRUE,
                 guide = "colorsteps")+
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],
                            tail(na.omit(coeff_df)$Date,n=1))) +
    theme_bw() +
    coord_flip() +
    xlab("Driver") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    ggtitle(paste0(Site," Lag Weights - ",k," clusters - Filtered"))
}

################################################################################
plotheatmap = function(Site,flux,k,km,name,include_wind=F){
  # Function to plot a heatmap of cluster centers to see the regimes
  # Inputs:
  #       Site = site name
  #       k = number of clusters
  #       km = output of KMeans_rcpp function
  #       name = name of model runs

  # Create dataframe of cluster centres
  cluster <- c(1:k)
  center <- km$centroids 
  colnames(center) = c("LAI",
                       "SWdown",
                       "Tair",
                       c("Wind")[include_wind],
                       "Precip")
  center_df <- data.frame(cluster, center)
  # Order so that more extreme cluster centers appear at the top of the plot
  center_df_sort = order(apply(center_df[,2:(5+include_wind)],1,function(x) max(abs(x))))
  
  # Reshape the data
  center_reshape <- gather(center_df, features, values, LAI:Precip)
  center_reshape$features = factor(center_reshape$features,levels = unique(center_reshape$features))
  center_reshape$cluster = factor(center_reshape$cluster,levels = center_df_sort)
  # Plot the heat map of cluster centres
  limits = max(abs(center_reshape$values)) * c(-1, 1)
  heatmap = ggplot(data = center_reshape, aes(x = features, y = cluster, fill = values)) +
    scale_y_discrete(breaks = seq(1, k, by = 1)) +
    geom_tile() +
    coord_equal() +
    scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=-90))
  
  # Save the heatmap
  ggsave(filename = paste0("outputs/kmeans/",flux,"/",name,"/images/heatmaps/",Site,"_",k,"c.png"), 
         plot = heatmap,
         width = 5, height = 49)
}

################################################################################
plotlinear = function(Site){
  
  # A function to plot the coefficients for a linear regression over the entire
  # time series (a.k.a. a single cluster). Each variable is coloured by the 
  # coefficient value and labeled with the p-value
  # Inputs:
  #       Site = site name
  
  # Load the input file and extract required data
  load(paste0("inputs/FLUXNET_processed/",Site,".Rdata"))
  
  # Manipulate the rainfall data
  Climate = Data %>%
    na.omit() %>% # Remove na rows from start
    select("SWdown",
           "Tair",
           "Precip",# Select long-term rainfall
           "P1to30",
           "P31to90",
           "P91to180",
           "P181to365",
           "P366to730",
           "P731to1095",
           "P1096to1460") %>%
    scale() # Scale the data
  
  # Initialise the comparison dataframe
  NEE_df = data.frame("Date" = Data %>% na.omit() %>% select(Date),
                      "NEE" = Data %>% na.omit() %>% select(NEE) %>% scale())
  
  model = lm(NEE_df$NEE~Climate)
  summary = summary(model)
  r2 = signif(summary$adj.r.squared,2)
  
  pvals = data.frame("name"=row.names(summary$coefficients),
                     "pval"=summary$coefficients[,4]) %>% remove_rownames()
  
  coeff_df = data.frame("Date"=Data %>%
                          na.omit() %>%
                          select(Date))
  
  coeff = model$coefficients # Extract coefficients
  coeff_df[,names(coeff)] = coeff %>% unname() %>%
    rep(.,each=nrow(coeff_df))
  
  coeff_df = coeff_df %>% pivot_longer(cols="ClimateSWdown":"ClimateP1096to1460") %>% merge(pvals) %>% arrange(Date)
  
  
  coeff_df$name = coeff_df$name %>% str_replace(pattern = "Climate",
                                                replacement = "") %>%
    factor(levels=c("SWdown",
                    "Tair",
                    "Precip",
                    "P1to30",
                    "P31to90",
                    "P91to180",
                    "P181to365",
                    "P366to730",
                    "P731to1095",
                    "P1096to1460"))
  
  MidDate = coeff_df$Date[ceiling(length(coeff_df$Date)/2)]
  
  Plot = ggplot() +
    geom_linerange(data=coeff_df,aes(x=name,ymin=Date,ymax=Date+1,colour=value),size=10) +
    geom_label(data=coeff_df[1:10,],aes(x=name,y=rep(MidDate,10),label=round(pval,2))) +
    scale_colour_gradient2(low="red",mid="white",high="blue",midpoint=0) +
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],tail(na.omit(coeff_df)$Date,n=1))) +
    
    theme_bw() +
    coord_flip() +
    xlab("Driver") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(Site," - Linear Regression - R^2 = ",r2))
}


################################################################################
plotmetricbarplot = function(sites,k,flux,include_wind=F){
  
  # A function to plot boxplots of model metrics for (hard-coded) different lags
  # to ensure that models are improving the fit as expected
  # Inputs:
  #       sites = list of site names
  #       k = number of clusters
  #       flux = the name of the modelled flux
  #       include_wind = boolean, models including wind or leaving it out
  
  # Create dataframe to place metric values in
  df = data.frame("Site" = rep(sites,
                               each=15),
                  "Model" = rep(c("NoLag","LagOnly","4yrLag"),
                                each = 5,
                                times = length(sites)),
                  "Metric" = rep(c("R^2",
                                   "Mean Bias Error",
                                   "Normalised Mean Error",
                                   "Std. Dev. Difference",
                                   "Correlation Coeff."),
                                 times = 3*length(sites)),
                  "Value" = NA)
  
  # Load in the metric values
  for (Site in sites){
    
    load(paste0("outputs/kmeans/",flux,"/noPPTlag",c("_W/")[include_wind],Site,"_",k,"c.Rdata"))
    df$Value[df$Site == Site & df$Model == "NoLag"] = unlist(output[c(3:7)])
    rm(output)
    
    load(paste0("outputs/kmeans/",flux,"/4yearPPTlag",c("_W/")[include_wind],Site,"_",k,"c.Rdata"))
    df$Value[df$Site == Site & df$Model == "4yrLag"] = unlist(output[c(3:7)])
    rm(output)
    
    load(paste0("outputs/kmeans/",flux,"/only4yearPPTlag",c("_W/")[include_wind],Site,"_",k,"c.Rdata"))
    df$Value[df$Site == Site & df$Model == "LagOnly"] = unlist(output[c(3:7)])
    rm(output)
    
    
  }
  
  # Turn model and metrics into factors so they appear in correct order
  df$Model = factor(df$Model,
                    levels = c("NoLag","4yrLag","LagOnly"))
  
  
  df$Metric = factor(df$Metric, 
                     levels = c("Normalised Mean Error",
                                "Std. Dev. Difference",
                                "Mean Bias Error",
                                "R^2",
                                "Correlation Coeff."),
                     labels = c(expression("Normalised~Mean~Error"),
                                expression("Std.~Dev.~Difference"),
                                expression("Mean~Bias~Error"),
                                expression("R^2"),
                                expression("Correlation~Coeff.")))
  
  
  
  # Plot the metric values
  Plot = ggplot(df) +
    geom_boxplot(aes(x=Model,y=Value,fill=Model)) +
    facet_grid(.~Metric,
               scales="free",
               labeller = label_parsed) +
    theme_bw() +
    xlab("") +
    scale_fill_viridis_d(direction = 1, end = 0.9, begin = 0.4) +
    scale_y_continuous(expand=c(0.05,0)) +
    theme(panel.grid.major.x = element_blank(),
          text = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor.y = element_blank()) +
    ggtitle(paste0(flux," Metrics - ",k,"c"))

}

################################################################################
plotmetriclineplot = function(sites,k,flux){
  
  # A function to plot line and point plots of model metrics for (hard-coded) 
  # different lags to ensure that models are improving the fit as expected
  # Inputs:
  #       sites = list of site names
  #       k = number of clusters
  #       flux = the name of the modelled flux
  
  # Create dataframe to place metric values in
  df = data.frame("Site" = rep(sites,
                               each=30),
                  "Cluster" = rep(c("Met","Met+Lag"),
                                  each = 15,
                                  times = length(sites)),
                  "Regression" = rep(c("Met","Lag","Met+Lag"),
                                     each = 5,
                                     times = 2*length(sites)),
                  "Metric" = rep(c("R^2",
                                   "Mean Bias Error",
                                   "Normalised Mean Error",
                                   "Std. Dev. Difference",
                                   "Correlation Coeff."),
                                 times = 6*length(sites)),
                  "Value" = NA)
  
  # Load in the metric values
  for (Site in sites){
    
    # Met cluster and met regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metcluster_metregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met" & 
               df$Regression == "Met"] = unlist(output[c(3:7)])
    rm(output)
    
    # Met cluster and lag regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metcluster_lagregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met" & 
               df$Regression == "Lag"] = unlist(output[c(3:7)])
    rm(output)
    
    # Met cluster and metlag regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metcluster_metlagregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met" & 
               df$Regression == "Met+Lag"] = unlist(output[c(3:7)])
    rm(output)
    
    # Metlag cluster and met regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metlagcluster_metregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met+Lag" & 
               df$Regression == "Met"] = unlist(output[c(3:7)])
    rm(output)
    
    # Metlag cluster and lag regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metlagcluster_lagregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met+Lag" & 
               df$Regression == "Lag"] = unlist(output[c(3:7)])
    rm(output)
    
    # Metlag cluster and metlag regression
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metlagcluster_metlagregress_",
                k,
                "c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Cluster == "Met+Lag" & 
               df$Regression == "Met+Lag"] = unlist(output[c(3:7)])
    rm(output)
    
    
  }
  
  # Turn model and metrics into factors so they appear in correct order
  df$Cluster = factor(df$Cluster,
                      levels = c("Met","Met+Lag"))
  
  df$Regression = factor(df$Regression,
                         levels = c("Met","Met+Lag","Lag"))
  
  df$Metric = factor(df$Metric, 
                     levels = c("Normalised Mean Error",
                                "Std. Dev. Difference",
                                "Mean Bias Error",
                                "R^2",
                                "Correlation Coeff."),
                     labels = c(expression("Normalised~Mean~Error"),
                                expression("Std.~Dev.~Difference"),
                                expression("Mean~Bias~Error"),
                                expression("R^2"),
                                expression("Correlation~Coeff.")))
  
  # Plot the metric values
  Plot =  df %>% #filter(Cluster == "Met+Lag" & Metric == "R^2") %>%
    ggplot() +
    geom_line(aes(x=Regression,y=Value,group=Site,color=Regression)) +
    geom_point(aes(x=Regression,y=Value,fill=Regression,group=Site),size = 3, shape = 21) + 
    facet_grid(Cluster~Metric,
               scales="free",
               labeller = label_parsed) +
    theme_bw() +
    xlab("") +
    scale_color_viridis_d(name = "Regression",direction = 1, end = 0.9, begin = 0.4) +
    scale_fill_viridis_d(direction = 1, end = 0.9, begin = 0.4) +
    scale_y_continuous(expand=c(0.05,0)) +
    theme(panel.grid.major.x = element_blank(),
          text = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    ggtitle(paste0(flux," Metrics - ",k,"c"))
  
}