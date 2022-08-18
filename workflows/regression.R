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
# 

flux = "NEE"
sitecsv = "longsites.csv"
seq = c(729)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load sites
sites = read_csv(paste0("inputs/",
                        sitecsv),
                 col_names=FALSE,
                 show_col_types = F) %>% 
        unlist()

g <- expand.grid(sites=sites,
                 k=seq,
                 clustmet = T,
                 clustlag = c(T,F),
                 regmet = c(T,F),
                 reglag = c(T,F)) %>%
  filter(!(regmet == F & reglag == F))

worker = function(Site,
                  k,
                  cluster_include_met,
                  cluster_include_lags,
                  regress_include_met,
                  regress_include_lags){
  
  # Load the input file and extract required data
  load(paste0("inputs/FLUXNET_processed/",
              Site,
              ".Rdata"))
  
  # Load the clustering output and extract required data
  load(paste0("outputs/kmeans/clusterings/",
              Site,
              "_",
              c("met")[cluster_include_met],
              c("lag")[cluster_include_lags],
              "_",
              k,
              "c.Rdata"))
  
  # Define climate variables
  metvars = c("LAI_MOD",
              "SWdown",
              "VPD",
              "Wind",
              "Tair",
              "Precip")
  
  lagvars = c("T1",
              "T2to7",
              "T8to14",
              "T15to30",
              "P1to30",
              "P31to90",
              "P91to180",
              "P181to365",
              "P366to730",
              "P731to1095",
              "P1096to1460")
  
  # Combine into regression predictors based on requirements
  climvars = c(metvars[regress_include_met],
               lagvars[regress_include_lags])
  
  # Select and scale the requested data
  predictor_data = Data %>% 
                    na.omit() %>% # Remove na rows from start
                    select(all_of(climvars)) %>%
                    scale() # Scale the data
  
  # Extract the clusters and cluster size
  clusters = km$clusters
  clustersize = km$obs_per_cluster
  
  # Initialise the comparison dataframe
  Flux_df = data.frame("Date" = Data %>% na.omit() %>% pull(Date),
                       "Flux" = Data %>% na.omit() %>% pull(flux),
                       "Flux_pred" = 0, 
                       "cluster" = clusters)
  
  # Scale the Flux
  Flux_df$Flux = Flux_df$Flux %>% scale()
  
  # Initialise output list
  output = list()
  info = list()
  R2 = data.frame("Cluster" = 1:k,
                  "R2" = 0)
  
  # For each cluster
  for (i in 1:k){
    # If cluster has no elements, save a blank list
    if (clustersize[i]==0){
        message("Cluster ",i," has no elements")
        info[[paste0("cluster_",i)]] = list()
      
    # If it is less than 10x the predictors, bring in the closest ones for the regression
    } else if (clustersize[i] < 10*ncol(predictor_data)){
        # Perform the regression
        clusterinfo = smallclusters(clusters,
                                    clustersize,
                                    predictor_data,
                                    km$cluster_data,
                                    Flux_df,
                                    i,
                                    km)
        # Extract the info from the regression
        Flux_df$Flux_pred[clusters==i] = clusterinfo$Flux_pred
        info[[paste0("cluster_",i)]] = clusterinfo
        R2$R2[R2$Cluster==i] = clusterinfo$expanded.r.squared
      
    # Otherwise regress as normal
    } else {
        # Perform the regression
        clusterinfo = normalclusters(clusters,
                                     predictor_data,
                                     Flux_df,
                                     i)
        # Extract the info from the regression
        Flux_df$Flux_pred[clusters==i] = clusterinfo$Flux_pred
        info[[paste0("cluster_",i)]] = clusterinfo
        R2$R2[R2$Cluster==i] = clusterinfo$r.squared
        
    }
  }
  
  # Define the observed and predicted Flux
  Flux_obs = Flux_df$Flux
  Flux_pred = Flux_df$Flux_pred
  output[["Flux_df"]] = Flux_df
  output[["predictor_data"]] = predictor_data
  
  # Calculate some performance metrics
  output[["R2"]] = summary(lm(Flux_obs ~ Flux_pred))$r.squared
  output[["MBE"]] = sum(Flux_pred-Flux_obs,na.rm=TRUE)/length(Flux_pred)
  output[["NME"]] = sum(abs(Flux_pred-Flux_obs),na.rm=TRUE)/sum(abs(mean(Flux_obs,na.rm=TRUE)-Flux_obs),na.rm=TRUE)
  output[["SDD"]] = abs(1-sd(Flux_pred,na.rm=TRUE)/sd(Flux_obs,na.rm=TRUE))
  output[["CCO"]] = cor(Flux_pred,Flux_obs,use = "complete.obs", method = "pearson")
  output[["totwithinss"]] = sum(km$WCSS_per_cluster)
  output[["clustersize"]] = clustersize
  output[["clusters"]] = clusters
  output[["AIC"]] = AIC(lm(Flux_obs ~ Flux_pred))
  
  
  # Find the coefficients of each lag in each cluster
  wide_coeff_df = data.frame("Date"=Data %>% 
                               na.omit() %>% 
                               select(Date),
                             "Cluster"=clusters)
  # We also check their significance
  signif_df = data.frame("Date"=Data %>% 
                           na.omit() %>% 
                           select(Date),
                         "Cluster"=clusters)
  
  # Extract coefficients and place in dataframe
  for (i in 1:k){
    wide_coeff_df$R2[wide_coeff_df$Cluster==i] = R2$R2[i]
    if (clustersize[i] > 0){
      # Extract coefficients
      coeff = info[[i]]$model$coefficients 
      # Extract significance
      signif = summary(info[[i]]$model)$coefficients[,4]<0.05 
      
      # Place into dataframes
      wide_coeff_df[wide_coeff_df$Cluster==i,names(coeff)] = coeff %>% 
        unname() %>%
        rep(.,each=nrow(wide_coeff_df[wide_coeff_df$Cluster==i,]))
      signif_df[signif_df$Cluster==i,names(signif)] = signif %>% 
        unname() %>%
        rep(.,each=nrow(signif_df[signif_df$Cluster==i,]))
    }
  }
  
  # Place dataframe into output
  output[["wide_coeff_df"]] = wide_coeff_df
  
  # We perform masking on the coefficient values to ensure that plots only contain
  # relevant and significant information
  # Before applying the masks, we count how many individual values and timesteps
  # are affected by each mask
  masks = list()
  # The coefficients are non-significant in the linear regression
  masks[["not_signif_obs"]] = sum(!signif_df[,4:ncol(signif_df)])
  masks[["not_signif_timesteps"]] = sum(rowSums(!signif_df[,4:ncol(signif_df)])>0)
  # The |coefficients| are less than 0.1
  #maxval = max(abs(wide_coeff_df[,5:ncol(wide_coeff_df)]))
  masks[["low_val_obs"]] = sum(abs(wide_coeff_df[,4:ncol(wide_coeff_df)])<0.1 & abs(wide_coeff_df[,4:ncol(wide_coeff_df)])>0)
  masks[["low_val_timesteps"]] = sum(rowSums(abs(wide_coeff_df[,4:ncol(wide_coeff_df)])<0.1 & abs(wide_coeff_df[,4:ncol(wide_coeff_df)])>0)>0)
  # The R^2 in the cluster linear regression is too low
  masks[["low_R2_timesteps"]] = sum(wide_coeff_df$R2<0.2)
  output[["masks"]] = masks
  # Now that the impact of the masks are recorded, we actually mask
  # Mask non-significant coefficients from linear regression
  wide_coeff_df[,5:ncol(wide_coeff_df)][!signif_df[,4:ncol(signif_df)]] <- 0
  
  # Rearrange
  coeff_df = wide_coeff_df %>% pivot_longer(cols=5:ncol(wide_coeff_df))               
  coeff_df$name = coeff_df$name %>% str_replace(pattern = "climate_cluster",
                                                replacement = "") %>%
    factor(climvars)
  
  # Mask further by absolute value less than 0.1, R^2 lower than 0.2
  # and also remove NaNs
  coeff_df = coeff_df %>%
    mutate(value=replace(value, abs(value)<0.1,0)) %>%
    mutate(value=replace(value,R2 < 0.2,0)) %>%
    mutate(value=replace(value,is.na(value),0))
  
  # Place dataframe into output
  output[["coeff_df"]] = coeff_df
  
  # Plot the coefficient magnitudes and observations vs predictions
  coeff_plot = plotcoeffs(Site,k,coeff_df,Data)
  coeff_plot_Abs = plotcoeffs_Abs(Site,k,coeff_df,Data)
  obsvpred_plot = plotobsvpred(Site,k,Flux_df) 
  
  # Save function inputs
  output[["Site"]]=Site
  output[["k"]]=k
  
  # Save cluster info
#  output[["info"]] = info
  
  # Save outputs
  # Save Rdata
  save(output,file = paste0("outputs/kmeans/",
                            flux,
                            "/",              
                            Site,
                            "_",
                            c("met")[cluster_include_met],
                            c("lag")[cluster_include_lags],
                            "cluster_",
                            c("met")[regress_include_met],
                            c("lag")[regress_include_lags],
                            "regress_",
                            k,
                            "c.Rdata"))
  # Save 3 plots
  ggsave(filename = paste0("outputs/kmeans/",
                           flux,
                           "/images/timeseries/",                            
                           Site,
                           "_",
                           c("met")[cluster_include_met],
                           c("lag")[cluster_include_lags],
                           "cluster_",
                           c("met")[regress_include_met],
                           c("lag")[regress_include_lags],
                           "regress_",
                           k,
                           "c.png"),
         plot = obsvpred_plot,
         device = "png",
         dpi = 320,
         width = 10, height = 5)

  ggsave(filename = paste0("outputs/kmeans/",
                           flux,
                           "/images/coeff/",                            
                           Site,
                           "_",
                           c("met")[cluster_include_met],
                           c("lag")[cluster_include_lags],
                           "cluster_",
                           c("met")[regress_include_met],
                           c("lag")[regress_include_lags],
                           "regress_",
                           k,
                           "c.png"),
         plot = coeff_plot,
         device = "png",
         dpi = 320,
         width = 10, height = 0.5*length(climvars))

  ggsave(filename = paste0("outputs/kmeans/",
                           flux,
                           "/images/coeff/abs/",                            
                           Site,
                           "_",
                           c("met")[cluster_include_met],
                           c("lag")[cluster_include_lags],
                           "cluster_",
                           c("met")[regress_include_met],
                           c("lag")[regress_include_lags],
                           "regress_",
                           k,
                           "c.png"),
         plot = coeff_plot_Abs,
         device = "png",
         dpi = 320,
         width = 10, height = 0.5*length(climvars))
  
}


# Perform regressions in parallel
foo = mcmapply(FUN = function(Site,
                              k,
                              cluster_include_met,
                              cluster_include_lags,
                              regress_include_met,
                              regress_include_lags) worker(Site,
                                                           k,
                                                           cluster_include_met,
                                                           cluster_include_lags,
                                                           regress_include_met,
                                                           regress_include_lags),
               g$sites,
               g$k,
               g$clustmet,
               g$clustlag,
               g$regmet,
               g$reglag,
               mc.cores = 1)

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