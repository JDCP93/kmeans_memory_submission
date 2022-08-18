################################################################################

### Modelling functions for use with k-means clustering

################################################################################

smallclusters = function(clusters,clustersize,Climate,Current,df,i,km){

  # clusters = kmeans() output
  # Climate = predictor matrix
  # df = dataframe of observations of predicted flux
  # i = cluster number in question

  # Find the climate in the cluster
  climate_cluster_actual = Climate[clusters==i,]
  # Find the Flux in the cluster
  Flux_cluster_actual = df[clusters==i,2]

  # Extract centre
  centre = km$centroids[i,]
  # Calculate the distances from each observation to the centre
  distances = rowSums((sweep(Current,2,centre))^2)
  # Find the closest
  closest = which(distances %in% sort(distances)[1:(10*ncol(Climate))])
  closestdistances = distances[closest]
  # Find the closest which aren't already in the cluster
  nocigar = closest[!(closest %in% clusters[clusters==i])]
  nocigardistances = distances[nocigar]
  # Add the closest ones to the cluster which aren't already in it
  toaddindex = which(nocigardistances %in% sort(nocigardistances)[1:(10*ncol(Climate)-clustersize[i])])
  toadd = nocigar[toaddindex]
  
  # Create the expanded clusters
  climate_cluster = rbind(climate_cluster_actual,Climate[toadd,])
  Flux_cluster = c(Flux_cluster_actual,df[toadd,2])

  # Fit the linear regression model in the cluster
  lin.mod = lm(Flux_cluster ~ climate_cluster, na.action = na.exclude)
  # Calculate the r.squared in the cluster
  r.squared = summary(lin.mod)$r.squared
  # Extract the predicted values for the actual cluster
  Flux_pred = fitted(lin.mod)[1:clustersize[i]]
  
  # Calculate information criterion
  AIC = AIC(lin.mod)
  BIC = BIC(lin.mod)
  # Calculate performance metrics
  R2 = summary(lm(Flux_cluster_actual ~ Flux_pred))$r.squared
  MBE = sum(Flux_pred-Flux_cluster_actual,na.rm=TRUE)/length(Flux_pred)
  NME = sum(abs(Flux_pred-Flux_cluster_actual),na.rm=TRUE)/sum(abs(mean(Flux_cluster_actual,na.rm=TRUE)-Flux_cluster_actual),na.rm=TRUE)
  SDD = abs(1-sd(Flux_pred,na.rm=TRUE)/sd(Flux_cluster_actual,na.rm=TRUE))
  CCO = cor(Flux_pred,Flux_cluster_actual,use = "complete.obs", method = "pearson")

  # Assign and output the cluster info
  output = list("climate" = climate_cluster_actual,
                "Flux_obs" = Flux_cluster_actual,
                "Flux_pred" = Flux_pred,
                "model" = lin.mod,
                "expanded.r.squared" = r.squared,
                "mod.r" = r.squared*nrow(climate_cluster_actual),
                "AIC" = AIC,
                "BIC" = BIC,
                "R2" = R2,
                "CCO" = CCO,
                "MBE" = MBE,
                "NME" = NME,
                "SDD" = SDD)

}

normalclusters = function(clusters,Climate,df,i){

  # clusters = kmeans() output
  # Climate = predictor matrix
  # df = dataframe of observations of predicted flux
  # i = cluster number in question

  # Find the climate in the cluster
  climate_cluster = Climate[clusters==i,]
  # Find the Flux in the cluster
  Flux_cluster = df$Flux[clusters==i]

  # Fit the linear regression model in the cluster
  lin.mod = lm(Flux_cluster ~ climate_cluster, na.action = na.exclude)
  # Calculate the r.squared in the cluster
  r.squared = summary(lin.mod)$r.squared
  # Place the k-means fitted Flux into the data frame
  Flux_pred = fitted(lin.mod)
  # Calculate information criterion
  AIC = AIC(lin.mod)
  BIC = BIC(lin.mod)

  R2 = summary(lm(Flux_cluster ~ Flux_pred))$r.squared
  MBE = sum(Flux_pred-Flux_cluster,na.rm=TRUE)/length(Flux_pred)
  NME = sum(abs(Flux_pred-Flux_cluster),na.rm=TRUE)/sum(abs(mean(Flux_cluster,na.rm=TRUE)-Flux_cluster),na.rm=TRUE)
  SDD = abs(1-sd(Flux_pred,na.rm=TRUE)/sd(Flux_cluster,na.rm=TRUE))
  CCO = cor(Flux_pred,Flux_cluster,use = "complete.obs", method = "pearson")

  # Assign and output the cluster info
  output = list("climate" = climate_cluster,
                "Flux_obs" = Flux_cluster,
                "Flux_pred" = Flux_pred,
                "model" = lin.mod,
                "r.squared" = r.squared,
                "mod.r" = r.squared*nrow(climate_cluster),
                "AIC" = AIC,
                "BIC" = BIC,
                "R2" = R2,
                "CCO" = CCO,
                "MBE" = MBE,
                "NME" = NME,
                "SDD" = SDD)

}

