### Workflow to plot Memory-Budyko Curve scatterplot

# We need to plot AET/PPT on y-axis and PET/PPT on the x-axis
# For each site we have 3 points for each rainfall memory impact
# Do we jitter these points or do we have 3 separate plots?

### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")
library(grid)
library(cowplot)
library(gridExtra)
library(ggpubr)

### End Setup ###
###############################################################################

# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# SET USER INPUTS FOR WORKFLOW
# 

sitecsv = "longsites.csv"

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load sites
sites = read_csv(paste0("inputs/",
                        sitecsv),
                 col_names=FALSE,
                 show_col_types = F) %>% 
  unlist()

### Plot for NEP 

# Create dataframe to place metric values in
NEE.df = data.frame()

# Load in the NEE metric values
for (Site in sites){
    
      # Load Historical model output
    load(paste0("outputs/kmeans/NEE/",              
                Site,
                "_metlagcluster_metlagregress_729c_nomask.Rdata"))
    
    coeff_df = output$coeff_df 
    
    # Define coefficient groupings
    summary = coeff_df %>% mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                                                  name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "LagT", 
                                                  name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Memory", 
                                                  name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Memory", 
                                                  name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Memory"))
    # Take negative of NEE to produce NEP
    summary$value = -summary$value
    # Calculate quantiles
    summary = summary %>% group_by(Cat) %>%
      summarise(Q0 = quantile(abs(value), probs = 0, na.rm = T), 
                Q5 = quantile(abs(value), probs = 0.05, na.rm = T), 
                Q10 = quantile(abs(value), probs = 0.1, na.rm = T), 
                Q25 = quantile(abs(value), probs = 0.25, na.rm = T), 
                Q50 = quantile(abs(value), probs = 0.5, na.rm = T), 
                Q75 = quantile(abs(value), probs = 0.75, na.rm = T), 
                Q90 = quantile(abs(value), probs = 0.90, na.rm = T), 
                Q95 = quantile(abs(value), probs = 0.95, na.rm = T), 
                Q100 = quantile(abs(value), probs = 1, na.rm = T), 
                Mean = mean(abs(value), na.rm = T), 
                SD = sd(abs(value), na.rm = T)) %>%
      ungroup() %>%
      mutate(Site = Site,
             Cluster = k)
    
    NEE.df = rbind(NEE.df, summary)
    rm(output)
    
}

# Load other data and merge into the dataframe
info = read_csv("inputs/siteinfo.csv",show_col_types = F) %>%
  select(site,MAPrecord,igbp,PPT_VPD_Qle)

NEE.df = merge(NEE.df,info,by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
NEE.df = merge(NEE.df,INFO,by="Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

# Define PFT groupings
NEE.df = NEE.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

# Turn values into factors so they appear in correct order
NEE.df$PFT = factor(NEE.df$PFT, levels = c("Forests", 
                                           "Shrub/Savanna", 
                                           "Grasses", 
                                           "Wetlands"))
NEE.df$Cat = factor(NEE.df$Cat, levels = c("Current",
                                           "LagT", 
                                           "Seasonal~Memory", 
                                           "`Mid-term`~Memory", 
                                           "`Long-term`~Memory"))

# Create a Budyko df
budyko.df = data.frame("AI" = seq(0,4,0.01))
budyko.df$EI = ((1-exp(-budyko.df$AI))*budyko.df$AI*tanh(1/budyko.df$AI))^(1/2)

# Fit regressions to AI
NEE.fit = data.frame("Cat" = NA,
                 "r.squared" = NA,
                 "p.value" = NA,
                 "intercept" = NA,
                 "coefficient" = NA)
NEE.Seasonal = NEE.df[NEE.df$Cat=="Seasonal~Memory",]
NEE.Seasonal.lin.mod = lm(NEE.Seasonal$Q50 ~ NEE.Seasonal$AI)
NEE.fit = rbind(NEE.fit,
                c("Seasonal~Memory",
                  summary(NEE.Seasonal.lin.mod)$r.squared,
                  summary(NEE.Seasonal.lin.mod)$coefficients[2,4],
                  NEE.Seasonal.lin.mod$coefficients[1],
                  NEE.Seasonal.lin.mod$coefficients[2]))


NEE.Midterm = NEE.df[NEE.df$Cat=="`Mid-term`~Memory",]
NEE.Midterm.lin.mod = lm(NEE.Midterm$Q50 ~ NEE.Midterm$AI)
NEE.fit = rbind(NEE.fit,
                c("`Mid-term`~Memory",
                  summary(NEE.Midterm.lin.mod)$r.squared,
                  summary(NEE.Midterm.lin.mod)$coefficients[2,4],
                  NEE.Midterm.lin.mod$coefficients[1],
                  NEE.Midterm.lin.mod$coefficients[2]))


NEE.Longterm = NEE.df[NEE.df$Cat=="`Long-term`~Memory",]
NEE.Longterm.lin.mod = lm(NEE.Longterm$Q50 ~ NEE.Longterm$AI)
NEE.fit = rbind(NEE.fit,
                c("`Long-term`~Memory",
                  summary(NEE.Longterm.lin.mod)$r.squared,
                  summary(NEE.Longterm.lin.mod)$coefficients[2,4],
                  NEE.Longterm.lin.mod$coefficients[1],
                  NEE.Longterm.lin.mod$coefficients[2])) %>%
  na.omit()

# Fit linear models between timescales
NEE.SL.lm = summary(lm(NEE.Longterm$Q50 ~ NEE.Seasonal$Q50))
NEE.SM.lm = summary(lm(NEE.Midterm$Q50 ~ NEE.Seasonal$Q50))
NEE.ML.lm = summary(lm(NEE.Longterm$Q50 ~ NEE.Midterm$Q50))

# Calculate the standard deviation of the timescale medians at each site
NEE.sd.df = NEE.df %>% filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
  group_by(Site,AI) %>%
  summarise(sym_sd = sd(Q50))

NEE.sd.lin.mod = lm(NEE.sd.df$sym_sd ~ NEE.sd.df$AI)
NEE.sd.fit = data.frame("r.squared" = summary(NEE.sd.lin.mod)$r.squared,
                    "p.value" = summary(NEE.sd.lin.mod)$coefficients[2,4],
                    "intercept" = NEE.sd.lin.mod$coefficients[1] ,
                    "coefficient" = NEE.sd.lin.mod$coefficients[2])

# Calculate the mean value for the timescale medians
NEE.Seasonal.Mean = NEE.df %>% filter(Cat == "Seasonal~Memory") %>% select (Q50) %>% colMeans()
NEE.Midterm.Mean = NEE.df %>% filter(Cat == "`Mid-term`~Memory") %>% select (Q50) %>% colMeans()
NEE.Longterm.Mean = NEE.df %>% filter(Cat == "`Long-term`~Memory") %>% select (Q50) %>% colMeans()

# Plot
NEE_main_plot = NEE.df %>% 
      filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
        arrange(Q50) %>%
        ggplot() +
        geom_vline(xintercept=1,linetype="dashed") +
        geom_point(data = NEE.df[NEE.df$Cat=="Seasonal~Memory",],
                   aes(x=AI-0.025,y=EI-0.01,fill=Q50,shape=Cat),
                   alpha=0.9,
                   size=5,
                   stroke=1) +
        geom_point(data = NEE.df[NEE.df$Cat=="`Mid-term`~Memory",],
                   aes(x=AI,y=EI+0.01,fill=Q50,shape=Cat),
                   alpha=0.9,
                   size=5,
                   stroke=1) +
        geom_point(data = NEE.df[NEE.df$Cat=="`Long-term`~Memory",],
                   aes(x=AI+0.025,y=EI-0.01,fill=Q50,shape=Cat),
                   alpha=0.9,
                   size=5,
                   stroke=1) +
        geom_text(data=NEE.fit[NEE.fit$Cat=="Seasonal~Memory",],
                  aes(x=2.75,y=1.25,label=paste0("list(Seasonal~Coeff~Mag == ",
                                                 round(as.numeric(intercept),2),
                                                 " + ",
                                                 round(as.numeric(coefficient),2),
                                                 " %*% AI, r^2 == ",
                                                 round(as.numeric(r.squared),2),
                                                 ",p-value == ",
                                                 round(as.numeric(p.value),2),")")),
                  parse=T,
                  size = 6) +
        geom_text(data=NEE.fit[NEE.fit$Cat=="`Mid-term`~Memory",],
                  aes(x=2.75,y=1.2,label=paste0("list(`Mid-term`~Coeff~Mag == ",
                                                round(as.numeric(intercept),2),
                                                " + ",
                                                round(as.numeric(coefficient),2),
                                                " %*% AI, r^2 == ",
                                                round(as.numeric(r.squared),2),
                                                ",p-value == ",
                                                round(as.numeric(p.value),2),")")),
                  parse=T,
                  size = 6) +
        geom_text(data=NEE.fit[NEE.fit$Cat=="`Long-term`~Memory",],
                  aes(x=2.75,y=1.15,label=paste0("list(`Long-term`~Coeff~Mag == ",
                                                 round(as.numeric(intercept),2),
                                                 " + ",
                                                 round(as.numeric(coefficient),2),
                                                 " %*% AI, r^2 == ",
                                                 round(as.numeric(r.squared),2),
                                                 ",p-value == ",
                                                 round(as.numeric(p.value),2),")")),
                  parse=T,
                  size = 6) +
        geom_text(data=NEE.sd.fit,
                  aes(x=2.75,y=1.1,label=paste0("list(Memory~Std~Dev == ",
                                                round(as.numeric(intercept),2),
                                                " + ",
                                                round(as.numeric(coefficient),2),
                                                " %*% AI, r^2 == ",
                                                round(as.numeric(r.squared),2),
                                                ",p-value == ",
                                                round(as.numeric(p.value),2),")")),
                  parse=T,
                  size = 6) +
        geom_line(data=budyko.df,
                  aes(x=AI,y=EI)) +
        geom_text(aes(x=0.25,y=0.9),
                  label="Energy\nlimited",
                  size=5) +
        geom_text(aes(x=1.5,y=0.1),
                  label="Water\nlimited",
                  size=5) +
        scale_fill_viridis_c(name = "Median Sensitivity Magnitude", 
                             limits = c(0,3),
                             breaks = c(0,0.5,1.5,3),
                             labels = c(0,0.5,1.5,3),
                             option = "inferno",
                             trans = "sqrt",
                             guide = guide_colourbar(direction = "horizontal",
                                                     title.position = "top",
                                                     barwidth = 20)) +
        scale_shape_manual(name = "Memory Length", 
                           labels = c("Long-term","Mid-term","Seasonal"),
                           values = c(21,22,23),
                           guide = guide_legend(title.position = "top",
                                                nrow = 1,
                                                reverse = TRUE)) +
        scale_y_continuous(limits = c(0,1.3),
                           breaks = c(0,0.5,1),
                           expand = c(0,0)) +
        scale_x_continuous(limits = c(0,4),
                           expand = c(0,0)) +
        ylab("Evaporative Index (AET/PPT)") +
        xlab("Aridity Index (PET/PPT)") +
        theme_bw() +
        theme(text = element_text(size = 20),
              legend.position = "top",
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"),
              legend.box.just = "center",
              legend.title.align=0.5,
              legend.box.margin=margin(5,10,5,10))

### Repeat for Latent Heat
# Create dataframe to place metric values in
Qle.df = data.frame()

# Load in the Qle metric values
for (Site in sites){
  for (k in seq){
    
    # Load Historical model output
    load(paste0("outputs/kmeans/Qle/",              
                Site,
                "_metlagcluster_metlagregress_729c_nomask.Rdata"))
    
    coeff_df = output$coeff_df 
    
    # Define coefficient groupings
    summary = coeff_df %>% mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                                                  name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "LagT", 
                                                  name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Memory", 
                                                  name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Memory", 
                                                  name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Memory"))

    # Calculate quantiles
    summary = summary %>% group_by(Cat) %>%
      summarise(Q0 = quantile(abs(value), probs = 0, na.rm = T), 
                Q5 = quantile(abs(value), probs = 0.05, na.rm = T), 
                Q10 = quantile(abs(value), probs = 0.1, na.rm = T), 
                Q25 = quantile(abs(value), probs = 0.25, na.rm = T), 
                Q50 = quantile(abs(value), probs = 0.5, na.rm = T), 
                Q75 = quantile(abs(value), probs = 0.75, na.rm = T), 
                Q90 = quantile(abs(value), probs = 0.90, na.rm = T), 
                Q95 = quantile(abs(value), probs = 0.95, na.rm = T), 
                Q100 = quantile(abs(value), probs = 1, na.rm = T), 
                Mean = mean(abs(value), na.rm = T), 
                SD = sd(abs(value), na.rm = T)) %>%
      ungroup() %>%
      mutate(Site = Site,
             Cluster = k)
    
    Qle.df = rbind(Qle.df, summary)
    rm(output)
    
  }
}

# Load other data and merge into the dataframe
info = read_csv("inputs/siteinfo.csv",show_col_types = F) %>%
  select(site,MAPrecord,igbp,PPT_VPD_Qle)

Qle.df = merge(Qle.df,info,by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
Qle.df = merge(Qle.df,INFO,by="Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

# Define PFT groupings
Qle.df = Qle.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

# Turn values into factors so they appear in correct order
Qle.df$PFT = factor(Qle.df$PFT, levels = c("Forests", 
                                           "Shrub/Savanna", 
                                           "Grasses", 
                                           "Wetlands"))

Qle.df$Cat = factor(Qle.df$Cat, levels = c("Current",
                                           "LagT",
                                           "Seasonal~Memory", 
                                           "`Mid-term`~Memory", 
                                           "`Long-term`~Memory"))




# Create a Budyko df
budyko.df = data.frame("AI" = seq(0,4,0.01))
budyko.df$EI = ((1-exp(-budyko.df$AI))*budyko.df$AI*tanh(1/budyko.df$AI))^(1/2)

# Fit regressions to AI
Qle.fit = data.frame("Cat" = NA,
                     "r.squared" = NA,
                     "p.value" = NA,
                     "intercept" = NA,
                     "coefficient" = NA)
Qle.Seasonal = Qle.df[Qle.df$Cat=="Seasonal~Memory",]
Qle.Seasonal.lin.mod = lm(Qle.Seasonal$Q50 ~ Qle.Seasonal$AI)
Qle.fit = rbind(Qle.fit,
                c("Seasonal~Memory",
                  summary(Qle.Seasonal.lin.mod)$r.squared,
                  summary(Qle.Seasonal.lin.mod)$coefficients[2,4],
                  Qle.Seasonal.lin.mod$coefficients[1],
                  Qle.Seasonal.lin.mod$coefficients[2]))


Qle.Midterm = Qle.df[Qle.df$Cat=="`Mid-term`~Memory",]
Qle.Midterm.lin.mod = lm(Qle.Midterm$Q50 ~ Qle.Midterm$AI)
Qle.fit = rbind(Qle.fit,
                c("`Mid-term`~Memory",
                  summary(Qle.Midterm.lin.mod)$r.squared,
                  summary(Qle.Midterm.lin.mod)$coefficients[2,4],
                  Qle.Midterm.lin.mod$coefficients[1],
                  Qle.Midterm.lin.mod$coefficients[2]))


Qle.Longterm = Qle.df[Qle.df$Cat=="`Long-term`~Memory",]
Qle.Longterm.lin.mod = lm(Qle.Longterm$Q50 ~ Qle.Longterm$AI)
Qle.fit = rbind(Qle.fit,
                c("`Long-term`~Memory",
                  summary(Qle.Longterm.lin.mod)$r.squared,
                  summary(Qle.Longterm.lin.mod)$coefficients[2,4],
                  Qle.Longterm.lin.mod$coefficients[1],
                  Qle.Longterm.lin.mod$coefficients[2])) %>%
  na.omit()

# Fit linear models between timescales
Qle.SL.lm = summary(lm(Qle.Longterm$Q50 ~ Qle.Seasonal$Q50))
Qle.SM.lm = summary(lm(Qle.Midterm$Q50 ~ Qle.Seasonal$Q50))
Qle.ML.lm = summary(lm(Qle.Longterm$Q50 ~ Qle.Midterm$Q50))

# Calculate the standard deviation of the timescale medians at each site
Qle.sd.df = Qle.df %>% filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
  group_by(Site,AI) %>%
  summarise(sym_sd = sd(Q50))

Qle.sd.lin.mod = lm(Qle.sd.df$sym_sd ~ Qle.sd.df$AI)
Qle.sd.fit = data.frame("r.squared" = summary(Qle.sd.lin.mod)$r.squared,
                        "p.value" = summary(Qle.sd.lin.mod)$coefficients[2,4],
                        "intercept" = Qle.sd.lin.mod$coefficients[1] ,
                        "coefficient" = Qle.sd.lin.mod$coefficients[2])

# Calculate the mean value for the timescale medians
Qle.Seasonal.Mean = Qle.df %>% filter(Cat == "Seasonal~Memory") %>% select (Q50) %>% colMeans()
Qle.Midterm.Mean = Qle.df %>% filter(Cat == "`Mid-term`~Memory") %>% select (Q50) %>% colMeans()
Qle.Longterm.Mean = Qle.df %>% filter(Cat == "`Long-term`~Memory") %>% select (Q50) %>% colMeans()

# Plot
Qle_main_plot = Qle.df %>% filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
  arrange(Q50) %>%
  ggplot() +
  geom_vline(xintercept=1,linetype="dashed") +
  geom_point(data = Qle.df[Qle.df$Cat=="Seasonal~Memory",],
             aes(x=AI-0.025,y=EI-0.01,fill=Q50,shape=Cat),alpha=0.9,size=5,stroke=1) +
  geom_point(data = Qle.df[Qle.df$Cat=="`Mid-term`~Memory",],
             aes(x=AI,y=EI+0.01,fill=Q50,shape=Cat),alpha=0.9,size=5,stroke=1) +
  geom_point(data = Qle.df[Qle.df$Cat=="`Long-term`~Memory",],
             aes(x=AI+0.025,y=EI-0.01,fill=Q50,shape=Cat),alpha=0.9,size=5,stroke=1) +
  geom_text(data=Qle.fit[Qle.fit$Cat=="Seasonal~Memory",],
            aes(x=2.75,y=1.25,label=paste0("list(Seasonal~Coeff~Mag == ",
                                         round(as.numeric(intercept),2),
                                         " + ",
                                         round(as.numeric(coefficient),2),
                                         " %*% AI, r^2 == ",
                                         round(as.numeric(r.squared),2),
                                         ",p-value == ",
                                         round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_text(data=Qle.fit[Qle.fit$Cat=="`Mid-term`~Memory",],
            aes(x=2.75,y=1.2,label=paste0("list(`Mid-term`~Coeff~Mag == ",
                                          round(as.numeric(intercept),2),
                                          " + ",
                                          round(as.numeric(coefficient),2),
                                          " %*% AI, r^2 == ",
                                          round(as.numeric(r.squared),2),
                                          ",p-value == ",
                                          round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_text(data=Qle.fit[Qle.fit$Cat=="`Long-term`~Memory",],
            aes(x=2.75,y=1.15,label=paste0("list(`Long-term`~Coeff~Mag == ",
                                            round(as.numeric(intercept),2),
                                            " + ",
                                            round(as.numeric(coefficient),2),
                                            " %*% AI, r^2 == ",
                                            round(as.numeric(r.squared),2),
                                            ",p-value == ",
                                            round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_text(data=Qle.sd.fit,
            aes(x=2.75,y=1.1,label=paste0("list(Memory~Std~Dev == ",
                                        round(as.numeric(intercept),2),
                                        " + ",
                                        round(as.numeric(coefficient),2),
                                        " %*% AI, r^2 == ",
                                        round(as.numeric(r.squared),2),
                                        ",p-value == ",
                                        round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_line(data=budyko.df,
            aes(x=AI,y=EI)) +
  geom_text(aes(x=0.25,y=0.9),
            label="Energy\nlimited",
            size=5) +
  geom_text(aes(x=1.5,y=0.1),
            label="Water\nlimited",
            size=5) +
  scale_fill_viridis_c(name = "Median Sensitivity Magnitude", 
                       limits = c(0,3),
                       breaks = c(0,0.5,1.5,3),
                       labels = c(0,0.5,1.5,3),
                       option = "inferno",
                       trans = "sqrt",
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = "top",
                                               barwidth = 15)) +
  scale_shape_manual(name = "Memory Length",
                     labels = c("Long-term","Mid-term","Seasonal"),
                     values=c(21,22,23),
                     guide = guide_legend(title.position = "top",
                                          nrow = 1,
                                          reverse = TRUE)) +
  scale_y_continuous(limits = c(0,1.3),
                     breaks=c(0,0.5,1),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,4),
                     expand=c(0,0)) +
  ylab("Evaporative Index (AET/PPT)") +
  xlab("Aridity Index (PET/PPT)") +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.box.just = "center",
        legend.title.align=0.5,
        legend.box.margin=margin(5,10,5,10))

### Combine into a single figure

# Load the inset data
load("outputs/Figure3_Qle_inset.Rdata")
load("outputs/Figure3_NEE_inset.Rdata")

# Combine the main and inset plots
NEE_plot = NEE_main_plot + 
      theme(legend.position = "none") + 
      annotation_custom(ggplotGrob(NEE_box),
                        xmin=2.12,
                        xmax = 3.95, 
                        ymin = 0.04,
                        ymax = 0.715)
Qle_plot = Qle_main_plot + 
      theme(legend.position = "none") + 
      annotation_custom(ggplotGrob(Qle_box),
                        xmin=2.12,
                        xmax = 3.95,
                        ymin = 0.04, 
                        ymax = 0.715)
legend = get_legend(NEE_main_plot)

# Save the plot
png(filename = "images/Figure3.png",
    width = 16,
    height = 16,
    units = "in",
    res = 320)
print(ggarrange(legend,
                NEE_plot,
                Qle_plot,
                nrow = 3,
                heights = c(1,5,5), 
                labels = c(NA,"(a)","(b)"),
                font.label = list(size = 20)))
dev.off()


# Check correlations across fluxes
Flux.S = summary(lm(NEE.df$Q50[NEE.df$Cat=="Seasonal~Memory"] ~ Qle.df$Q50[Qle.df$Cat=="Seasonal~Memory"]))
Flux.M = summary(lm(NEE.df$Q50[NEE.df$Cat=="`Mid-term`~Memory"] ~ Qle.df$Q50[Qle.df$Cat=="`Mid-term`~Memory"]))
Flux.L = summary(lm(NEE.df$Q50[NEE.df$Cat=="`Long-term`~Memory"] ~ Qle.df$Q50[Qle.df$Cat=="`Long-term`~Memory"]))

# Check correlations across timescales within fluxes
Corr.df = data.frame()
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Flux_S",
                           "coeff" = Flux.S$coefficients[2,1],
                           "r2" = Flux.S$r.squared,
                           "p" = Flux.S$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Flux_M",
                           "coeff" = Flux.M$coefficients[2,1],
                           "r2" = Flux.M$r.squared, 
                           "p" = Flux.M$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Flux_L",
                           "coeff" = Flux.L$coefficients[2,1],
                           "r2" = Flux.L$r.squared, 
                           "p" = Flux.L$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "NEE_SM",
                           "coeff" = NEE.SM.lm$coefficients[2,1],
                           "r2" = NEE.SM.lm$r.squared, 
                           "p" = NEE.SM.lm$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "NEE_SL",
                           "coeff" = NEE.SL.lm$coefficients[2,1],
                           "r2" = NEE.SL.lm$r.squared, 
                           "p" = NEE.SL.lm$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "NEE_ML",
                           "coeff" = NEE.ML.lm$coefficients[2,1],
                           "r2" = NEE.ML.lm$r.squared, 
                           "p" = NEE.ML.lm$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Qle_SM",
                           "coeff" = Qle.SM.lm$coefficients[2,1],
                           "r2" = Qle.SM.lm$r.squared, 
                           "p" = Qle.SM.lm$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Qle_SL",
                           "coeff" = Qle.SL.lm$coefficients[2,1],
                           "r2" = Qle.SL.lm$r.squared, 
                           "p" = Qle.SL.lm$coefficients[2,4]))
Corr.df = rbind(Corr.df,
                data.frame("Corr" = "Qle_ML",
                           "coeff" = Qle.ML.lm$coefficients[2,1],
                           "r2" = Qle.ML.lm$r.squared,
                           "p" = Qle.ML.lm$coefficients[2,4]))

# Check mean Q50 sensitivities at different timescales
NEE.df %>% filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
  summarise(NEPAllMean = mean(Q50),NEPAllMin = min(Q50),NEPAllMax = max(Q50))

NEE.df %>% filter(Cat %in% c("Seasonal~Memory")) %>%
  summarise(NEPSeasonalMean = mean(Q50),NEPSeasonalMin = min(Q50),NEPSeasonalMax = max(Q50))
NEE.df %>% filter(Cat %in% c("`Mid-term`~Memory")) %>%
  summarise(NEPMidtermMean = mean(Q50),NEPMidtermMin = min(Q50),NEPMidtermMax = max(Q50))
NEE.df %>% filter(Cat %in% c("`Long-term`~Memory")) %>%
  summarise(NEPLongtermMean = mean(Q50),NEPLongtermMin = min(Q50),NEPLongtermMax = max(Q50))

Qle.df %>% filter(Cat %in% c("Seasonal~Memory","`Mid-term`~Memory","`Long-term`~Memory")) %>%
  summarise(LEAllMean = mean(Q50),LEAllMin = min(Q50),LEAllMax = max(Q50))

Qle.df %>% filter(Cat %in% c("Seasonal~Memory")) %>%
  summarise(LESeasonalMean = mean(Q50),LESeasonalMin = min(Q50),LESeasonalMax = max(Q50))
Qle.df %>% filter(Cat %in% c("`Mid-term`~Memory")) %>%
  summarise(LEMidtermMean = mean(Q50),LEMidtermMin = min(Q50),LEMidtermMax = max(Q50))
Qle.df %>% filter(Cat %in% c("`Long-term`~Memory")) %>%
  summarise(LELongtermMean = mean(Q50),LELongtermMin = min(Q50),LELongtermMax = max(Q50))







