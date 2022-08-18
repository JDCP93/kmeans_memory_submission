### Workflow to plot Legacy-PFT boxplots

# We need to plot legacy categories on the x-axis, group by PFT and have 
# coefficient values on the y-axis

### Begin Setup ###
# Tidy up
rm(list = ls())

# Message
Start = Sys.time()
message("Starting workflow at ", Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")

### End Setup ###
###############################################################################

# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# SET USER INPUTS FOR WORKFLOW
# 

sitecsv = "longsites.csv"
seq = c(729)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load sites
sites = read_csv(paste0("inputs/", 
                        sitecsv), 
                 col_names = FALSE, 
                 show_col_types = F) %>% 
  unlist()

# Create dataframe to place metric values in
NEE.df = data.frame()

# Load in the metric values
for (Site in sites){
  for (k in seq){
    
    # Load model outputs
    load(paste0("outputs/kmeans/NEE/",               
                Site, 
                "_metlagcluster_metlagregress_", 
                k, 
                "c.Rdata"))
    
    summary = output$coeff_df %>% 
      mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                             name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "LagT", 
                             name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Legacy", 
                             name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Legacy", 
                             name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Legacy"), 
             Site = Site, 
             Cluster = k)
    
    # Invert value for NEE so it is now NEP
    summary$value = -summary$value
    
    # Bind into single dataframe
    NEE.df = rbind(NEE.df, summary)
    rm(output,summary)
    
  }
}

# Load and merge site info
info = read_csv("inputs/siteinfo.csv", show_col_types = F) %>%
  select(site, MAPrecord, igbp, PPT_VPD_Qle)

NEE.df = merge(NEE.df, info, by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
NEE.df = merge(NEE.df, INFO, by = "Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

# Group PFTs
NEE.df = NEE.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

# Create factors so that plots are in correct order
NEE.df$IGBP = factor(NEE.df$IGBP, levels = c("dbf", "ebf", "enf", "mf", "sav", "shr", "gra", "cro", "wet"))
NEE.df$PFT = factor(NEE.df$PFT, levels = c("Forests", "Shrub/Savanna", "Grasses", "Wetlands"))
NEE.df$Cat = factor(NEE.df$Cat, levels = c("Current", "LagT", "Seasonal~Legacy", "`Mid-term`~Legacy", "`Long-term`~Legacy"))

# Plot boxplots
NEE_box = NEE.df %>% filter(!(value==0)) %>%
  group_by(Cat) %>%
  summarise(Q0 = quantile(value, probs = 0), 
            Q1 = quantile(value, probs = 0.01), 
            Q5 = quantile(value, probs = 0.05), 
            Q10 = quantile(value, probs = 0.1), 
            Q25 = quantile(value, probs = 0.25), 
            Q50 = quantile(value, probs = 0.5), 
            Q75 = quantile(value, probs = 0.75), 
            Q90 = quantile(value, probs = 0.90), 
            Q95 = quantile(value, probs = 0.95), 
            Q99 = quantile(value, probs = 0.99), 
            Q100 = quantile(value, probs = 1), 
            Mean = mean(value), 
            SD = sd(value)) %>% 
  filter(Cat %in% c("Seasonal~Legacy", "`Mid-term`~Legacy", "`Long-term`~Legacy")) %>%
  ggplot() +
  geom_errorbar(aes(x = Cat, ymin = Q90, ymax = Q99, group = Cat, color = Cat), position = "dodge", width = 0.25) +
  geom_errorbar(aes(x = Cat, ymin = Q1, ymax = Q10, group = Cat, color = Cat), position = "dodge", width = 0.25) +
  geom_crossbar(aes(x = Cat, ymin = Q10, y = Q50, ymax = Q90, group = Cat, fill = Cat), position = "dodge", alpha = 0.6, width = 0.25) +
  geom_crossbar(aes(x = Cat, ymin = Q25, y = Q50, ymax = Q75, group = Cat, fill = Cat), position = "dodge", width = 0.5) +
  geom_point(aes(x = Cat, y = Mean, group = Cat, fill = Cat), shape = 4, stroke = 2, show.legend = FALSE) +
  scale_fill_viridis_d(name = "Category",begin = 0.2, end = 0.8, direction = -1, option = "viridis") +
  scale_color_viridis_d(name = "Category",begin = 0.2, end = 0.8, direction = -1, option = "viridis") +
  scale_x_discrete(labels = c("Seasonal","Mid-term","Long-term")) +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0), trans = "pseudo_log") +
  ylab("Sensitivity (-)") +
  xlab("Rainfall Category") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        text = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(c(10, 10, 10, 10), "pt"),
        plot.background = element_rect(colour = "black", size=1))


# REPEAT FOR LATENT HEAT

# Create dataframe to place metric values in
Qle.df = data.frame()

# Load in the metric values
for (Site in sites){
  for (k in seq){
    
    # Load model outputs
    load(paste0("outputs/kmeans/Qle/",               
                Site, 
                "_metlagcluster_metlagregress_", 
                k, 
                "c.Rdata"))
    
    summary = output$coeff_df %>% 
      mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                             name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "LagT", 
                             name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Legacy", 
                             name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Legacy", 
                             name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Legacy"), 
             Site = Site, 
             Cluster = k)
    
    Qle.df = rbind(Qle.df, summary)
    rm(output,summary)
    
  }
}



Qle.df = merge(Qle.df, info, by.x = "Site", by.y = "site")

Qle.df = merge(Qle.df, INFO, by = "Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))


Qle.df = Qle.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

Qle.df$IGBP = factor(Qle.df$IGBP, levels = c("dbf", "ebf", "enf", "mf", "sav", "shr", "gra", "cro", "wet"))
Qle.df$PFT = factor(Qle.df$PFT, levels = c("Forests", "Shrub/Savanna", "Grasses", "Wetlands"))
Qle.df$Cat = factor(Qle.df$Cat, levels = c("Current", "LagT", "Seasonal~Legacy", "`Mid-term`~Legacy", "`Long-term`~Legacy"))

Qle_box = Qle.df %>% filter(!(value==0)) %>%
  group_by(Cat) %>%
  summarise(Q0 = quantile(value, probs = 0),
            Q1 = quantile(value, probs = 0.01), 
            Q5 = quantile(value, probs = 0.05), 
            Q10 = quantile(value, probs = 0.1), 
            Q25 = quantile(value, probs = 0.25), 
            Q50 = quantile(value, probs = 0.5), 
            Q75 = quantile(value, probs = 0.75), 
            Q90 = quantile(value, probs = 0.90), 
            Q95 = quantile(value, probs = 0.95), 
            Q99 = quantile(value, probs = 0.99), 
            Q100 = quantile(value, probs = 1), 
            Mean = mean(value), 
            SD = sd(value)) %>% 
  filter(Cat %in% c("Seasonal~Legacy", "`Mid-term`~Legacy", "`Long-term`~Legacy")) %>%
  ggplot() +
  geom_errorbar(aes(x = Cat, ymin = Q90, ymax = Q99, group = Cat, color = Cat), position = "dodge", width = 0.25) +
  geom_errorbar(aes(x = Cat, ymin = Q1, ymax = Q10, group = Cat, color = Cat), position = "dodge", width = 0.25) +
  geom_crossbar(aes(x = Cat, ymin = Q10, y = Q50, ymax = Q90, group = Cat, fill = Cat), position = "dodge", alpha = 0.6, width = 0.25) +
  geom_crossbar(aes(x = Cat, ymin = Q25, y = Q50, ymax = Q75, group = Cat, fill = Cat), position = "dodge", width = 0.5) +
  geom_point(aes(x = Cat, y = Mean, group = Cat, fill = Cat), shape = 4, stroke = 2, show.legend = FALSE) +
  scale_fill_viridis_d(name = "Category",begin = 0.2, end = 0.8, direction = -1, option = "viridis") +
  scale_color_viridis_d(name = "Category",begin = 0.2, end = 0.8, direction = -1, option = "viridis") +
  scale_x_discrete(labels = c("Seasonal","Mid-term","Long-term")) +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0), trans = "pseudo_log") +
  ylab("Sensitivity (-)") +
  xlab("Rainfall Category") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        text = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(c(10, 10, 10, 10), "pt"),
        plot.background = element_rect(colour = "black", size=1))

# Save plots and data

ggsave("Qle_boxplot_inset.png",plot=Qle_box)
ggsave("NEE_boxplot_inset.png",plot=NEE_box)

save(Qle_box,file = "Qle_boxplot_inset.Rdata")
save(NEE_box,file = "NEE_boxplot_inset.Rdata")
