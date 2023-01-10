### Workflow to plot Memory-PFT boxplots

# We need to plot Memory categories on the x-axis, group by PFT and have 
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

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load sites
sites = read_csv(paste0("inputs/", 
                        sitecsv), 
                 col_names = FALSE, 
                 show_col_types = F) %>% 
  unlist()

### Plot for NEP

# Create dataframe to place metric values in
NEE.df = data.frame()

# Load in the metric values
for (Site in sites){
    
      # Load Historical model output
    load(paste0("outputs/kmeans/NEE/",               
                Site, 
                "_metlagcluster_metlagregress_729c.Rdata"))
    
    coeff_df = output$coeff_df 
    
    # Define coefficient groupings
    summary = coeff_df %>% mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                                                  name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "Temperature~Memory", 
                                                  name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Precipitation~Memory", 
                                                  name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Precipitation~Memory", 
                                                  name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Precipitation~Memory"), 
                                  Site = Site, 
                                  Cluster = k) 
    
    # Take negative of NEE to produce NEP
    summary$value = -summary$value
    
    NEE.df = rbind(NEE.df, summary)
    rm(output)
}

# Load other data and merge into the dataframe
info = read_csv("inputs/siteinfo.csv", show_col_types = F) %>%
  select(site, MAPrecord, igbp, PPT_VPD_Qle)

NEE.df = merge(NEE.df, info, by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
NEE.df = merge(NEE.df, INFO, by = "Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

# Define PFT groupings
NEE.df = NEE.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

# Turn values into factors so they appear in correct order
NEE.df$IGBP = factor(NEE.df$IGBP, levels = c("dbf", 
                                             "ebf", 
                                             "enf",
                                             "mf", 
                                             "sav", 
                                             "shr",
                                             "gra", 
                                             "cro", 
                                             "wet"))

NEE.df$PFT = factor(NEE.df$PFT, levels = c("Forests", 
                                           "Shrub/Savanna",
                                           "Grasses", 
                                           "Wetlands"))

NEE.df$Cat = factor(NEE.df$Cat, levels = c("Current", 
                                           "Temperature~Memory", 
                                           "Seasonal~Precipitation~Memory", 
                                           "`Mid-term`~Precipitation~Memory", 
                                           "`Long-term`~Precipitation~Memory"))

# Find the number of sites and site-years for each PFT
siteyears = merge(info,INFO,by.x = "site", by.y= "Site" ) %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr")) %>%
  group_by(igbp) %>% 
  summarise(nsites = length(Yrs),
            nyears=sum(Yrs-4))

NEE.df = merge(NEE.df,siteyears,by.x="IGBP",by.y="igbp")

# Create plot
NEE_box = NEE.df %>% group_by(IGBP, Cat, PFT, nyears, nsites) %>%
  summarise(Q0 = quantile(value, probs = 0, na.rm = T), 
            Q1 = quantile(value, probs = 0.01, na.rm = T), 
            Q5 = quantile(value, probs = 0.05, na.rm = T), 
            Q10 = quantile(value, probs = 0.1, na.rm = T), 
            Q25 = quantile(value, probs = 0.25, na.rm = T), 
            Q50 = quantile(value, probs = 0.5, na.rm = T), 
            Q75 = quantile(value, probs = 0.75, na.rm = T), 
            Q90 = quantile(value, probs = 0.90, na.rm = T), 
            Q95 = quantile(value, probs = 0.95, na.rm = T), 
            Q99 = quantile(value, probs = 0.99, na.rm = T), 
            Q100 = quantile(value, probs = 1, na.rm = T), 
            Mean = mean(value, na.rm = T), 
            SD = sd(value, na.rm = T)) %>% 
  filter(Cat %in% c("Temperature~Memory",
                    "Seasonal~Precipitation~Memory", 
                    "`Mid-term`~Precipitation~Memory", 
                    "`Long-term`~Precipitation~Memory")) %>%
  ggplot() +
  geom_errorbar(aes(x = IGBP, ymin = Q90, ymax = Q95, group = IGBP, color = PFT),
                position = "dodge",
                width = 0.5) +
  geom_errorbar(aes(x = IGBP, ymin = Q5, ymax = Q10, group = IGBP, color = PFT),
                position = "dodge",
                width = 0.5) +
  geom_crossbar(aes(x = IGBP, ymin = Q10, y = Q50, ymax = Q90, group = IGBP, fill = PFT),
                position = "dodge", 
                alpha = 0.6, 
                width = 0.5) +
  geom_crossbar(aes(x = IGBP, ymin = Q25, y = Q50, ymax = Q75, group = IGBP, fill = PFT), 
                position = "dodge") +
  geom_point(aes(x = IGBP, y = Mean, group = IGBP, fill = PFT), 
             shape = 4,
             stroke = 2, 
             show.legend = FALSE) +
  facet_grid(.~Cat,
             labeller = label_parsed) +
  scale_fill_viridis_d(name = "PFT Group",
                       begin = 0.2,
                       end = 0.95, 
                       direction = -1, 
                       option = "inferno") +
  scale_color_viridis_d(name = "PFT Group",
                        begin = 0.2,
                        end = 0.95, 
                        direction = -1,
                        option = "inferno") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        text = element_text(size = 20),
        axis.text.x = element_text(angle = -45, vjust = 0),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.box.just = "center",
        legend.title.align=0.5,
        legend.box.margin=margin(5,10,5,10)) +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0)) +
  ylab("Sensitivity (-)")

### Repeat for Latent Heat

# Create dataframe to place metric values in
Qle.df = data.frame()

# Load in the metric values
for (Site in sites){
    
      # Load Historical model output
    load(paste0("outputs/kmeans/Qle/",               
                Site, 
                "_metlagcluster_metlagregress_729c.Rdata"))
    
    coeff_df = output$coeff_df 
    # Define coefficient groupings
    summary = coeff_df %>% mutate(Cat = case_when(name %in% c("LAI_MOD", "SWdown", "VPD", "Wind", "Tair","Precip") ~ "Current", 
                                                  name %in% c("T1", "T2to7", "T8to14", "T15to30") ~ "Temperature~Memory", 
                                                  name %in% c("P1to30", "P31to90", "P91to180") ~ "Seasonal~Precipitation~Memory", 
                                                  name %in% c("P181to365", "P366to730") ~ "`Mid-term`~Precipitation~Memory", 
                                                  name %in% c("P731to1095", "P1096to1460") ~ "`Long-term`~Precipitation~Memory"), 
                                  Site = Site, 
                                  Cluster = k)
    
    Qle.df = rbind(Qle.df, summary)
    rm(output)
}

# Load other data and merge into the dataframe
info = read_csv("inputs/siteinfo.csv", show_col_types = F) %>%
  select(site, MAPrecord, igbp, PPT_VPD_Qle)

Qle.df = merge(Qle.df, info, by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
Qle.df = merge(Qle.df, INFO, by = "Site") %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

# Define PFT groupings
Qle.df = Qle.df %>% mutate(PFT = case_when(igbp %in% c("dbf", "ebf", "enf", "mf") ~ "Forests", 
                                           igbp %in% c("cro", "gra") ~ "Grasses", 
                                           igbp %in% c("sav","shr") ~ "Shrub/Savanna", 
                                           igbp %in% c("wet") ~ "Wetlands")) %>%
  rename(IGBP = igbp)

# Turn values into factors so they appear in correct order
Qle.df$IGBP = factor(Qle.df$IGBP, levels = c("dbf", 
                                             "ebf", 
                                             "enf", 
                                             "mf", 
                                             "sav", 
                                             "shr",
                                             "gra", 
                                             "cro", 
                                             "wet"))

Qle.df$PFT = factor(Qle.df$PFT, levels = c("Forests", 
                                           "Shrub/Savanna", 
                                           "Grasses", 
                                           "Wetlands"))

Qle.df$Cat = factor(Qle.df$Cat, levels = c("Current", 
                                           "Temperature~Memory",
                                           "Seasonal~Precipitation~Memory", 
                                           "`Mid-term`~Precipitation~Memory",
                                           "`Long-term`~Precipitation~Memory"))

# Find the number of sites and site-years for each PFT
siteyears = merge(info,INFO,by.x = "site", by.y= "Site" ) %>% 
  mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
  mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr")) %>%
  group_by(igbp) %>% 
  summarise(nsites = length(Yrs),
            nyears=sum(Yrs-4))

Qle.df = merge(Qle.df,siteyears,by.x="IGBP",by.y="igbp")

# Create plot
Qle_box = Qle.df %>% group_by(IGBP, Cat, PFT, nyears, nsites) %>%
  summarise(Q0 = quantile(value, probs = 0, na.rm = T), 
            Q1 = quantile(value, probs = 0.01, na.rm = T), 
            Q5 = quantile(value, probs = 0.05, na.rm = T), 
            Q10 = quantile(value, probs = 0.1, na.rm = T), 
            Q25 = quantile(value, probs = 0.25, na.rm = T), 
            Q50 = quantile(value, probs = 0.5, na.rm = T), 
            Q75 = quantile(value, probs = 0.75, na.rm = T), 
            Q90 = quantile(value, probs = 0.90, na.rm = T), 
            Q95 = quantile(value, probs = 0.95, na.rm = T), 
            Q99 = quantile(value, probs = 0.99, na.rm = T), 
            Q100 = quantile(value, probs = 1, na.rm = T), 
            Mean = mean(value, na.rm = T), 
            SD = sd(value, na.rm = T)) %>% 
  filter(Cat %in% c("Temperature~Memory",
                    "Seasonal~Precipitation~Memory", 
                    "`Mid-term`~Precipitation~Memory",
                    "`Long-term`~Precipitation~Memory")) %>%
  ggplot() +
  geom_errorbar(aes(x = IGBP, ymin = Q90, ymax = Q95, group = IGBP, color = PFT), 
                position = "dodge",
                width = 0.5) +
  geom_errorbar(aes(x = IGBP, ymin = Q5, ymax = Q10, group = IGBP, color = PFT), 
                position = "dodge",
                width = 0.5) +
  geom_crossbar(aes(x = IGBP, ymin = Q10, y = Q50, ymax = Q90, group = IGBP, fill = PFT),
                position = "dodge",
                alpha = 0.6,
                width = 0.5) +
  geom_crossbar(aes(x = IGBP, ymin = Q25, y = Q50, ymax = Q75, group = IGBP, fill = PFT), 
                position = "dodge") +
  geom_point(aes(x = IGBP, y = Mean, group = IGBP, fill = PFT), 
             shape = 4, 
             stroke = 2, 
             show.legend = FALSE) +
  geom_label(aes(x = IGBP, y = -4,label=nyears),
             size = 5) +
  geom_label(aes(x = IGBP, y = -3,label=nsites),
             size = 5) +
  facet_grid(.~Cat,
             labeller = label_parsed) +
  scale_fill_viridis_d(name = "PFT Group",
                       begin = 0.2, 
                       end = 0.95, 
                       direction = -1, 
                       option = "inferno") +
  scale_color_viridis_d(name = "PFT Group",
                        begin = 0.2, 
                        end = 0.95,
                        direction = -1,
                        option = "inferno") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        text = element_text(size = 20),
        axis.text.x = element_text(angle = -45, vjust = 0),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.box.just = "center",
        legend.title.align=0.5,
        legend.box.margin=margin(5,10,5,10)) +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0)) +
  ylab("Sensitivity (-)")


### Combine into a single figure
library(cowplot)
library(ggpubr)

legend = get_legend(NEE_box)
NEE_box = NEE_box + theme(legend.position = "none")
Qle_box = Qle_box + theme(legend.position = "none")

png(filename = "images/Figure4.png",
    width = 16,
    height = 16,
    units = "in",
    res = 320)

print(
      ggarrange(legend,
                NEE_box,
                Qle_box,
                nrow = 3,
                heights = c(1,7,7), 
                labels = c(NA,"(a)","(b)"), 
                font.label = list(size = 20))
      )

dev.off()

