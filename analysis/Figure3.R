### Workflow to plot Legacy-Budyko Curve scatterplot

# We need to plot AET/PPT on y-axis and PET/PPT on the x-axis
# For each site we have 3 points for each rainfall legacy impact

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
 
# Create dataframe to place metric values in
NEE.df = data.frame()

# Load in the NEE metric values
for (Site in sites){
  for (k in seq){
    
    # Metlag cluster and metlag regression
    load(paste0("outputs/kmeans/NEE/",              
                Site,
                "_metlagcluster_metlagregress_",
                k,
                "c.Rdata"))
    
    # Categorise the predictor variables and summarise them
    summary = output$coeff_df %>% 
      mutate(Cat = case_when(name %in% c("LAI_MOD","SWdown","VPD","Wind","Tair","Precip") ~ "Current",
                             name %in% c("T1","T2to7","T8to14","T15to30") ~ "LagT",
                             name %in% c("P1to30","P31to90","P91to180") ~ "SeasonalP",
                             name %in% c("P181to365","P366to730") ~ "MidtermP",
                             name %in% c("P731to1095","P1096to1460") ~ "LongtermP")) %>%
      filter(!(value==0)) %>%
      group_by(Cat) %>%
      summarise(Q0 = quantile(abs(value),probs=0),
                Q5 = quantile(abs(value),probs=0.05),
                Q10 = quantile(abs(value),probs=0.1),
                Q25 = quantile(abs(value),probs=0.25),
                Q50 = quantile(abs(value),probs=0.5),
                Q75 = quantile(abs(value),probs=0.75),
                Q90 = quantile(abs(value),probs=0.90),
                Q95 = quantile(abs(value),probs=0.95),
                Q100 = quantile(abs(value),probs=1),
                Mean = mean(abs(value)),
                SD = sd(value)) %>%
      ungroup() %>%
      mutate(Site = Site,
             Cluster = k)
    
    # Bind into a single df
    NEE.df = rbind(NEE.df,summary)
    # Tidy up
    rm(output,summary)
    
  }
}

# Read in the site information
info = read_csv("inputs/siteinfo.csv",show_col_types = F) %>%
  select(site,MAPrecord,igbp,PPT_VPD_Qle)

# Merge into dataframe
NEE.df = merge(NEE.df,info,by.x = "Site", by.y = "site")

# Load other site info and merge
load("inputs/siteinfo.Rdata")
NEE.df = merge(NEE.df,INFO,by="Site")

# Group PFTs
NEE.df = NEE.df %>% mutate(PFT = case_when(igbp %in% c("dbf","ebf","enf","mf") ~ "Foresty",
                                   igbp %in% c("cro","gra","sav") ~ "Grassy",
                                   igbp %in% c("csh","osh","wsa") ~ "Woody",
                                   igbp %in% c("wet") ~ "Wetland"))

# Define factors to ensure correct plot order
NEE.df$PFT = factor(NEE.df$PFT, levels = c("Foresty","Woody","Grassy","Wetland"))
NEE.df$Cat = factor(NEE.df$Cat,levels = c("Current","LagT","SeasonalP","MidtermP","LongtermP"))

# Create a Budyko df
budyko.df = data.frame("AI" = seq(0,4,0.01))
budyko.df$EI = ((1-exp(-budyko.df$AI))*budyko.df$AI*tanh(1/budyko.df$AI))^(1/2)

# Fit regressions to AI for each timescale
NEE.fit = data.frame("Cat" = NA,
                 "r.squared" = NA,
                 "p.value" = NA,
                 "intercept" = NA,
                 "coefficient" = NA)
NEE.Seasonal = NEE.df[NEE.df$Cat=="SeasonalP",]
NEE.Seasonal.lin.mod = lm(NEE.Seasonal$Q50 ~ NEE.Seasonal$AI)
NEE.fit = rbind(NEE.fit,
                c("SeasonalP",
                  summary(NEE.Seasonal.lin.mod)$r.squared,
                  summary(NEE.Seasonal.lin.mod)$coefficients[2,4],
                  NEE.Seasonal.lin.mod$coefficients[1],
                  NEE.Seasonal.lin.mod$coefficients[2]))


NEE.Midterm = NEE.df[NEE.df$Cat=="MidtermP",]
NEE.Midterm.lin.mod = lm(NEE.Midterm$Q50 ~ NEE.Midterm$AI)
NEE.fit = rbind(NEE.fit,
                c("MidtermP",
                  summary(NEE.Midterm.lin.mod)$r.squared,
                  summary(NEE.Midterm.lin.mod)$coefficients[2,4],
                  NEE.Midterm.lin.mod$coefficients[1],
                  NEE.Midterm.lin.mod$coefficients[2]))


NEE.Longterm = NEE.df[NEE.df$Cat=="LongtermP",]
NEE.Longterm.lin.mod = lm(NEE.Longterm$Q50 ~ NEE.Longterm$AI)
NEE.fit = rbind(NEE.fit,c("LongtermP",
                          summary(NEE.Longterm.lin.mod)$r.squared,
                          summary(NEE.Longterm.lin.mod)$coefficients[2,4],
                          NEE.Longterm.lin.mod$coefficients[1],
                          NEE.Longterm.lin.mod$coefficients[2])) %>%
  na.omit()

# Fit linear regressions between timescales
NEE.SL.lm = summary(lm(NEE.Longterm$Q50 ~ NEE.Seasonal$Q50))
NEE.SM.lm = summary(lm(NEE.Midterm$Q50 ~ NEE.Seasonal$Q50))
NEE.ML.lm = summary(lm(NEE.Longterm$Q50 ~ NEE.Midterm$Q50))

# Calculate standard deviation between timescales
NEE.sd.df = NEE.df %>% filter(Cat %in% c("SeasonalP","MidtermP","LongtermP")) %>%
  group_by(Site,AI) %>%
  summarise(sym_sd = sd(Q50))

NEE.sd.lin.mod = lm(NEE.sd.df$sym_sd ~ NEE.sd.df$AI)
NEE.sd.fit = data.frame("r.squared" = summary(NEE.sd.lin.mod)$r.squared,
                    "p.value" = summary(NEE.sd.lin.mod)$coefficients[2,4],
                    "intercept" = NEE.sd.lin.mod$coefficients[1] ,
                    "coefficient" = NEE.sd.lin.mod$coefficients[2])


NEE.Seasonal.Mean = NEE.df %>% 
                    filter(Cat == "SeasonalP") %>% 
                    select (Q50) %>% 
                    colMeans()
NEE.Midterm.Mean = NEE.df %>% 
                    filter(Cat == "MidtermP") %>% 
                    select (Q50) %>% 
                    colMeans()
NEE.Longterm.Mean = NEE.df %>% 
                    filter(Cat == "LongtermP") %>% 
                    select (Q50) %>% 
                    colMeans()

# Plot
NEE_main_plot = NEE.df %>% 
  filter(Cat %in% c("SeasonalP","MidtermP","LongtermP")) %>%
  arrange(Q50) %>%
  ggplot() +
  geom_vline(xintercept=1,linetype="dashed") +
  geom_point(data = NEE.df[NEE.df$Cat=="SeasonalP",],
             aes(x=AI-0.025,y=EI-0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_point(data = NEE.df[NEE.df$Cat=="MidtermP",],
             aes(x=AI,y=EI+0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_point(data = NEE.df[NEE.df$Cat=="LongtermP",],
             aes(x=AI+0.025,y=EI-0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_text(data=NEE.fit[NEE.fit$Cat=="SeasonalP",],
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
  geom_text(data=NEE.fit[NEE.fit$Cat=="MidtermP",],
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
  geom_text(data=NEE.fit[NEE.fit$Cat=="LongtermP",],
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
            aes(x=2.75,y=1.1,label=paste0("list(Legacy~Std~Dev == ",
                                          round(as.numeric(intercept),2),
                                          " + ",
                                          round(as.numeric(coefficient),2),
                                          " %*% AI, r^2 == ",
                                          round(as.numeric(r.squared),2),
                                          ",p-value == ",
                                          round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_line(data=budyko.df,aes(x=AI,y=EI)) +
  geom_text(aes(x=0.25,y=0.9),label="Energy\nlimited",size=5) +
  geom_text(aes(x=1.5,y=0.1),label="Water\nlimited",size=5) +
  scale_fill_viridis_c(name = "Coefficient Magnitude", 
                       limits = c(0.1,3),
                       breaks = c(0.1,0.5,1,2,3),
                       labels = c(0,0.5,1,2,3),
                       option = "inferno",
                       trans = "log",
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = "top",
                                               barwidth = 15)) +
  scale_shape_manual(name = "Category", 
                     labels = c("Long-term Legacy","Mid-term Legacy","Seasonal Legacy"),
                     values=c(21,22,23),
                     guide = guide_legend(title.position = "top",
                                          nrow = 1,
                                          reverse = TRUE)) +
  scale_y_continuous(limits = c(0,1.3),breaks=c(0,0.5,1),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,4),expand=c(0,0)) +
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


# REPEAT THE ABOVE FOR LATENT HEAT

# Create dataframe to place metric values in
Qle.df = data.frame()

# Load in the Qle metric values
for (Site in sites){
  for (k in seq){
    
    # Metlag cluster and metlag regression
    load(paste0("outputs/kmeans/Qle/",              
                Site,
                "_metlagcluster_metlagregress_",
                k,
                "c.Rdata"))
    
    summary = output$coeff_df %>% 
      mutate(Cat = case_when(name %in% c("LAI_MOD","SWdown","VPD","Wind","Tair","Precip") ~ "Current",
                             name %in% c("T1","T2to7","T8to14","T15to30") ~ "LagT",
                             name %in% c("P1to30","P31to90","P91to180") ~ "SeasonalP",
                             name %in% c("P181to365","P366to730") ~ "MidtermP",
                             name %in% c("P731to1095","P1096to1460") ~ "LongtermP")) %>% 
      filter(!(value==0)) %>%
      group_by(Cat) %>%
      summarise(Q0 = quantile(abs(value),probs=0),
                Q5 = quantile(abs(value),probs=0.05),
                Q10 = quantile(abs(value),probs=0.1),
                Q25 = quantile(abs(value),probs=0.25),
                Q50 = quantile(abs(value),probs=0.5),
                Q75 = quantile(abs(value),probs=0.75),
                Q90 = quantile(abs(value),probs=0.90),
                Q95 = quantile(abs(value),probs=0.95),
                Q100 = quantile(abs(value),probs=1),
                Mean = mean(abs(value)),
                SD = sd(value)) %>%
      ungroup() %>%
      mutate(Site = Site,
             Cluster = k)
    
    Qle.df = rbind(Qle.df,summary)
    rm(output)
    
  }
}


Qle.df = merge(Qle.df,info,by.x = "Site", by.y = "site")
Qle.df = merge(Qle.df,INFO,by="Site")

Qle.df = Qle.df %>% mutate(PFT = case_when(igbp %in% c("dbf","ebf","enf","mf") ~ "Foresty",
                                   igbp %in% c("cro","gra","sav") ~ "Grassy",
                                   igbp %in% c("csh","osh","wsa") ~ "Woody",
                                   igbp %in% c("wet") ~ "Wetland"))

Qle.df$PFT = factor(Qle.df$PFT, levels = c("Foresty","Woody","Grassy","Wetland"))
Qle.df$Cat = factor(Qle.df$Cat,levels = c("Current","LagT","SeasonalP","MidtermP","LongtermP"))

# Fit regressions to AI
Qle.fit = data.frame("Cat" = NA,
                 "r.squared" = NA,
                 "p.value" = NA,
                 "intercept" = NA,
                 "coefficient" = NA)
Qle.Seasonal = Qle.df[Qle.df$Cat=="SeasonalP",]
Qle.Seasonal.lin.mod = lm(Qle.Seasonal$Q50 ~ Qle.Seasonal$AI)
Qle.fit = rbind(Qle.fit,
                c("SeasonalP",
                  summary(Qle.Seasonal.lin.mod)$r.squared,
                  summary(Qle.Seasonal.lin.mod)$coefficients[2,4],
                  Qle.Seasonal.lin.mod$coefficients[1],
                  Qle.Seasonal.lin.mod$coefficients[2]))


Qle.Midterm = Qle.df[Qle.df$Cat=="MidtermP",]
Qle.Midterm.lin.mod = lm(Qle.Midterm$Q50 ~ Qle.Midterm$AI)
Qle.fit = rbind(Qle.fit,
                c("MidtermP",
                  summary(Qle.Midterm.lin.mod)$r.squared,
                  summary(Qle.Midterm.lin.mod)$coefficients[2,4],
                  Qle.Midterm.lin.mod$coefficients[1],
                  Qle.Midterm.lin.mod$coefficients[2]))


Qle.Longterm = Qle.df[Qle.df$Cat=="LongtermP",]
Qle.Longterm.lin.mod = lm(Qle.Longterm$Q50 ~ Qle.Longterm$AI)
Qle.fit = rbind(Qle.fit,
                c("LongtermP",
                  summary(Qle.Longterm.lin.mod)$r.squared,
                  summary(Qle.Longterm.lin.mod)$coefficients[2,4],
                  Qle.Longterm.lin.mod$coefficients[1],
                  Qle.Longterm.lin.mod$coefficients[2])) %>%
  na.omit()

Qle.SL.lm = summary(lm(Qle.Longterm$Q50 ~ Qle.Seasonal$Q50))
Qle.SM.lm = summary(lm(Qle.Midterm$Q50 ~ Qle.Seasonal$Q50))
Qle.ML.lm = summary(lm(Qle.Longterm$Q50 ~ Qle.Midterm$Q50))

Qle.sd.df = Qle.df %>% filter(Cat %in% c("SeasonalP","MidtermP","LongtermP")) %>%
  group_by(Site,AI) %>%
  summarise(sym_sd = sd(Q50))

Qle.sd.lin.mod = lm(Qle.sd.df$sym_sd ~ Qle.sd.df$AI)
Qle.sd.fit = data.frame("r.squared" = summary(Qle.sd.lin.mod)$r.squared,
                    "p.value" = summary(Qle.sd.lin.mod)$coefficients[2,4],
                    "intercept" = Qle.sd.lin.mod$coefficients[1] ,
                    "coefficient" = Qle.sd.lin.mod$coefficients[2])


Qle.Seasonal.Mean = Qle.df %>% filter(Cat == "SeasonalP") %>% select (Q50) %>% colMeans()
Qle.Midterm.Mean = Qle.df %>% filter(Cat == "MidtermP") %>% select (Q50) %>% colMeans()
Qle.Longterm.Mean = Qle.df %>% filter(Cat == "LongtermP") %>% select (Q50) %>% colMeans()

# Plot
Qle_main_plot = Qle.df %>% filter(Cat %in% c("SeasonalP","MidtermP","LongtermP")) %>%
  arrange(Q50) %>%
  ggplot() +
  geom_vline(xintercept=1,linetype="dashed") +
  geom_point(data = Qle.df[Qle.df$Cat=="SeasonalP",],
             aes(x=AI-0.025,y=EI-0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_point(data = Qle.df[Qle.df$Cat=="MidtermP",],
             aes(x=AI,y=EI+0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_point(data = Qle.df[Qle.df$Cat=="LongtermP",],
             aes(x=AI+0.025,y=EI-0.01,fill=Q50,shape=Cat),
             alpha=0.8,
             size=5,
             stroke=1) +
  geom_text(data=Qle.fit[Qle.fit$Cat=="SeasonalP",],
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
  geom_text(data=Qle.fit[Qle.fit$Cat=="MidtermP",],
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
  geom_text(data=Qle.fit[Qle.fit$Cat=="LongtermP",],
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
            aes(x=2.75,y=1.1,label=paste0("list(Legacy~Std~Dev == ",
                                          round(as.numeric(intercept),2),
                                          " + ",
                                          round(as.numeric(coefficient),2),
                                          " %*% AI, r^2 == ",
                                          round(as.numeric(r.squared),2),
                                          ",p-value == ",
                                          round(as.numeric(p.value),2),")")),
            parse=T,
            size = 6) +
  geom_line(data=budyko.df,aes(x=AI,y=EI)) +
  geom_text(aes(x=0.25,y=0.9),label="Energy\nlimited",size=5) +
  geom_text(aes(x=1.5,y=0.1),label="Water\nlimited",size=5) +
  scale_fill_viridis_c(name = "Coefficient Magnitude", 
                       limits = c(0.1,3),
                       breaks = c(0.1,0.5,1,2,3),
                       labels = c(0,0.5,1,2,3),
                       option = "inferno",
                       trans = "log",
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = "top",
                                               barwidth = 15)) +
  scale_shape_manual(name = "Category", 
                     labels = c("Long-term Legacy","Mid-term Legacy","Seasonal Legacy"),
                     values=c(21,22,23),
                     guide = guide_legend(title.position = "top",
                                          nrow = 1,
                                          reverse = TRUE)) +
  scale_y_continuous(limits = c(0,1.3),breaks=c(0,0.5,1),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,4),expand=c(0,0)) +
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

# Load the inset boxplots created with "Figure3_Inset.R"
load("Qle_boxplot_inset.Rdata")
load("NEE_boxplot_inset.Rdata")

# Remove legends and add insets to main plots
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
# Extract the legend
legend = get_legend(NEE_main_plot)

# Save the plot!
png(filename = paste0(k,"_Figure3.png"),
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




