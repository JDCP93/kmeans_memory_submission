## NICE PLOT OF AI/EI/PFT AND METRIC IMPROVEMENT

# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")
library(Hmisc)

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

# Create dataframe to place metric values in
df = data.frame("Site" = rep(sites,
                             each=20),
                "Flux" = rep(c("NEE","Qle"),each=10,times = length(sites)),
                "Clustering" = rep(c("Met","Met+Lag"),
                                   each = 5,
                                   times = 2*length(sites)),
                "Regression" = rep(c("Met","Met+Lag"),
                                   each = 5,
                                   times = 2*length(sites)),
                "Metric" = rep(c("r^2",
                                 "Mean Bias Error",
                                 "Normalised Mean Error",
                                 "Std. Dev. Difference",
                                 "Correlation Coeff."),
                               times = 2*length(sites)),
                "Value" = NA)

# Load in the metric values
for (Site in sites){
  for (flux in c("NEE","Qle")){
    
    # Load met cluster and met regression output (i.e. Instantaneous model)
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metcluster_metregress_729c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Flux == flux &
               df$Clustering == "Met" & 
               df$Regression == "Met"] = unlist(output[c(3:7)])
    rm(output)
    
    # Load met+lag cluster and met+lag regression output (i.e. Historical model)
    load(paste0("outputs/kmeans/",
                flux,
                "/",              
                Site,
                "_metlagcluster_metlagregress_729c.Rdata"))
    
    df$Value[df$Site == Site & 
               df$Flux == flux &
               df$Clustering == "Met+Lag" & 
               df$Regression == "Met+Lag"] = unlist(output[c(3:7)])
    rm(output)
    
  }
}

# Load other data and merge into the dataframe
info = read_csv("inputs/siteinfo.csv") %>%
  select(site,MAPrecord,igbp,PPT_VPD_Qle)

df = merge(df,info,by.x = "Site", by.y = "site")

load("inputs/siteinfo.Rdata")
df = merge(df,INFO,by="Site")

# Turn model and metrics into factors so they appear in correct order
df$Clustering = factor(df$Clustering,
                       levels = c("Met","Met+Lag"))

df$Regression = factor(df$Regression,
                       levels = c("Met","Met+Lag"))

df$Metric = factor(df$Metric, 
                   levels = c("Normalised Mean Error",
                              "Std. Dev. Difference",
                              "Mean Bias Error",
                              "r^2",
                              "Correlation Coeff."),
                   labels = c(expression("Normalised~Mean~Error"),
                              expression("Std.~Dev.~Difference"),
                              expression("Mean~Bias~Error"),
                              expression("r^2"),
                              expression("Correlation~Coeff.")))

df = df %>% mutate(igbp = replace(igbp,igbp=="wsa","sav")) %>%
            mutate(igbp = replace(igbp,igbp%in% c("csh","osh"),"shr"))

df = df %>% mutate(PFT = case_when(igbp %in% c("dbf","ebf","enf","mf") ~ "Foresty",
                                   igbp %in% c("cro","gra","sav") ~ "Grassy",
                                   igbp %in% c("csh","osh","wsa") ~ "Woody",
                                   igbp %in% c("wet") ~ "Wetland"))

df$PFT = factor(df$PFT, levels = c("Foresty","Woody","Grassy","Wetland"))


# Calculate the relative and absolute improvement between the two models
metrics = df %>% pivot_wider(names_from = Clustering:Regression, values_from = Value) %>%
  mutate("RelImpMet" = case_when(Metric == "Normalised~Mean~Error" ~ (Met_Met-`Met+Lag_Met+Lag`)/Met_Met,
                                 TRUE ~ (`Met+Lag_Met+Lag`-Met_Met)/Met_Met),
         "AbsImpMet" = case_when(Metric =="Normalised~Mean~Error" ~ (Met_Met-`Met+Lag_Met+Lag`),
                                 TRUE ~ `Met+Lag_Met+Lag`-Met_Met)) %>%
  filter(Metric %in% c("r^2","Normalised~Mean~Error")) %>%
  arrange(MAPrecord)

# Create a Budyko df
budyko.df = data.frame("AI" = seq(0,4,0.01))
budyko.df$EI = ((1-exp(-budyko.df$AI))*budyko.df$AI*tanh(1/budyko.df$AI))^(1/2)

# Calculate the r squared for RelImp ~ AI
fit.Rel = data.frame("Metric" = NA,
                 "Flux" = NA,
                 "r.squared" = NA,
                 "p.value" = NA,
                 "intercept" = NA,
                 "coefficient" = NA)
NME.NEE = metrics[metrics$Metric=="Normalised~Mean~Error" & metrics$Flux=="NEE",]
NME.NEE.lin.mod = lm(NME.NEE$RelImpMet ~ NME.NEE$AI)
fit.Rel = rbind(fit.Rel,
                c("Normalised~Mean~Error",
                  "NEP",
                  summary(NME.NEE.lin.mod)$r.squared,
                  summary(NME.NEE.lin.mod)$coefficients[2,4],
                  NME.NEE.lin.mod$coefficients[1],
                  NME.NEE.lin.mod$coefficients[2]))

NME.Qle = metrics[metrics$Metric=="Normalised~Mean~Error" & metrics$Flux=="Qle",]
NME.Qle.lin.mod = lm(NME.Qle$RelImpMet ~ NME.Qle$AI)
fit.Rel = rbind(fit.Rel,
                c("Normalised~Mean~Error",
                  "LE",
                  summary(NME.Qle.lin.mod)$r.squared,
                  summary(NME.Qle.lin.mod)$coefficients[2,4],
                  NME.Qle.lin.mod$coefficients[1],
                  NME.Qle.lin.mod$coefficients[2]))

R2.NEE = metrics[metrics$Metric=="r^2" & metrics$Flux=="NEE",]
R2.NEE.lin.mod = lm(R2.NEE$RelImpMet ~ R2.NEE$AI)
fit.Rel = rbind(fit.Rel,
                c("r^2",
                  "NEP",
                  summary(R2.NEE.lin.mod)$r.squared,
                  summary(R2.NEE.lin.mod)$coefficients[2,4],
                  R2.NEE.lin.mod$coefficients[1],
                  R2.NEE.lin.mod$coefficients[2]))

R2.Qle = metrics[metrics$Metric=="r^2" & metrics$Flux=="Qle",]
R2.Qle.lin.mod = lm(R2.Qle$RelImpMet ~ R2.Qle$AI)
fit.Rel = rbind(fit.Rel,
                c("r^2",
                  "LE",
                  summary(R2.NEE.lin.mod)$r.squared,
                  summary(R2.Qle.lin.mod)$coefficients[2,4],
                  R2.NEE.lin.mod$coefficients[1],
                  R2.NEE.lin.mod$coefficients[2])) %>%
  na.omit()

# Calculate the r squared for RelImp ~ EI
EI.fit.Rel = data.frame("Metric" = NA,
                     "Flux" = NA,
                     "r.squared" = NA,
                     "p.value" = NA,
                     "intercept" = NA,
                     "coefficient" = NA)
EI.NME.NEE = metrics[metrics$Metric=="Normalised~Mean~Error" & metrics$Flux=="NEE",]
EI.NME.NEE.lin.mod = lm(EI.NME.NEE$RelImpMet ~ EI.NME.NEE$EI)
EI.fit.Rel = rbind(EI.fit.Rel,
                   c("Normalised~Mean~Error",
                     "NEP",
                     summary(EI.NME.NEE.lin.mod)$r.squared,
                     summary(EI.NME.NEE.lin.mod)$coefficients[2,4],
                     EI.NME.NEE.lin.mod$coefficients[1],
                     EI.NME.NEE.lin.mod$coefficients[2]))

EI.NME.Qle = metrics[metrics$Metric=="Normalised~Mean~Error" & metrics$Flux=="Qle",]
EI.NME.Qle.lin.mod = lm(EI.NME.Qle$RelImpMet ~ EI.NME.Qle$EI)
EI.fit.Rel = rbind(EI.fit.Rel,
                   c("Normalised~Mean~Error",
                     "LE",
                     summary(EI.NME.Qle.lin.mod)$r.squared,
                     summary(EI.NME.Qle.lin.mod)$coefficients[2,4],
                     EI.NME.Qle.lin.mod$coefficients[1],
                     EI.NME.Qle.lin.mod$coefficients[2]))

EI.R2.NEE = metrics[metrics$Metric=="r^2" & metrics$Flux=="NEE",]
EI.R2.NEE.lin.mod = lm(EI.R2.NEE$RelImpMet ~ EI.R2.NEE$EI)
EI.fit.Rel = rbind(EI.fit.Rel,
                   c("r^2",
                     "NEP",
                     summary(EI.R2.NEE.lin.mod)$r.squared,
                     summary(EI.R2.NEE.lin.mod)$coefficients[2,4],
                     EI.R2.NEE.lin.mod$coefficients[1],
                     EI.R2.NEE.lin.mod$coefficients[2]))

EI.R2.Qle = metrics[metrics$Metric=="r^2" & metrics$Flux=="Qle",]
EI.R2.Qle.lin.mod = lm(EI.R2.Qle$RelImpMet ~ EI.R2.Qle$EI)
EI.fit.Rel = rbind(EI.fit.Rel,
                   c("r^2",
                     "LE",
                     summary(EI.R2.NEE.lin.mod)$r.squared,
                     summary(EI.R2.Qle.lin.mod)$coefficients[2,4],
                     EI.R2.NEE.lin.mod$coefficients[1],
                     EI.R2.NEE.lin.mod$coefficients[2])) %>%
  na.omit()

# Get some summary statistics
summaries = metrics %>% 
      group_by(Flux,Metric) %>%
      summarise_at(c("Met_Met","Met+Lag_Met+Lag","RelImpMet","AbsImpMet"),
                   list(min = min, mean = mean,median = median,max = max))

# Perform corrleation tests
relImpCorr = metrics %>% 
      select(Site,Flux,Metric,RelImpMet) %>% 
      pivot_wider(names_from = c("Flux","Metric"), values_from = "RelImpMet")
rcorr(as.matrix(relImpCorr[,2:5]))

cor.test(metrics$Met_Met[metrics$Flux=="NEE" & metrics$Metric == "r^2"],
         metrics$Met_Met[metrics$Flux=="Qle" & metrics$Metric == "r^2"])
cor.test(metrics$Met_Met[metrics$Flux=="NEE" & metrics$Metric == "Normalised~Mean~Error"],
         metrics$Met_Met[metrics$Flux=="Qle" & metrics$Metric == "Normalised~Mean~Error"])

# Rename fluxes
metrics$Flux[metrics$Flux=="NEE"] = "NEP"
metrics$Flux[metrics$Flux=="Qle"] = "LE"

metrics$Flux = factor(metrics$Flux,levels=c("NEP","LE"))

# Create the plot
plot = metrics %>%
  arrange(RelImpMet) %>%
  ggplot() +
  geom_vline(xintercept=1,
             linetype="dashed") +
  geom_point(aes(x=AI,y=EI,fill=RelImpMet),
             shape=23,
             alpha=0.9,
             size=5,
             stroke=1) +
  geom_line(data=budyko.df,
            aes(x=AI,y=EI)) +
  geom_text(data=fit.Rel,
            aes(x=3.25,y=0.275,label=paste0("Imp == ",
                                            signif(as.numeric(intercept),2),
                                            " + ",
                                            signif(as.numeric(coefficient),2),
                                            " %*% AI")),
            parse=T,
            size = 20*0.36) +
  geom_text(data=fit.Rel,
            aes(x=3.25,y=0.2,label=paste0("r^2 == ",
                                          signif(as.numeric(r.squared),2))),
            parse=T,
            size = 20*0.36) +
  geom_text(data=fit.Rel,
            aes(x=3.25,y=0.125,label=paste0("p-value == ",
                                            signif(as.numeric(p.value),2))),
            parse=T,
            size = 20*0.36) +
  geom_text(aes(x=0.25,y=1),
            label="Energy\nlimited",
            size=20*0.36) +
  geom_text(aes(x=1.75,y=0.1),
            label="Water\nlimited",
            size=20*0.36) +
  scale_fill_viridis_c(name = "Relative Improvement", 
                       limits = c(0,0.6),
                       breaks = c(0,0.2,0.4,0.6),
                       option = "inferno",
                       labels = scales::percent,
                       direction = -1,
                       trans = "sqrt",
                       guide = guide_colourbar(direction = "vertical",
                                               title.position = "right",
                                               barheight = 20)) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  facet_grid(factor(Flux,levels=c('NEP','LE'))~Metric,labeller = label_parsed) +
  ylab("Evaporative Index (AET/PPT)") +
  xlab("Aridity Index (PET/PPT)") +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = "right",
        legend.box.just = "center",
        legend.title.align=0.5,
        legend.box.margin=margin(5,0,5,0),
        legend.title = element_text(angle = -90))

ggsave(filename = "images/Figure2.png",
       plot = plot,
       device = "png",
       dpi = 320,
       width = 16, height = 12)





