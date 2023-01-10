### Workflow to plot timeseries of coefficients on site-specific scales
### Begin Setup ###
# Tidy up
rm(list=ls())

# Message
Start = Sys.time()
message("Starting workflow at ",Start)

# Make sure directories and libraries are set up
source("workflows/setup.R")
library(RColorBrewer)
library(RcppRoll)
library(colorspace)
library(ggpubr)

### End Setup ###

# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# SET USER INPUTS FOR WORKFLOW

flux = "NEE"
sitecsv = "longsites.csv"

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

# Load sites
sites = read_csv(paste0("inputs/",
                        sitecsv),
                 col_names=FALSE,
                 show_col_types = F) %>% 
  unlist()

for (Site in sites){
  
      # Load Historical model output
  load(paste0("outputs/kmeans/",
              flux,
              "/",              
              Site,
              "_metlagcluster_metlagregress_729c.Rdata"))
  
  coeff_df = output$coeff_df
  
  # Filter for the lagged precipitation
  coeff_df = coeff_df %>%
    filter(name %in% c("P1to30", 
                       "P31to90",
                       "P91to180",
                       "P181to365",
                       "P366to730",
                       "P731to1095",
                       "P1096to1460"))
  
  # Take negative of NEE to produce NEP
  if (flux == "NEE"){coeff_df$value = -coeff_df$value}
  
  # Define the maximum coefficient value seen
  lim = ceiling(max(abs(coeff_df$value),na.rm=T))
  
  # Load other data and manipulate into our requirements
  load("inputs/sitePPTdeficits.Rdata")
  
  PPTdef = PPTdef %>% mutate(Date = as.Date(paste0(PPTdef$YearMon,"-01"),format="%Y-%m-%d")) %>%
    filter(site == Site) %>%
    mutate(ms_deficitmean = roll_sumr(deficitmean, 6, fill = NA),
           ms_deficitmedian = roll_sumr(deficitmedian, 6, fill = NA))
  
  wide_coeff_df = output$wide_coeff_df
  
  PPTSeries = merge(wide_coeff_df[,1:4],PPTdef[,c("Date","Precip","deficitmedian","ms_deficitmedian")], by.x="Date",by.y="Date",all=T) %>%
    fill(Precip,deficitmedian,ms_deficitmedian,.direction="downup") %>%
    mutate(Rainfall=Precip,
           Deficit=deficitmedian,
           RollDeficit=ms_deficitmedian) %>%
    select(-Precip,-deficitmedian,-ms_deficitmedian) %>% 
    pivot_longer(cols="Rainfall":"RollDeficit") %>%
    na.omit()
  
  coeff_df = rbind(coeff_df,PPTSeries)
  
  # Give the coefficients nice names
  coeff_df = coeff_df %>%
    mutate(name = as.character(name),
           name = str_replace_all(name, 
                                  pattern = "to", replacement = " to "),
           name = str_replace_all(name, 
                                  pattern = "P", replacement = ""),
           name=replace(name, name=="Deficit", "Month Deficit"),
           name=replace(name, name=="RollDeficit", "6 Month Deficit"))
  
  # Turn coefficient names into a factor
  coeff_df$name = coeff_df$name %>% 
    factor(levels=c("6 Month Deficit",
                    "Month Deficit",
                    "Rainfall",
                    "1 to 30", 
                    "31 to 90",
                    "91 to 180",
                    "181 to 365",
                    "366 to 730",
                    "731 to 1095",
                    "1096 to 1460"))

  # Create a plot from which to extract the legend
  legend_plot = ggplot() +
    geom_linerange(data=coeff_df[!coeff_df$name %in% c("Month Deficit","6 Month Deficit","Rainfall"),],
                   aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       colour=value),
                   size=10) +
    scale_colour_continuous_diverging(name = "Coefficient Value (-)",
                           palette = "Red-Green",
                           mid = 0,
                           na.value = "white",
                           trans = "pseudo_log",
                           limits = c(-1,1)*lim,
                           guide = guide_colourbar(direction = "horizontal",
                                                   title.position = "top",
                                                   barwidth = 20,
                                                   order = 2)) +
    ggnewscale::new_scale_color() +
    geom_linerange(data=coeff_df[coeff_df$name %in% c("Month Deficit", "6 Month Deficit"),],
                  aes(x=name,
                      ymin=Date,
                      ymax=Date+1,
                      color=value),
                  size=10) +
    scale_color_continuous_diverging(name = "Rainfall\nDeficit (mm)",
                         palette = "Vik",
                         rev = TRUE,
                         mid = 0,
                         guide = guide_colourbar(direction = "horizontal",
                                                 title.position = "top",
                                                 barwidth = 10,
                                                 order = 3)) +
    
    ggnewscale::new_scale_color() +
    geom_linerange(data=coeff_df[coeff_df$name %in% c("Rainfall"),],
                  aes(x=name,
                      ymin=Date,
                      ymax=Date+1,
                      color=value),
                  size=10) +
    scale_color_continuous_sequential(name = "Rainfall (mm)",
                        palette = "Blues 3",
                        rev = TRUE,
                        guide = guide_colourbar(direction = "horizontal",
                                                title.position = "top",
                                                barwidth = 10,
                                                order = 1)) +
    
    geom_vline(aes(xintercept = 3.5), color = "black") +
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],
                            tail(na.omit(coeff_df)$Date,n=1))) +
    scale_x_discrete(drop = FALSE) +
    theme_bw() +
    coord_flip() +
    xlab("Rainfall Lag (Days)") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 20),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.minor.x = element_line(colour = "black", linetype = "1F"),
          legend.background = element_blank(),
          legend.spacing = unit(2,"cm"),
          legend.position = "bottom",
          legend.box.just = "center",
          legend.title.align=0.5,
          legend.box.margin=margin(0,10,0,10),
          legend.text = element_text(angle=-45,hjust = 0,vjust = 1),
          panel.ontop = TRUE,
          panel.background = element_rect(color = NA, fill = NA)) +
    ggtitle(paste0(Site," Lag Weights"),
            paste0(flux," - ",k,"c"))
  
  # Extract the legend
  legend = get_legend(legend_plot)
  
  # Create the plot of lagged driver coefficients
  lag_plot = ggplot() +
    geom_linerange(data=coeff_df[!coeff_df$name %in% c("Month Deficit","6 Month Deficit","Rainfall"),],
                   aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       colour=value),
                   size=10) +
    scale_colour_continuous_diverging(name = "Coefficient Value (-)",
                                      palette = "Red-Green",
                                      mid = 0,
                                      na.value = "white",
                                      trans = "pseudo_log",
                                      limits = c(-1,1)*lim,
                                      guide = "none") +
    
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],
                            tail(na.omit(coeff_df)$Date,n=1))) +
    scale_x_discrete() +
    theme_bw() +
    coord_flip() +
    xlab("Rainfall Lag (Days into Past)") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 20),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.minor.x = element_line(colour = "black", linetype = "1F"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.spacing = unit(2,"cm"),
          legend.position = "bottom",
          legend.box.just = "center",
          legend.title.align=0.5,
          legend.box.margin=margin(0,10,0,10),
          legend.text = element_text(angle=-45,hjust = 0,vjust = 1),
          panel.ontop = TRUE,
          panel.background = element_rect(color = NA, fill = NA))
  
  # Create the plot of the rainfall deficits
  deficit_plot = ggplot() +
    geom_linerange(data=coeff_df[coeff_df$name %in% c("Month Deficit", "6 Month Deficit"),],
                   aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       color=value),
                   size=10) +
    scale_color_continuous_diverging(name = "Rainfall\nDeficit (mm)",
                                     palette = "Vik",
                                     rev = TRUE,
                                     mid = 0,
                                     guide = "none") +
    
    ggnewscale::new_scale_color() +
    geom_linerange(data=coeff_df[coeff_df$name %in% c("Rainfall"),],
                   aes(x=name,
                       ymin=Date,
                       ymax=Date+1,
                       color=value),
                   size=10) +
    scale_color_continuous_sequential(name = "Rainfall (mm)",
                                      palette = "Blues 3",
                                      rev = TRUE,
                                      guide = "none") +
    scale_y_date(date_breaks = "1 year",
                 date_labels = "%b '%y",
                 expand = c(0,0),
                 limits = c(na.omit(coeff_df)$Date[1],
                            tail(na.omit(coeff_df)$Date,n=1))) +
    scale_x_discrete() +
    theme_bw() +
    coord_flip() +
    xlab("") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 20),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.minor.x = element_line(colour = "black", linetype = "1F"),
          legend.background = element_blank(),
          legend.spacing = unit(2,"cm"),
          legend.position = "bottom",
          legend.box.just = "center",
          legend.title.align=0.5,
          legend.box.margin=margin(0,10,0,10),
          legend.text = element_text(angle=-45,hjust = 0,vjust = 1),
          panel.ontop = TRUE,
          panel.background = element_rect(color = NA, fill = NA))
  
  # Arrange the coefficient and deficit plots together
  plot = ggarrange(lag_plot,deficit_plot, nrow = 2, heights = c(7,3), align = "v")
  
  # Save the plot data
  coeffplots = ggarrange(plot,legend,nrow=2,heights = c(5,1),labels = c("(c)"), font.label = list(size = 20), hjust = -0.1, vjust = 1.1)
  save(coeffplots,file = paste0("outputs/",Site,"_Figure5to8Style.Rdata"))
  
  # Save the plot
  png(filename = paste0("images/",Site,"_",flux,"_Figure5to8Style.png"),
      width = 20,
      height = 10,
      units = "in",
      res = 640)
  print(ggarrange(plot,legend,nrow=2,heights = c(5,1)))
  dev.off()
  
}

