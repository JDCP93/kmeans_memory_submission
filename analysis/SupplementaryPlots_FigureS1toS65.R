
rm(list=ls())

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
  
  
  load(paste0(Site,"_timeseries_plot.Rdata"))
  
  load(paste0(Site,"_coeffplots.Rdata"))
  
  
  
  ggarrange(plots[[4]],coeffplots, nrow = 2, align = "h", heights=c(2,3))

  
  png(filename = paste0(Site,"_supplement_",flux,"_729c.png"),
      width = 20,
      height = 22,
      units = "in",
      res = 320)
  print(ggarrange(plots[[4]]+
                        theme(plot.background=element_rect(color="black")),
                  coeffplots+
                        theme(plot.background=element_rect(color="black")), 
                    nrow = 2, 
                    align = "h", 
                    heights=c(2,3)))
  dev.off()
  
}