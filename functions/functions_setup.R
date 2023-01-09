################################################################################

### Setup functions to make sure directories and similar exist

################################################################################

################################################################################
directorysetup = function(){
  # Function to create all the high-level directories
  dir.create(file.path("outputs/"), showWarnings = FALSE)
  dir.create(file.path("outputs/kmeans/"), showWarnings = FALSE)
  dir.create(file.path("outputs/kmeans/clusterings"), showWarnings = FALSE)
  dir.create(file.path("outputs/kmeans/clusterings/LOO"), showWarnings = FALSE)
  dir.create(file.path("inputs/"), showWarnings = FALSE)
  dir.create(file.path("inputs/FLUXNET_processed/"), showWarnings = FALSE)
  dir.create(file.path("inputs/ELI/"), showWarnings = FALSE)
  dir.create(file.path("inputs/images/"), showWarnings = FALSE)
}


################################################################################
outputfolders = function(flux){
  # Function to create output folders for a certain model output
  # Inputs:
  #         flux = flux being modelled
  dir.create(file.path(paste0("outputs/kmeans/",flux)), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/images")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/images/coeff")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/images/coeff/abs")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/images/heatmaps")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/images/timeseries")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO/images")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO/images/coeff")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO/images/coeff/abs")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO/images/heatmaps")), showWarnings = FALSE)
  dir.create(file.path(paste0("outputs/kmeans/",flux,"/LOO/images/timeseries")), showWarnings = FALSE)
  dir.create(file.path("images/"), showWarnings = FALSE)
  dir.create(file.path("images/Daily"), showWarnings = FALSE)
  dir.create(file.path("images/Monthly"), showWarnings = FALSE)
  dir.create(file.path("images/Yearly"), showWarnings = FALSE)
  dir.create(file.path("images/ModelGrid"), showWarnings = FALSE)
  dir.create(file.path("images/ObsPred"), showWarnings = FALSE)
  dir.create(file.path("images/YearAvgs"), showWarnings = FALSE)
  dir.create(file.path("images/Rolltime"), showWarnings = FALSE)
}