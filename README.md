# kmeans_memory_submission

Repo to store example code for the submission of Cranko Page et al., 2022 - currently submitted to JGR: Biogeosciences.

## Repo Structure

* Workflows - Code here is used to prepare input files and run the model
  - cluster.R - Clusters the site data (the "k-means clustering" part)
  - regression.R - performs the robust regressions within each cluster (the "plus regression" part)
* Functions - Code here is called by scripts in workflows/, in an attempt to compartmentalise and improve code readability
* Analysis - Scripts used to create Figure 3 from the paper

Please be aware that scripts are hardcoded to point at certain directories - you may need to amend the scripts to suit your local workspace.

### Admin

Apologies for the standard of the code.

Please contact me if you have any questions, issues, or requests.
Code to recreate other figures will be shared upon request, and added here slowly, but requires some clean up to be suitable for public consumption. 
