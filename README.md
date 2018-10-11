# Soil_Property_Mapping_2d_v_3D
Repository of code for different approaches of predictive mapping of soil properties

Two scripts are provided that develop a full worflow taking field soil observations, extracting environmental covariates, building a random forest regression models, creating prediction and uncertainty rasters, and validating the model. The workflows allow automated prediction of each soil property for what ever depths a user may choose by looping through desired depths. Scripts were developed for both a workflow where the depth of the soil is included as a covariate in model building (3D approach), and for a workflow where a different model is built for each depth (2D) workflow. The scripts include outputs for cross validation accuracy and global prediction interval summarization. The scripts are set up to run pH models as are, but can be modified following internal commenting to creat models for soil organic carbon, % fine sand + very fine sand, and 1:2 H2O soil electrical conductivity as documentation for the following paper to be submitted to the journal Geoderma:

Nauman, T.W., Duniway, M.C., In Prep. Cautionary notes on uncertainty and accuracy patterns for predictive soil property mapping by depth. For: Geoderma.

Files:
Scripts:
QuantRFmodel_2D.R: 2D approach workflow script
QuantRFmodel_3D.R: 3D approach workflow script

Data:
cop_ncss17SOC__covarsc.txt: Training points for soil organic carbon model with environmental covariates already extracted.
cop_ncss17_FS_VFS_pct_covarsc.txt: Training points for % wt fine sand + very fine sand model with environmental covariates already extracted.
cop_ncss17pH_h20_covarsc.txt: Training points for soil 1:1 H2O pH with environmental covariates already extracted.
cop_ncss17wLIMS_ec12_covarsc.txt: Training points for soil 1:2 H2O electrical conductivity with environmental covariates already extracted.

The environmental covariate data used for predictions is quite large and available at the following link:

 https://www.dropbox.com/sh/bdbpzok2rm2sybn/AACqvp1Mh4kn0ckFoql3JIMga?dl=0
 
 Travis Nauman, PhD
 Soil Scientist,
 Moab, UT
 tnauman@usgs.gov,
 naumi421@gmail.com
