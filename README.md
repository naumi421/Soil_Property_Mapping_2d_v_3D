# Soil_Property_Mapping_2d_v_3D
Repository of R statistical programing language code for different approaches of predictive mapping of soil properties.

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

Two scripts are provided that develop a full worflow taking field soil observations, extracting environmental covariates, building a random forest regression models, creating prediction and uncertainty rasters, and validating the model. The workflows allow automated prediction of each soil property for what ever depths a user may choose by looping through desired depths. Scripts were developed for both a workflow where the depth of the soil is included as a covariate in model building (3D approach), and for a workflow where a different model is built for each depth (2D) workflow. The scripts include outputs for cross validation accuracy and global prediction interval summarization. The scripts are set up to run pH models as are, but can be modified following internal commenting to creat models for soil organic carbon, % fine sand + very fine sand, and 1:2 H2O soil electrical conductivity as documentation for the following paper submitted to the journal Geoderma:

Nauman, T.W., Duniway, M.C., In Revision. Relative prediction intervals reveal larger uncertainty in 3D approaches to predictive digital soil mapping of soil properties with legacy data. For: Geoderma.

Files:

Scripts:

QuantRFmodel_2D.R: 2D approach workflow script

QuantRFmodel_3D.R: 3D approach workflow script

Data:

cop_ncss17SOC__covarsc.txt: Training points for soil organic carbon model with environmental covariates already extracted.

cop_ncss17_FS_VFS_pct_covarsc.txt: Training points for % wt fine sand + very fine sand model with environmental covariates already extracted.

cop_ncss17pH_h20_covarsc.txt: Training points for soil 1:1 H2O pH with environmental covariates already extracted.

cop_ncss17wLIMS_ec12_covarsc.txt: Training points for soil 1:2 H2O electrical conductivity with environmental covariates already extracted.

Example maps for review produced from these scripts for the paper are available at the following repository:

Version 1: http://doi.org/10.5281/zenodo.1458273

Version 2: http://doi.org/10.5281/zenodo.2545882
   - Includes updated 2D EC RPI raster after correcting an error with the 95% interquantile range used in calculations


 
 Travis Nauman, PhD,
 Soil Scientist,
 Moab, UT
 tnauman@usgs.gov,
 naumi421@gmail.com
