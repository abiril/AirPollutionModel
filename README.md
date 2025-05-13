# A Bayesian Multisource Fusion Model for Spatiotemporal PM2.5 in an Urban Setting

Code and complementary material for the paper 'A Bayesian Multisource Fusion Model for Spatiotemporal PM2.5 in an Urban Setting' submitted to Environmetrics journal.

## Abstract

We present a Bayesian spatiotemporal fusion model to estimate monthly PM$_{2.5}$ concentrations over Greater London (2014 to 2019) at a 1km resolution. The model integrates multiple air pollution data sources with predictive variables, such as vegetation indices and satellite-derived aerosol optical depth. We deal with issues of spatial misalignment and use an upscaling method to predict across the whole grid area. 

Building on current spatiotemporal stochastic partial differential equations (SPDE) models using the integrated nested Laplace approximations (INLA) approach, our method introduces spatially- and temporally-varying coefficients to flexibly calibrate data sources and better capture fine-scale patterns.

We balance model performance and complexity using predictive metrics such as the PMCC and through thorough cross-validation. The best performing model shows excellent fit (R$^2$ = 0.94, PMCC = 185) and strong predictive performance (temporal CV R$^2$ = 0.72, Cov = 0.90; spatial CV R$^2$ = 0.68, Cov = 0.91), outperforming the baseline INLA-SPDE and related models.

This work is particularly valuable for policy assessment and future exposure studies, and model uncertainty can be explored directly from the outputs of the Bayesian modelling approach.

## GitHub Structure

```md
├── 1. DATA                            # Creating, cleaning, joining and extracting all model data and full potential covariates
│  ├── 1-AURN.Rmd                      # Air pollution monitoring site data from the Automatic Urban and Rural Network (AURN)
│  ├── 2- London Air.Rmd               # Air pollution monitoring site data from the London Air network
│  ├── 3-AURN and LA Data.Rmd          # Joining monitoring site data
│  ├── 4-Climate.Rmd                   # Processing HAD-UK climate data
│  ├── 5-AOD.Rmd                       # (OLD) Initial processing of NASA MAIAC AOD satellite data with Inverse Distance Weighting for gap filling
│  ├── 6-NO2 OMI.Rmd                   # (OLD) Initial processing of NASA OMI/AURA NO2 column data with Inverse Distance Weighting for gap filling
│  ├── 7-Boundary Layer Height.Rmd     # Processing and reprojection of Boundary Layer Height (BLH) from ERA5 Land
│  ├── 8-PCM.Rmd                       # Processing of Pollution Climate Mapping (PCM) model data
│  ├── 9-NDVI.Rmd                      # Processing MODIS NDVI and reprojection with Inverse Distance Weight
│  ├── 10-Rural-Urban.Rmd              # Rural - Urban classification of UK output areas
│  ├── 11-GHSL.Rmd                     # Extraction of GHSL data for study domain and reprojection
│  ├── 12-Population.Rmd               # World Pop population data processing
│  ├── 13-Wind and Evapo.Rmd           # Processing and projection of Wind and Evapotranspiration/Evaporation from ERA5 Land
│  ├── 14-AQUM.Rmd                     # Processing of AQUM data and reprojection by Inverse Distance Weighting
│  ├── 15-Corine.Rmd                   # Processing for CORINE Land Cover data and categorisation 
│  ├── 16-LAEI.Rmd                     # Processing London Atmospheric Emissions Inventory
│  ├── 17-OSM_Roads.Rmd                # Download roads data using 'osmdata' and creating measures of nearby roads
│  ├── 17-Road_Density-.Rmd            # Processing GRIP global road density data
│  ├── 18-Greenspace.QGIS              # Processing OSM greenspace
│  ├── 18-Greenspace.Rmd               # Processing greenspace data from QGIS
│  ├── 19-Roads.Rmd                    # Download roads using 'osmdata', creating distance to road variables
│  ├── 20-Met Office Reanalysis.Rmd    # Processing the Met Office Reanalysis of the AQUM model
│  ├── Data Processing
│  │  ├── A-All.Rmd                    # Collating all datsets: rasterisation and extract at monitoring site locations
│  │  ├── B-All and Sat.Rmd            # Joining satellite data: Using Inverse Distance Weighting to gap fill, reproject and extract at monitoring site locations
│  │  ├── C-Joining All Data.Rmd       # Joining site data and grid data and creating various different variables
│  │  └── D-Final Data Prep.Rmd        # Creating more variables and scaling covariates
│  ├──Variable Selection
│  │  ├── Correlation Analysis.Rmd
│  │  ├── Stepwise Forward Selection.Rmd
│  │  └── VIF and Autocorrelation.Rmd
├── 2. MODELS
│  ├── Code
│  │  ├── Covariates.R                 # Model 0a
│  │  ├── Covariates_wNDVI.R           # Model 0b - with NDVI
│  │  ├── SpatialField.R               # Model 0c
│  │  ├── Baseline.R                   # Model 1a
│  │  ├── Baseline_wNDVI.R             # Model 1b - with NDVI
│  │  ├── SVC_AQR.R                    # Model 2a - SVC for AQR
│  │  ├── SVC_NDVI.R                   # Model 2b - SVC for NDVI
│  │  ├── TVC_AQR.R                    # Model 3a - TVC for AQR with AR1 Time
│  │  ├── TVC_NDVI.R                   # Model 3b - TVC for NDVI with AR1 Time
│  │  ├── SVC_PCM_TVC_AQR.R            # Model 4a - SVC for PCM and TVC for AQR with AR1 Time
│  │  ├── SVC_AQR_TVC_NDVI.R           # Model 4b - SVC for AQR and TVC for NDVI with AR1 Time
│  │  ├── SVC_PCM_NDVI_TVC_AQR.R       # Model 4c - SVC for PCM and NDVI, and TVC for AQR with AR1 Time
│  │  ├── SVC_PCM_AQR_TVC_NDVI.R       # Model 4d - SVC for PCM and AQR, and TVC for NDVI with AR1 Time
│  │  └── FINAL_SVC_AQR_TVC_NDVI.R     # Re-run of Final Model: MODEL 4b - SVC for AQR and TVC for NDVI with AR1 Time with 'Laplace' method
│  ├── ShellFile
│  │  ├── Covariates.sh
│  │  ├── Covariates_wNDVI.sh
│  │  ├── SpatialField.sh
│  │  ├── Baseline.sh                  
│  │  ├── Baseline_wNDVI.sh 
│  │  ├── SVC_AQR.sh
│  │  ├── SVC_NDVI.sh 
│  │  ├── TVC_AQR.sh 
│  │  ├── TVC_NDVI.sh
│  │  ├── SVC_PCM_TVC_AQR.sh 
│  │  ├── SVC_AQR_TVC_NDVI.sh 
│  │  ├── SVC_PCM_NDVI_TVC_AQR.sh  
│  │  ├── SVC_PCM_AQR_TVC_NDVI.sh 
│  │  └── FINAL_SVC_AQR_TVC_NDVI.sh 
├── 3. RESULTS
│  ├── Model Outputs and Cross-Validation
│  │  ├── Covariates.ROut
│  │  ├── Covariates_wNDVI.ROut
│  │  ├── SpatialField.ROut
│  │  ├── Baseline.ROut                  
│  │  ├── Baseline_wNDVI.ROut 
│  │  ├── SVC_AQR.ROut
│  │  ├── SVC_NDVI.ROut 
│  │  ├── TVC_AQR.ROut 
│  │  ├── TVC_NDVI.ROut
│  │  ├── SVC_PCM_TVC_AQR.ROut 
│  │  ├── SVC_AQR_TVC_NDVI.ROut 
│  │  ├── SVC_PCM_NDVI_TVC_AQR.ROut  
│  │  ├── SVC_PCM_AQR_TVC_NDVI.ROut 
│  │  └── FINAL_SVC_AQR_TVC_NDVI.ROut 
│  ├── Final Model
│  │  ├── INLA_OUTPUT_SVC_AQR_TVC_NDVI.ROut   # Final Model INLA Output File
│  │  ├── OUTPUT_SVC_AQR_TVC_NDVI.RData       # Final Model INLA Output Data
│  │  ├── TINE_CV_SVC_AQR_TVC_NDVI.RData      # Final Model Time Cross-Validation INLA Output Data
│  │  ├── SPATIAL_CV_SVC_AQR_TVC_NDVI.RData   # Final Model Spatial Cross-Validation INLA Output Data
│  │  ├── ALL_SVC_AQR_TVC_NDVI.RData          # ALL MODEL DATA with Outputs: Dataframe of all data with final model output
│  │  └── FINAL_PM25.RData                    # Final predicted PM25 data: observed, predicted, sd, quantiles, etc



