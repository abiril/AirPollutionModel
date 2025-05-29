# A Bayesian Multisource Fusion Model for Spatiotemporal PM2.5 in an Urban Setting

Code and complementary material for the paper 'A Bayesian Multisource Fusion Model for Spatiotemporal PM2.5 in an Urban Setting' submitted to Environmetrics journal.

## Abstract

Airborne particulate matter (PM2.5) is a major public health concern in urban environments, where population density and emission sources exacerbate exposure risks. We present a novel Bayesian spatiotemporal fusion model to estimate monthly PM2.5 concentrations over Greater London (2014-2019) at 1km resolution. The model integrates multiple PM2.5 data sources, including outputs from two atmospheric air quality dispersion models and predictive variables, such as vegetation and satellite aerosol optical depth, while explicitly modelling a latent spatiotemporal field. Spatial misalignment of the data is addressed through an upscaling approach to predict across the entire area. 
Building on stochastic partial differential equations (SPDE) within the integrated nested Laplace approximations (INLA) framework, our method introduces spatially- and temporally-varying coefficients to flexibly calibrate datasets and capture fine-scale variability.

Model performance and complexity are balanced using predictive metrics such as the predictive model choice criterion and thorough cross-validation. The best performing model shows excellent fit and solid predictive performance, enabling reliable high-resolution spatiotemporal mapping of PM2.5 concentrations with the associated uncertainty.
Furthermore, the model outputs, including full posterior predictive distributions, can be used to map exceedance probabilities of regulatory thresholds, supporting air quality management and targeted interventions in vulnerable urban areas, as well as providing refined exposure estimates of PM2.5 for epidemiological applications.

## GitHub Structure

```md
┌── 1. DATA                            # Creating, cleaning, joining and extracting all model data and full potential covariates
│  ┌── 1-AURN.Rmd                      # Air pollution monitoring site data from the Automatic Urban and Rural Network (AURN)
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
│  │  ┌── A-All.Rmd                    # Collating all datsets: rasterisation and extract at monitoring site locations
│  │  ├── B-All and Sat.Rmd            # Joining satellite data: Using Inverse Distance Weighting to gap fill, reproject and extract at monitoring site locations
│  │  ├── C-Joining All Data.Rmd       # Joining site data and grid data and creating various different variables
│  │  └── D-Final Data Prep.Rmd        # Creating more variables and scaling covariates
│  ├──Variable Selection
│  │  ┌── Correlation Analysis.Rmd          # Calculating correlations to PM2.5 and seasonality correlations
│  │  ├── Stepwise Forward Selection.Rmd    # Stepwise Forward Model Selection by the PMCC statistic
│  │  ├── Spike-and-Slab Prior.Rmd          # Spike-and-Slab Prior Regression for variable selection
│  │  └── VIF and Autocorrelation.Rmd       # Variance Inflation Factor (VIF) and Autocorrelation for variable selection
├── 2. MODELS
│  ├── Code
│  │  ┌── M0a.R                 # Model 0a - Covariates Only
│  │  ├── M0b.R                 # Model 0b - Covariates only with  NDVI
│  │  ├── M0c.R                 # Model 0c - Spatial Field Only
│  │  ├── M1a.R                 # Model 1a - Baseline model with covariates and spatial field
│  │  ├── M1b.R                 # Model 1b - Baseline model with covariates and spatial field with NDVI
│  │  ├── M2a.R                 # Model 2a - SVC for AQR
│  │  ├── M2b.R                 # Model 2b - SVC for NDVI
│  │  ├── M3a.R                 # Model 3a - TVC for AQR with AR1 Time
│  │  ├── M3b.R                 # Model 3b - TVC for NDVI with AR1 Time
│  │  ├── M4a.R                 # Model 4a - SVC for PCM and TVC for AQR with AR1 Time
│  │  ├── M4b.R                 # Model 4b - SVC for AQR and TVC for NDVI with AR1 Time
│  │  ├── M4c.R                 # Model 4c - SVC for PCM and NDVI, and TVC for AQR with AR1 Time
│  │  ├── M4d.R                 # Model 4d - SVC for PCM and AQR, and TVC for NDVI with AR1 Time
│  │  └── M4b_Laplace.R         # Re-run of Final Model: MODEL 4b - SVC for AQR and TVC for NDVI with AR1 Time with 'Laplace' method
│  ├── ShellFile
│  │  ┌── M0a.sh
│  │  ├── M0b.sh
│  │  ├── M0c.sh
│  │  ├── M1a.sh
│  │  ├── M1b.sh
│  │  ├── M2a.sh
│  │  ├── M2b.sh
│  │  ├── M3a.sh
│  │  ├── M3b.sh
│  │  ├── M4a.sh
│  │  ├── M4b.sh
│  │  ├── M4c.sh
│  │  ├── M4d.sh
│  │  └── M4b_Laplace.sh
├── 3. RESULTS
│  ├── Model Outputs and Cross-Validation
│  │  ┌── M0a.ROut
│  │  ├── M0b.ROut
│  │  ├── M0c.ROut
│  │  ├── M1a.ROut
│  │  ├── M1b.ROut
│  │  ├── M2a.ROut
│  │  ├── M2b.ROut
│  │  ├── M3a.ROut
│  │  ├── M3b.ROut
│  │  ├── M4a.ROut
│  │  ├── M4b.ROut
│  │  ├── M4c.ROut
│  │  ├── M4d.ROut
│  │  └── M4b_Laplace.ROut
│  ├── Final Model (To be uploaded)
│  │  ┌── INLA_OUTPUT_M4b.ROut   # Final Model INLA Output File
│  │  ├── OUTPUT_M4b.RData       # Final Model INLA Output Data
│  │  ├── TINE_CV_M4b.RData      # Final Model Time Cross-Validation INLA Output Data
│  │  ├── SPATIAL_CV_M4b.RData   # Final Model Spatial Cross-Validation INLA Output Data
│  │  ├── ALL_M4b.RData          # All Model Data with Outputs: Dataframe of all data with final model output
└──└──└── FINAL_PM25.RData       # Final predicted PM25 data: observed, predicted, sd, quantiles, etc



