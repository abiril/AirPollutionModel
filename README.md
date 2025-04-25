# A Bayesian Multisource Fusion Model for Spatiotemporal PM$2.5 in an Urban Setting

Code and complementary material for the paper 'A Bayesian Multisource Fusion Model for Spatiotemporal PM2.5 in an Urban Setting' submitted to Environmetrics journal.

## Abstract

We present a Bayesian spatiotemporal fusion model to estimate monthly PM$_{2.5}$ concentrations over Greater London (2014 to 2019) at a 1km resolution. The model integrates multiple air pollution data sources with predictive variables, such as vegetation indices and satellite-derived aerosol optical depth. We deal with issues of spatial misalignment and use an upscaling method to predict across the whole grid area. 

Building on current spatiotemporal stochastic partial differential equations (SPDE) models using the integrated nested Laplace approximations (INLA) approach, our method introduces spatially- and temporally-varying coefficients to flexibly calibrate data sources and better capture fine-scale patterns.

We balance model performance and complexity using predictive metrics such as the PMCC and through thorough cross-validation. The best performing model shows excellent fit (R$^2$ = 0.94, PMCC = 185) and strong predictive performance (temporal CV R$^2$ = 0.72, Cov = 0.90; spatial CV R$^2$ = 0.68, Cov = 0.91), outperforming the baseline INLA-SPDE and related models.

This work is particularly valuable for policy assessment and future exposure studies, and model uncertainty can be explored directly from the outputs of the Bayesian modelling approach.

## GitHub Structure

```md

├── 1. DATA                             # Creating, cleaning, joining and extracting all model data and full potential covariates
│  ├── 1-AURN.Rmd                       # Air pollution monitoring site data from the Automatic Urban and Rural Network
│  ├── 2- London Air.Rmd                # Air pollution monitoring site data from the London Air network
│  ├── 3-AURN and LA Data.Rmd           # Joining monitoring site data
│  ├── 4-Climate.Rmd                    # Processing HAD-UK climate data
│  ├── 5-AOD.Rmd                        # Initial processing of NASA MAIAC AOD satellite data with Inverse Distance Weighting for gap filling
│  ├── 6-NO2 OMI.Rmd                    # 
│  ├── 1-AURN.Rmd
│  ├── 1-AURN.Rmd
