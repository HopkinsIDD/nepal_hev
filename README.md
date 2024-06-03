# nepal_hev
This repository contains the code and serosurvey data used to map the predicted seroprevalence of hepatitis E virus across Nepal in 2021, as described in the associated manuscript (add link when preprint/paper available). The code is published exactly as it was used to generate the figures and estimates presented in the manuscript, and therefore depends on the household-level GPS data which we cannot make public. However, we do include here the deidentified serosurvey data, alongside the spatial covariates and the administrative boundaries used.

The code folder includes the following R scripts:
1. "dependencies.R" which includes R packages and functions used in the analysis.
2. "get_rasters.R" which extracts, tranforms, and plots the publicaly available spatial data on population density, travel time to nearest city and elevation, used in the spatial regression.
3. "generate_grid_data.R" which generates the dataset of covariates and merges these with the serosurvey data for the spatial analyses.
4. "implement_spatial_analysis.R" which fits hierarchical logistic spatial regression models to the spatial data within a Bayesian framework using Laplace approximation as implemented in the R package "INLA".

The data folder includes:
1. The shapefiles used to specify the adminsitrative boundaries used in mapping.
2. The spatial covariates (population density, elevation and travel time to the nearest city)
3. The deidentified serosurvey data
4. The generated data - the predicted seroprevalence across Nepal.
