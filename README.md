# FUSE-MMComparison-paper

This repository allows for the reproduction of results from the paper:

**"Comparing multi-model mosaic and multi-model combination methods to simulate streamflow across the contiguous USA"** 
*Hydrology and Earth System Sciences, HESS*
*Cyril Thébault, Wouter J. M. Knoben, Nans Addor, Andrew J. Newman, Diana Spieler, Nicolás A. Vásquez, Yalan Song, Gaby J. Gründemann, Shaun Carney, Mukesh Kumar, Katie van Werkhoven, Chaopeng Shen, Andrew W. Wood, and Martyn P. Clark (2026)* 

## Repository Structure

The repository is organized as follows:

- **`00_DATA/`**: Contains the evaluation metrics and sampling uncertainty from individual FUSE models, multi-model mosaics, static combinations, and dynamic combination approaches. The lists of FUSE decisions, parameters for the dynamic combination and catchments informations are also provided here. The data can be download directly from: (https://www.hydroshare.org/resource/53063da8fc894149b5b25f3201d1e8e8/)

- **`01_FUSEscripts/`**: Includes scripts to run individual FUSE models. These scripts help set up, execute and evaluate the 78 hydrological models used in the study.

- **`02_MOSAscripts/`**: Provides scripts to reproduce the multi-model mosaics, which was derived from the paper "Technical note: How many models do we need to simulate hydrologic processes across large geographical domains?" (Knoben et al., 2025).

- **`03_SCscripts/`**: Contains scripts to conduct and evaluate the static combination approaches.

- **`04_DCscripts/`**: Contains scripts to run and evaluate the dynamic combination approach.

- **`05_SUscripts/`**: Contains scripts to calculate sampling uncertainty for the tested multi-model approaches.

- **`80_Visualisation/`**: Includes scripts to generate the figures used in the manuscript, leveraging the data stored in `00_DATA/` (outputs of the previous stages).

- **`99_Figures/`**: Include the figures generated with `80_Visualisation/`.

- **`Shp/`**: Include USA boundaries shapefile (from “North American Atlas - Political Boundaries” (Commission for Environmental Cooperation, 2022)) and catchment shapefiles (derived from CAMELS dataset (Addor et al, 2017)) used in `04_Visualisation/`. The folder can be download directly from: (https://www.hydroshare.org/resource/53063da8fc894149b5b25f3201d1e8e8/)

Metrics.R file includes various functions to calculate metrics for streamflow evaluation. 

InstallRPackages.R installs the various R packages used in this work.

## Reproducing the Results

To reproduce the results presented in the paper, follow these steps:

1. **Run FUSE models**  
   The scripts in `01_FUSEscripts/` generate individual model outputs.

2. **Create multi-model mosaics**  
   The script in `02_MOSAscripts/` creates the mosaic approaches, based on performance and performance-equivalence.

3. **Perform static combinations**  
   The scripts in `03_SCscripts/` generate combinations of individual model outputs of two or three models with a simple average.

4. **Apply dynamic combination approach**  
   The scripts in `04_DCscripts/` implement the time-varying model combination methodology.

5. **Calculate sampling uncertainty**  
   The script in `05_SUscripts/` calculates the sampling uncertainty surrounding performance score for each multi-model approach.

6. **Generate visualizations**  
   Run the scripts in `80_Visualisation/` to reproduce the figures presented in the manuscript.

## Requirements

To run the scripts, ensure you have the following dependencies installed:

- R (version 4.3.1 used here)  
- Required R packages: run InstallRPackages.R script
- On a HPC: load the different module needed to run FUSE and to install R packages (e.g. gcc, R, netcdf, proj, gdal, geos, hdf5, openblas -- HPC dependant)

## Citation

If you use this repository or its outputs in your research, please cite:

> Thébault, C., Knoben W. J. M., Addor, N., Newman, A. J., Spieler, D., Vásquez, N. A., Song, Y., Gründemann, G., J., Carney, S., Kumar, M., van Werkhoven, K., Shen, C., Wood, A., W., Clark, M. P. (2026). *Comparing multi-model mosaic and multi-model combination methods to simulate streamflow across the contiguous USA*. Hydrology and Earth System Sciences.
