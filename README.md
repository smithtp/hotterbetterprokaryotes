# Hotter is better prokaryote meta-analysis

###################################

## Table of Contents
1. Overview
2. Repository Contents
3. Software Requirements
    1. Python
    2. R
4. Instructions for Use
    1. Demo
    2. Main Analysis
    3. Further Analysis and Figures
5. Citation

## Overview

Data and code for prokaryote thermal sensitivity project, bioRxiv pre-print doi: https://doi.org/10.1101/524264 

## Repository Contents

Repository contains the full raw prokaryote TPC dataset, alongside the thermal fitting code (python) and analysis scripts (R) required to reproduce the main results and figures for our bioRxiv pre-print.

  * /Data/ - contains the raw datasets.
  * /Code/ - contains the python wrapper scripts to run the fitting code.
  * /Code/Thermal_Models/ - the underlying thermal fitting code. A thorough README.md file describes the code and functions in detail.
  * /Code/R_scripts/ - analysis code written in R to produce results figures.
  * /Results/ - summaries of fitting results, TPC fit figures and results figures produced for manuscript.
  * /Demo/ - short demo of the fitting code


## Software Requirements

Analysis code has been tested in Windows 10 and Ubuntu 15.04.  

### Python

The thermal model fitting code was written using python version 3.5.1 and has also been tested in version 3.4.3. The thermal fitting code has a number of dependencies:

  * lmfit (0.9.2)
  * Numpy (1.11.0)
  * Pandas (0.17.0)
  * Matplotlib (1.5.1)
  * Seaborn (0.4.0) 
  * Scipy (0.16.0)    
  * progress (1.2)

Packages that the code was written with in parentheses; code has not been tested with other versions of these packages.

### R

The bulk of the code to analyse the fitting results and produce the manuscript figures was written in R version 3.2.2. The analysis code also has a number of package requirements:

  * Segmented (0.5-1.4)
  * ggplot2 (2.2.1)
  * gridExtra (2.2.1)
  * gtable (0.2.0)

Code has again not been tested with other versions of these packages.

## Instructions for Use

Clone the git repository

    > git clone https://github.com/smithtp/hotterbetterprokaryotes.git

or download zip

    https://github.com/smithtp/hotterbetterprokaryotes/archive/master.zip

The repository requires no installation or compilation of code once downloaded. Ensure python, R and required packages are installed (see software requirements), and fitting and analysis scripts can then be run. Installation time - less than a minute to download and extra package. Installation of python and R dependancies will increase this time.

### Demo

A short demo of the Schoolfield TPC fitting code and result outputs - try this first. Navigate to /Demo/ - Folder contains own README.md instructions.

### Main Analysis

First navigate to /Code/. To run through main fitting code and analysis:

    > python 01_fit_from_database.py
    > python 02_clean_summary.py

These read the database and attempt to fit the two-factor Schoolfield model (or other models as specified) to each TPC ID in the database. The models are plotted and parameter estimates collected in a summary csv. The summary is then cleaned to provide only a subset of the results fields and only the growth rate data for further analysis.

    R_Scripts/03_segmented_analysis.R
    R_Scripts/04_reattach_temp_prefs.R

These R scripts perform the segmented analysis that we used to determine the temperature break-points in the data (see manuscript methods). These break-points are used to assign strains as mesophiles or thermophiles and this information added to the summary data.

Note that these R scripts set a working directory (location of the script within the file system) and then use relative paths to read from /Data/ and write to /Results/. Users should edit the "setwd()" line to the location of the /Code/ directory within their own file architecture. Relative paths are written for linux directory structure, so Windows users would need to alter these within the scripts (forward-slashes to backslashes).

    > python 05_calc_group_means_no_agg.py
    > python 06_fit_to_summary_no_agg.py

Here we calculate the group means (short term responses, E<sub>S</sub> in the manuscript) before fitting boltzmann-arrhenius to the TPC peaks to calculate the evolutionary thermal sensitivity (long term response,E<sub>G</sub> in the code, now E<sub>L</sub> in the manuscript!).

    R_Scripts/07_Final_analysis_temppref.R
    R_Scripts/08_temperature_figures.R

Finally, we compare E<sub>S</sub> and E<sub>L</sub> for our taxonomic and temperature groupings to ask whether "hotter is better" (see manuscript methods).

### Further Analysis and Figures

First we need to remove pseudoreplicates from the autotroph dataset:

    R_Scripts/Autotroph_data.R

Subsequently, we compare E of growth vs E of fluxes in our dataset:

    > python flux_clean_summary.py
    R_Scripts/flux_clean_data.R
    > python flux_group_means.py
    R_Scripts/Flux_vs_Growth.R

These perform a new clean of the original summary and calculate the group means, E<sub>L</sub>, before producing the associated results figure.

pathogens.R performs a subsequent categorisation of the microbial strains by whether they are pathogenic and if so, what their pathogen carriers are - see manuscript methods for detais.

Figure production:

    Figure 2 produced by combining plot 1_Temp_Overlap_Plot_2.png (from Code/07_Final_analysis_temppref.R) with archaea_TempPref_weighted_axis_flip.png and bacteria_TempPref_weighted_axis_flip (from Code/08_temperature_figures.R)
    Figure 3 produced by combining ES_plot.png and ES_plot_funct.png from Code/ES_figure.R
    Figure 4 produced directly in Flux_vs_Growth.R
    Figure 5 produced by combining the output figures and additional legends from Code/Q10_figure.R


## Citation

If using this thermal fitting code, please cite our bioRxiv pre-print (https://doi.org/10.1101/524264)


