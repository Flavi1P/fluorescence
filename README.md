# Influence of Phytoplankton Communities on In-Situ Fluorescence 

## Project Objective
The purpose of this project is to investigate the influence of phytoplankton community composition on the in situ fluorescence signal. The phytoplankton community composition is evaluate through HPLC data.


### Data
* Argo Floats Output - all data needed for argo fluo vs hplc analysis
* Boussole - all data needed for the boussole analysis
* Longhurst - Longhurst bioregion data
* Raw HPLC Argo - raw HPLC from BGC argo campaign
* absorption - data for absorption analysis

### Scripts
* abs_data - Exploration of aphy dataset (deprecated)
* aphy_boussole_shaping - Shaping of the raw boussole dataset
* argo_table - format argo ncdf into tabular data
* bouss_acp - create ACPs from Boussole HPLC/Fluo data
* bouss_merge - merge abs and pig data at Boussole
* bplr_pig - format pigments from bplr region
* cluster_abs - create clusters from absorption data (deprecated)
* data_viz - some exploration (deprecated)
* glo_aphy_viz data exploration of aphy (deprecated)
* hplc_shaping - format hplc raw data
* match_dist - match hplc and argo unsupervised
* match_onebyone - match hplc and argo, float by float, best to use
* merge - matchup visualization
* new_merge - complete the bgc/hplc matchup
* plot_profiles - function to plot profiles
* read_hplc - function to read raw HPLC xsls data
* region_plot - compute figures of argo/hplc dataset
* soclim - exploration on soclim campagne data (deprecated)
* soclim_abs - exploration on soclim campagne data (deprecated)
* table_pigment - create a table of pigment and absorption
* ternary_diagramm - visualisation of the slope factor in a ternary diagram
* zeu_moma - function to compute Zeu

### Output

The different figures for the analysis.


## Contact
* If you have any question you can contact me on my github profile or at flavien.petit@imev-mer.fr  
