# Dust Constraints from joint Observational-Modelling-experiMental analysis – DustCOMM Version 1

If you have any question, please email the corresponding author, Adeyemi Adebiyi. at aadebiyi@ucla.edu.

==========================

This repository contains the codes used to generate the version-1 of DustCOMM dataset (Adebiyi et al. Geoscientific Model Development), which are annual and seasonal climatologies of constrained dust aerosol properties, including spatially-varying dust size distribution, mass extinction efficiency, and atmospheric dust loading.
This code uses modeling constraints from six global models, namely GISS, WRFChem, CESM, GEOSChem, CNRM (ARPEGE), IMPACT.

The consolidated DustCOMM (v1) output data can be downloaded on Zenodo using this link http://zenodo.org/record/2620475

In addition to this ReadMe file, three other files (dustcomm_main.m, define_global_variables.m, and run_dustcomm) and three folders (model_data, functions, use_script) are provided in this repository.
Below, we give a brief description of each file/folder.

--------
The three (3) files
--------
1.  dustcomm_main.m           -- This is the main MATLAB file that produce DustCOMM's annually-averaged and seasonally-averaged estimates.
2.  define_global_variables.m -- Global parameters and variables are defined here. You can also define the  directories for both input and output datasets.
3.  run_dustcomm              -- simple Unix code to run "dustcomm_main.m" in the background and parse the standard output & error into a file "matlabjob.txt". This is perhaps better not suitable for a desktop computer.

--------
The three (3) Folders
--------
1.  functions   -- This folder contains main functions used in DustCOMM
2.  use_script  -- This folder contains auxiliary functions used in DustCOMM
3.  input_data  -- This folder contains the input data into dustCOMM.

Below, we give a brief description of the files in each folder.

--------
The "functions" folder
--------
1.  constrain_biasCorr_dust_fraction_clim_seas_3d.m   -- Function that corrects and constrain the model dust size distribution for each location
2.  calc_annu_constrained_PSD_bootstrap.m             -- Function that estimate the constrained annually-averaged dust size distribution
3.  collect_annu_constrained_PSD_bootstrap.m          -- Function that is used to collect a defined number of realizations of the constrained annually-averaged dust size distribution
4.  calc_annu_constrained_MEE_Load_bootstrap.m        -- Function that estimates constrained dust mass extinction efficiency and atmospheric loading using the constrained dust size distribution
5.  consolidate_DustCOMM_annu.m                       -- Function that consolidate the annually-averaged size distribution, dust mass extinction efficiency and dust loading...in a few, easy to read files.
6.  calc_seas_constrained_PSD_MEE_Load_bootstrap.m    -- Function that estimate the constrained seasonally-averaged dust size distribution, dust mass extinction efficiency and dust loading
7.  consolidate_DustCOMM_seas.m                       -- Function that consolidate the seasonally-aveaged size distribution, dust mass extinction efficiency and dust loading...in a few, easy to read files.
8.  bias_correction_all_3d_v1.m                       -- FUnctioon that bias-correct the six modelled dust size distribution
9.  scale_bicor_dustfraction_model_3d.m               -- Envelop function to call scale_model_no_3d.m
10. scale_model_no_3d.m                               -- Function that scale bias-corrected size distribution
8.  generalized_PSD_nonlinear_function.m              -- Generalized size distribution function used to constrain the model dust mass fraction.
9.  get_model_dust_concentration_data.m               -- Function that collects the model 3-D dust concentration
10. get_model_dust_Vweight_data.m                     -- Function that collects the model vertical dust weight for each location. This is calculated for eah model

--------
The "use_script" folder
--------
1.  chi_square_function_subbin.m -- Function that defines the Chi-square equation which is minimum to determine the constrained dust size distribution
2.  fminsearchbnd.m             -- Modified MATLAB function to find the minimum of unconstrained multivariable function using derivative-free method
3.  Gaussian.m                  -- Simple function to generate a random value of a Gaussian distribution
4.  Gaussian2D.m                -- Simple function to generate 2-D random values for a Gaussian distribution for each location.
5.  PSD_Plaw.m                  -- Function calculates power law distribution in log-space, using a minimum of two pointwise values.
6.  PSD_V_integral.m            -- Simple function to calculate the intergral between two diameter for dV/dD, Using the dV/dlnD
7.  calculate_logmean_error.m   -- Function used to calculate the log mean error. For a constrained PSD with error more than a threshold, an secondary procedure is taken to minimum the error. See generalized_PSD_nonlinear_function.m
6.  delete_file_ifexist.m       -- Simple function to delete an external file if it exist.
7.  write_netcdf.m              -- Simple function to write it in netcdf output file.

--------
The "input_data" folder (data can be downloaded at http://zenodo.org/record/2620547)
--------
1.  model_dust_concentration          -- Simulated 3-D dust concentration (ug/m3) from 6 model simulations. This folder contains both annually-averaged and seasonally-averaged model simulations. 
2.  model_dust_Vweight                -- 3-D model dust vertical weight used in DustCOME. This is calculated from the modeled dust concentration. This folder contains both annually-averaged and seasonally-averaged model simulations.
3.  DAOD_clim_reanalysis_all.nc       -- Constrained annually-averaged dust aerosol optical depth. See Adebiyi et al. Geoscientific Model Development for details
4.  DAOD_seas_reanalysis_all.nc       -- Constrained seasonally-averaged dust aerosol optical depth. See Adebiyi et al. Geoscientific Model Development for details
5.  DustCOMM_lat.nc                   -- Latitude array used in DustCOMM. Dimension = 96.
6.  DustCOMM_lon.nc                   -- Longitude array used in DustCOMM. Dimension = 144.
7.  DustCOMM_lev.nc                   -- Altitude (hPa) array used in DustCOMM. Dimension = 48.
8.  Globally_averaged_dust_PSD.mat    -- Constrained Globally-averaged dust size distribution from Kok et al, Nature Geoscience, 2017.
9.  Size_resolved_ext_eff.mat         -- Constrained Globally-averaged dust extinction efficiency from Kok et al, Nature Geoscience, 2017.

==========================

