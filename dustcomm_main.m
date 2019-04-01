% ==============
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles
% ==============
% These are the codes used to produce DustCOMM (v1) data published in Adebiyi
% et al (GMD, 2019).
% ==============
% first clear memory pane.
  clear all;clc;
% ==============

% -----------
% Model simulations to use.
  modelnames = {'GISS','WRFChem','CESM','GEOSChem','CNRM','IMPACT'};
  
% -----------
% number of different realizations to calculate...
  nboot = 1500;
  
%  necessary parameters and directories are defined in
%  "define_global_variables.m"
  define_global_variables (true);
  global  clat clon nThreads_Max  ...
          use_desktop_computer use_parallel_pool_in_function ...
          out_dir_fit_annu dustcomm_dir_psd_annu dustcomm_dir_seas ...
          dustcomm_dir_psd dustcomm_dir_mee out_dir_fit_seas

%% ===============================
% Annually-averaged climatology
% -----------

% -----------
% Step 1 -- Estimate the constraints over each location.
% -----------
% Procedure:
%   1. Collect the uncorrected dust mass fraction from all the model
%   simulations
%   2. Bias correct these dust mass fractions for each location, using the globally-averaged PSD from Kok et al, 2017.
%   3. Regularize the model data, and scale the sizes to between 0.2 and 20 microns
%   3. Calculate the constrained PSD by fitting a generalized analytical function to the bias-corrected dust mass fraction over each location and height
%   4. Save the result for each location.

% -----------
  [~] = constrain_biasCorr_dust_fraction_clim_seas_3d (modelnames,0,'annu',true);
  disp('STEP 1 COMPLETE.....');
% -----------
% You dont have to collect the output because they are saved for each
% location.
% -----------
% Step 2 -- Use the annually-averaged constrained PSD estimated from the
% different model to calculate the most likely constrained PSD over each
% location.
% -----------
% Procedure:
% 1. Collect the constrained PSD data for all the models
% 2. Randomly scale the mean distribution based on the realizations of the
% globally-averaged PSD from Kok et al, 2017
% 3. Use that with the constrained PSD dat randonly selecting from
% one of the models
% -----------
  calc_annu_constrained_PSD_bootstrap (modelnames,nboot,true);
  disp('STEP 2 COMPLETE.....');
% -----------
% Step 3 --  Now calculate the dust extinction and dust loading using the constrained PSD
% -----------
  [~] = calc_annu_constrained_MEE_Load_bootstrap (modelnames,nboot,true);
  disp('STEP 3 COMPLETE.....');
% -----------
% Step 4 --  Consolidate all data for all locations.
% Only do this if you have calculated for every location on the globe
% -----------
  consolidate_DustCOMM_annu (nboot,true);
  disp('STEP 4 COMPLETE.....');
% -----------
% Clean up -- Delete location Files
% -----------
% % where each model annully-averaged PSD are stored
%   system(char(strcat('/bin/rm -rf',{' '},out_dir_fit_annu))); 

% % where each location annually-averaged MEE are stored
%   system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_mee_annu))); 
% -----------
  
%% ===============================
% Seasonally-averaged climatology
% ---------------------------------

% -----------
% Step 1 -- Estimate the constraints over each location.
% -----------
% % get all the locations
  nlat = max(size(clat)); %
  nlon = max(size(clon)); %
% -----------
  % Define if you want to use many worrkers or not
  if (use_parallel_pool_in_function)
    if (nThreads_Max ~= 0)
      if (use_desktop_computer)
      % If a desktop computer is in use, this is set to > 1 to trigger
      % default later in the code.
        nThreads = 2; % number of workers to use.
      else
        nThreads = ceil(nlat/2); % number of workers to use.
      end
    else
      nThreads = nThreads_Max;
    end %nThreads_Max
  else
    nThreads = -1; % use only one worker
  end  %use_parallel_pool_in_function


  % First calculate
%   parfor ilat =1:nlat
  for ilat =1:nlat
    
    if ((nThreads > 1) && ~use_desktop_computer)
      pool = parpool(nThreads,'IdleTimeout', Inf); %Start the parallel pool using the 'local' profile.
%           The pool is not defined here for desktop computer. When Parfor
%           command is called later, the code defaults to the standard
%           number of workers allowed on the computer.
    end
    
    parfor ilon =1:nlon
%     for ilon =1:nlon
      [ilat nlat ilon nlon]
      [exitFlag] = calc_seas_constrained_PSD_MEE_Load_bootstrap (ilon,ilat,modelnames,nboot,false);
    end %ilon
    delete(gcp('nocreate')); % First close any existing pool
    
  end % ilat
  
  delete(gcp('nocreate')); % First close any existing pool
  disp('STEP 1 COMPLETE.....');
 
 % -----------
% Step 4 --  Consolidate all data for all locations.
% Only do this if you have calculated for every location on the globe
% -----------
  consolidate_DustCOMM_seas (true);
  disp('STEP 2 COMPLETE.....');
% ----------- 
  
% ===============================

% % -----------
% % Clean up Delete folders
% % -----------

% % where each model seasonally-averaged PSD are stored
%   system(char(strcat('/bin/rm -rf',{' '},out_dir_fit_seas))); 

% % where bootstrapped annully-averaged PSD are stored
%   system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_psd_annu))); 

% % where the seasonally-averaged MEE are stored
%   system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_mee_seas))); 

% % where the combined seasonally-averaged PSD-MEE-Load are stored
%   system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_seas))); 

% ------------
%   listing = dir(dustcomm_dir_psd); check_f = numel(listing) == 2 && isequal({listing.name},{'.','..'});
%   if check_f; system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_psd))); end
%   
%   listing = dir(dustcomm_dir_mee); check_f = numel(listing) == 2 && isequal({listing.name},{'.','..'});
%   if check_f; system(char(strcat('/bin/rm -rf',{' '},dustcomm_dir_mee))); end


  disp('All Done...')
