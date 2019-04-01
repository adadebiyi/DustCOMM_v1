function define_global_variables (define_param)
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This code defines the global variables used in DustCOMM

% Input variable:
% define_param - A boolean to define global parameters .

  global out_dir_vwgt ...
        out_dir_fit_annu out_dir_fit_seas out_dir_conc ...
        dustcomm_dir_psd_annu dustcomm_dir_psd_seas ...
        mean_rho_d error_rho_d ...
        fitPSD_filename fit_error_threshold ...
        fitPSD_filename_annu fitPSD_filename_seas ...
        gQe_OBS gPSD_OBS gPSD_OBS_median D_min D_max D_OBS ...
        clat clon clev max_level_in_hPa ...
        use_parallel_pool_in_function nThreads_Max use_desktop_computer ...
        DOAD_data_annu DOAD_data_seas ...
        dustcomm_dir_mee_annu mee_filename_annu ...
        dustcomm_dir_mee_seas mee_filename_seas dustcomm_dir_seas...
        dustcomm_dir_all dustcomm_dir_psd dustcomm_dir_mee seas_filename

% -------------
%   Version of DustCOMM to use
  expname = 'v1';
% -------------
%  working directory. Where the "ductcomm_main.m" document is located
  my_pwd        = pwd;
  dustcomm_dir  = char(strcat('~/data1/DustCOMM_',expname,'_data/'));

% Define the path for scripts, functions and data...
  script_dir = char(strcat(my_pwd,'/use_script/'));  if (exist(script_dir, 'dir') ); addpath(script_dir); end %external scripts used in the analysis
  func_dir = char(strcat(my_pwd,'/functions/'));  if (exist(func_dir, 'dir') ); addpath(func_dir); end %written functions used in the analysis
  in_data_dir = char(strcat(my_pwd,'/input_data/'));  if (exist(in_data_dir, 'dir') ); addpath(in_data_dir); end % model data used in this analysis

% % COMMENT OUT...OR DELETE
% % ===========
  system(char(strcat('/bin/mkdir -p',{' '},'./input_data/'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/bigdata/model_dust_concentration',{' '},in_data_dir))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/bigdata/model_dust_Vweight',{' '},in_data_dir))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/bigdata/model_dust_Vweight',{' '},in_data_dir))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/postprocess/CESM_grid_DAOD_clim_reanalysis_all.nc',{' '},in_data_dir,'DAOD_clim_reanalysis_all.nc'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/postprocess/CESM_grid_DAOD_seas_reanalysis_all.nc',{' '},in_data_dir,'DAOD_seas_reanalysis_all.nc'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/postprocess/lat.nc',{' '},in_data_dir,'lat.nc'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/postprocess/lon.nc',{' '},in_data_dir,'lon.nc'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/postprocess/lev.nc',{' '},in_data_dir,'lev.nc'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/jfk_data/Globally_averaged_dust_PSD.mat',{' '},in_data_dir,'Globally_averaged_dust_PSD.mat'))); % 
  system(char(strcat('/bin/ln -sf ',{' '},'~/data/jfk_data/Size_resolved_ext_eff.mat',{' '},in_data_dir,'Size_resolved_ext_eff.mat'))); % 
% % ===========

  out_dir_conc = char(strcat(in_data_dir,'/model_dust_concentration/'));
  out_dir_vwgt = char(strcat(in_data_dir,'/model_dust_Vweight/'));
  
%   Directory to store constrained PSD
  dustcomm_dir_psd      = char(strcat(dustcomm_dir,'/constrained_PSD/'));
  dustcomm_dir_mee      = char(strcat(dustcomm_dir,'/constrained_MEE_Load/'));
  dustcomm_dir_all      = char(strcat(dustcomm_dir,'/DustCOMM_',expname,'_final_output/'));
  dustcomm_dir_seas = char(strcat(dustcomm_dir,'/seas_all_constrained/'));
  
  out_dir_fit_annu      = char(strcat(dustcomm_dir_psd,'/constrained_model_dust_PSD_annu/')); 
  dustcomm_dir_psd_annu = char(strcat(dustcomm_dir_psd,'/annu/'));
  dustcomm_dir_mee_annu = char(strcat(dustcomm_dir_mee,'/annu/'));

  out_dir_fit_seas      = char(strcat(dustcomm_dir_psd,'/constrained_model_dust_PSD_seas/'));
  dustcomm_dir_psd_seas = char(strcat(dustcomm_dir_psd,'/seas/'));
%   dustcomm_dir_mee_seas = char(strcat(dustcomm_dir_mee,'/seas/'));
  
% the filename for fit
  fitPSD_filename = 'constrained_dust_PSD'; % name for to store the constrained.
  fitPSD_filename_annu = 'constrained_dust_PSD_annu'; % name for to store the constrained.
  fitPSD_filename_seas = 'constrained_dust_PSD_seas'; % name for to store the constrained.
  
  mee_filename_annu = 'constrained_dust_MEE_Load_annu'; % name for to store the constrained.
  mee_filename_seas = 'constrained_dust_MEE_Load_seas'; % name for to store the constrained.
  seas_filename = 'seas_constrained';
  
  %Error threshold for the sub-bin calculation
  fit_error_threshold = 1; % % maximum fitting logmean error allowed for the constrained PSD

% the diameter limits
  D_min = 0.2; % the minimum diameter for DustCOMM
  D_max = 20; % the minimum diameter for DustCOMM

% define if you'd like to use the parrallel pool
  use_parallel_pool_in_function = true; % 0 is false and 1 is true

  % Specify if the desktop computer is in use, or a more powerful machine.  
%   use_desktop_computer = false;
  use_desktop_computer = ~isempty(strfind(pwd,'home'));
  
%Specify the  maximum number of workers allowed for the matlab code.
% The maximum usable cores by the code is 48.
% If you want to use default...set nThreads_Max to 0
% Note that if use_parallel_pool_in_function only one worker is used
  nThreads_Max = 48; 
  
% maximum pressure level needed in hPa
% DustCOMM can constrain up to 10hPa
  max_level_in_hPa = 100;
  
% Define the value of dust density to use in calculating the extinction
  mean_rho_d = 2500;  %density of dust particles, in kg/m3
  error_rho_d = 200;  %error in dust density

% specify the path to constrained dust AOD 
% the mean and error are stored as "DAOD_m" and "DAOD_sd"
  daod_filename = 'DAOD_clim_reanalysis_all.nc';
  DOAD_data_annu = strcat(in_data_dir,daod_filename);
  daod_filename = 'DAOD_seas_reanalysis_all.nc';
  DOAD_data_seas = strcat(in_data_dir,daod_filename);

% -------------
  if (define_param)
  % =================
    % Calculate the median of the PSD from Kok et al, 2017
    PSDfile = strcat(in_data_dir,'Globally_averaged_dust_PSD.mat');
    load (PSDfile, 'PSD_lnV_load_norm','D');   %
    ibin = find(D>=D_min & D<=D_max);
    gPSD_OBS = PSD_lnV_load_norm(:,ibin);
    D_OBS = D(1,ibin);
    gPSD_OBS_median = zeros(1,size(D_OBS,2));
    for i = 1:size(D_OBS,2)
      A=sort(PSD_lnV_load_norm(:,i));
      gPSD_OBS_median(i) = A(round(0.5*size(gPSD_OBS,1))); %mean bootstrap PSD, since it is gaussian for each bin
    end
    
    % =================
    % Collect the extinction at 550 nm
    Qefile = strcat(in_data_dir,'Size_resolved_ext_eff.mat');
    load (Qefile, 'Qe');   %
    gQe_OBS = Qe(:,ibin);

    clear A ibin D PSD_lnV_load_norm Qe
    
% =================
  % Longitide and latitude used in DustCOMM
    clat = ncread(strcat(in_data_dir,'lat.nc'),'lat');
    clon = ncread(strcat(in_data_dir,'lon.nc'),'lon');
    clev = ncread(strcat(in_data_dir,'lev.nc'),'lev');

  end %define_true

 

%   DOne...
end
