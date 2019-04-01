function [constrained_dust_PSD,param_min_latlon,range_warning_latlon,biCorr_dust_fraction_norm] = constrain_biasCorr_dust_fraction_clim_seas_3d(model_name,reg_lonlat,annu_or_seas,save_data_only_netcdf)
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This functions  fit the generalized PSD function to the 3D global
% bias-corrected dust mass fraction.

% ==========
% Steps:
%   1. the uncorrected global dust mass fraction is collected
%   2. bias correct the dust mass fraction for each location, based on the globally-averaged PSD from Kok et al, 2017.
%   2a. regularize the data, and scale to between 0.2 and 20 microns
%   3. A generalized analytical function is fit into each mass fraction over
%   each location and height using the generalized function
%   4. the result is saved for each model.

% ==========
% Variables:
% constrained_dust_PSD . -- contains new fitted bias-corrected dust mass distribution
% param_min_latlon = [Ds_min,sigma_s_min,a_min,b_min,Cv_min]; -- fitting parameters for each location
% range_warning_latlon = [Ds_range_warning,sigma_s_range_warning,a_range_warning,b_range_warning]; flag showing if the parameter is at the end of the chi-square
% biCorr_dust_fraction_norm  -- the new normalized bias-corrected dust fraction
%
% % ==========
% % Example
% clear all;clc;
% % model_name = 'GISS'; %
% model_name = {'GISS', 'WRFChem', 'CESM','GEOSChem','CNRM','IMPACT'};
% % reg_lonlat = 0;
% reg_lonlat = [0 2; 0,1];
% annu_or_seas='seas';
% save_data_only_netcdf = true;
% % ==========


% Begin here....
% =================

%  define some global variables
  define_global_variables (true);
  global  out_dir_fit_annu out_dir_fit_seas ...
          fitPSD_filename fit_error_threshold ...
          D_OBS clat clon clev ...
          nThreads_Max use_desktop_computer use_parallel_pool_in_function

%  use the appropriate one
  if (annu_or_seas == char('seas'))
     out_dir = out_dir_fit_seas;
  else
    out_dir = out_dir_fit_annu;
  end
  system(char(strcat('/bin/mkdir -p',{' '},out_dir))); % make the directory if not available

% general names of seasons
  seas_str = {'DJF', 'MAM', 'JJA','SON'}; % seasons
%  convert char, if not...
  annu_or_seas = char(annu_or_seas);
% copy for the different workers
  xD_OBS = D_OBS;
  xfit_error_threshold = fit_error_threshold;

% =================
 % information for the spatial range
  if (reg_lonlat == 0)
  % default
    indlat = 1:max(size(clat));
    indlon = 1:max(size(clon));
  else
    % make it a box, even if it is one lon and lat
    if (size(reg_lonlat,2) == 1 || size(reg_lonlat,1) == 1)
      [d, indlat] = min(abs( clat-reg_lonlat(2)));
      [d, indlon] = min(abs( clon-reg_lonlat(1)));
    else
%       check if regional or individual location
      if (size(reg_lonlat,2) == 2 && size(reg_lonlat,1) == 2)
        %     reg_lonlat = [xlonW xlonE; xlatS xlatN];
        xlatN = reg_lonlat(2,2);
        xlatS = reg_lonlat(2,1);
        xlonE = reg_lonlat(1,2);
        xlonW = reg_lonlat(1,1);
        indlat = find(clat >= xlatS & clat <= xlatN);
        indlon = find(clon >= xlonW & clon <= xlonE);
      else

    % ; get the indices for the logitude and latitude  and height
        indlat = NaN(1,max(size(reg_lonlat)));
        indlon = indlat;
        for i=1:max(size(reg_lonlat))
          [d, indlat(i)] = min(abs( clat-reg_lonlat(i,1)));
          [d, indlon(i)] = min(abs( clon-reg_lonlat(i,2)));
        end
      end
    end
  end

  xlon = clon(indlon);
  xlat = clat(indlat);

  clear clon clat;

% =================
%dimension
  nlat = max(size(indlat));
  nlon = max(size(indlon));
  nlev = max(size(clev));
  nD = size(D_OBS,2);
  nD_model = 8; % this is the maximum bin for models so far
  nseas = max(size(seas_str));
  nmodel = max(size(model_name));

% =================
%   number of matlab workers to use.
%  Check the work station the code will be run. If it is super computer
%  usually having "home", then use higher number.

  if (use_parallel_pool_in_function)
    if (nThreads_Max ~= 0)
%       if (isempty(strfind(pwd,'home')))
      if (use_desktop_computer)
        % If a desktop computer is in use, this is set to > 1 to trigger
        % default later in the code.
         nThreads = 2; % number of workers to use.
      else
         nThreads = nlev; % number of workers to use.
      end
    else
      nThreads = nThreads_Max;
    end %nThreads_Max
  else
    nThreads = -1; % use only one worker
  end  %use_parallel_pool_in_function

% =================
%  Check if any of the output files alreadt exist.
%  boolean to check if the file already exist
  exist_models = zeros(1,nmodel);
  if (annu_or_seas == char('annu'))
    for ix=1:nmodel
      for ilat=1:nlat
        for ilon=1:nlon
          xlon_str = sprintf('%0.2f',xlon(ilon));
          xlat_str = sprintf('%0.2f',xlat(ilat));
          conc_file = char(strcat(out_dir,model_name(ix),'_',fitPSD_filename,'_',xlon_str,'_',xlat_str,'_annu.nc'));
          if (exist(conc_file, 'file'))
            exist_models(ix) = exist_models(ix) + ix;
          else
            exist_models(ix) = exist_models(ix) + NaN;
          end
        end % ilon
      end % ilat
    end % ix
  elseif (annu_or_seas == char('seas'))
    for is=1:nseas
      for ix=1:nmodel
        for ilat=1:nlat
          for ilon=1:nlon
            xlon_str = sprintf('%0.2f',xlon(ilon));
            xlat_str = sprintf('%0.2f',xlat(ilat));
            conc_file = char(strcat(out_dir,model_name(ix),'_',fitPSD_filename,'_',xlon_str,'_',xlat_str,'_',seas_str(is),'.nc'));
            if (exist(conc_file, 'file'))
              exist_models(ix) = exist_models(ix) + ix;
            else
              exist_models(ix) = exist_models(ix) + NaN;
            end
          end % ilon
        end % ilat
      end % ix
    end
  end % annu_or_seas

  do_calc = any(isnan(exist_models)); %to know if to do calculation or not.
%    =================
%
  if (annu_or_seas == char('annu'))
    disp(' ')
    disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Running....For annually-averaged calculation...')
    disp(' ')

%     where to store the data...
    nparam = 6; % number of parameters

    if (~save_data_only_netcdf)
      disp('constrain_biasCorr_dust_fraction_clim_seas_3d: save_data_only_netcdf is false --> This code may be very slow.')
      disp('Consider making save_data_only_netcdf to only do the calcution and store the data...')
      constrained_dust_PSD = NaN(nD,nlon,nlat,nlev,nmodel);
      param_min_latlon = NaN(nparam,nlon,nlat,nlev,nmodel);
      range_warning_latlon = NaN(nparam,nlon,nlat,nlev,nmodel);
      biCorr_dust_fraction_norm = NaN(nD_model,nlon,nlat,nlev,nmodel);
    end %save_data_only_netcdf

  % =================
    if (do_calc)
  %  Load the uncorrected dust mass fraction
      disp(' ');
      disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Collecting the model dust loading...');
%       evalc('[model_column_loading_orig,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_all_3d(annu_or_seas,model_name,false)');
%       evalc('[model_column_loading_orig,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_data(annu_or_seas,model_name)');
      [model_column_loading_orig,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_data(annu_or_seas,model_name);
      clear llon llat llev;
    end % do_calc

%     loop over the models
    for m=1:nmodel
      disp(['constrain_biasCorr_dust_fraction_clim_seas_3d: Continue only with...',model_name(m)]);
      disp('');

      if (do_calc && isnan(exist_models(m)))
        % convert to normalized fraction
        xdum = squeeze(model_column_loading_orig(m,1:Nbin(m),:,:,:));
        sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),Nbin(m),1,1,1,1));
        column_loading  = xdum./sum_xdum;

        clear xdum sum_xdum;

        % get the diameter limits
        Dlower = D_lower(m,1:Nbin(m));
        Dupper = D_upper(m,1:Nbin(m));

        % apply the bias correction here...
        [biCorr_dust_fraction] = bias_correction_all_3d_v1(column_loading,0,Dlower,Dupper,model_name(m),annu_or_seas);

        %  scale the bias corrected version and extend to 20 micron
        [biCorr_dust_fraction_scale,Dlower,Dupper,model_Nbin] = scale_bicor_dustfraction_model_3d(biCorr_dust_fraction,Dlower,Dupper,model_name(m),annu_or_seas);

        % apply the bias correction here again...
        [biCorr_dust_fraction_scale] = bias_correction_all_3d_v1(biCorr_dust_fraction_scale,0,Dlower,Dupper,model_name(m),annu_or_seas);

        ind_model_Nbin = 1:model_Nbin;
        %normalize
        sum_xdum = squeeze(repelem(sum(biCorr_dust_fraction_scale(ind_model_Nbin,:,:,:),1,'omitnan'),model_Nbin,1,1,1));
        biCorr_dust_fraction_scale_norm  = biCorr_dust_fraction_scale./sum_xdum;

        clear biCorr_dust_fraction biCorr_dust_fraction_scale column_loading sum_xdum

        if ((nThreads > 1) && ~use_desktop_computer)
          pool = parpool(nThreads,'IdleTimeout', Inf); %Start the parallel pool using the 'local' profile.
%           The pool is not defined here for desktop computer. When Parfor
%           command is called later, the code defaults to the standard
%           number of workers allowed on the computer.
        end

      end % do_calc

      % now loop over the locations
      for ilat=1:nlat
        for ilon=1:nlon
%           Check if the location data is already saved
          xlon_str = sprintf('%0.2f',xlon(ilon));
          xlat_str = sprintf('%0.2f',xlat(ilat));
          conc_file = char(strcat(out_dir,model_name(m),'_',fitPSD_filename,'_',xlon_str,'_',xlat_str,'_annu.nc'));

        %  Check if the output file alreadt exist.
          if (~exist(conc_file, 'file'))
%           data to work with - pick them one at a time.
            loc_biCorr  = squeeze(biCorr_dust_fraction_scale_norm(:,indlon(ilon),indlat(ilat),:));

%           Check if atleast one level has all the bias corrected values.
%           Or if all is missing
            if (any(all(loc_biCorr)))

% Check if you only want to save data or want to do all the calculation together
              if (~save_data_only_netcdf)
% loop through the levels
                parfor ilev=1:nlev
  %                 for ilev=1:nlev
                    % disp([m nmodel ilon nlon ilat nlat ilev nlev])
                  if all((~isnan(loc_biCorr(:,ilev))))
                    [constrained_dust_PSD(:,ilon,ilat,ilev,m),param_min_latlon(:,ilon,ilat,ilev,m),range_warning_latlon(:,ilon,ilat,ilev,m),biCorr_dust_fraction_norm(ind_model_Nbin,ilon,ilat,ilev,m)] =  generalized_PSD_nonlinear_function(loc_biCorr(:,ilev), Dlower, Dupper, xD_OBS,0,0,xfit_error_threshold,false);
                  end
                end % ilev

                clear param_min range_warning chi_square C1 C2 C3 C4 C5;
  % now save the file here...
  %STore to a file
                x = delete_file_ifexist(conc_file);
                f = write_netcdf (conc_file, 2, {'nD','nlev'}, char('constrained_dust_PSD'), squeeze(constrained_dust_PSD(:,ilon,ilat,:,m)));
                f = write_netcdf (conc_file, 2, {'nD_model','nlev'}, char('biCorr_dust_fraction_norm'), squeeze(biCorr_dust_fraction_norm(ind_model_Nbin,ilon,ilat,:,m)));
                f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('param_min_latlon'), squeeze(param_min_latlon(:,ilon,ilat,:,m)));
                f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('range_warning_latlon'), squeeze(range_warning_latlon(:,ilon,ilat,:,m)));

              else

% if only data is to be saved...
% Redefine where to put them
                constrained_dust_PSD = NaN(nD,nlev);
                param_min_latlon = NaN(nparam,nlev);
                range_warning_latlon = NaN(nparam,nlev);
                biCorr_dust_fraction_norm = NaN(model_Nbin,nlev);

% disp([m nmodel ilon nlon ilat nlat ])
% now loop over the levels
                parfor ilev=1:nlev
  %                 for ilev=1:nlev
                    % disp([m nmodel ilon nlon ilat nlat ilev nlev])
                  if all((~isnan(loc_biCorr(:,ilev))))
                    [constrained_dust_PSD(:,ilev),param_min_latlon(:,ilev),range_warning_latlon(:,ilev),biCorr_dust_fraction_norm(ind_model_Nbin,ilev)] =  generalized_PSD_nonlinear_function(loc_biCorr(:,ilev), Dlower, Dupper, xD_OBS,0,0,xfit_error_threshold,false);
                  end
                end % ilev
                clear param_min range_warning dust_frac_norm;

  % now save the file here...
  %STore to a file
                x = delete_file_ifexist(conc_file);
                f = write_netcdf (conc_file, 2, {'nD','nlev'}, char('constrained_dust_PSD'), constrained_dust_PSD);
                f = write_netcdf (conc_file, 2, {'nD_model','nlev'}, char('biCorr_dust_fraction_norm'), biCorr_dust_fraction_norm(ind_model_Nbin,:));
                f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('param_min_latlon'), param_min_latlon);
                f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('range_warning_latlon'), range_warning_latlon);

              end %save_data_only_netcdf

% Add other important information
              f = write_netcdf (conc_file, 1, {'nD_model'}, char('Dupper'), Dupper(ind_model_Nbin));
              f = write_netcdf (conc_file, 1, {'nD_model'}, char('Dlower'), Dlower(ind_model_Nbin));
              f = write_netcdf (conc_file, 1, {'nbins'}, char('D'), D_OBS);
              f = write_netcdf (conc_file, 1, {'nlev'}, char('lev'), clev(1:nlev));
              f = write_netcdf (conc_file, 1, {'no_of_model_bins'}, char('N'), model_Nbin);

            end %loc_biCorr
            clear loc_biCorr;

          else
            % disp([m nmodel ilon nlon ilat nlat])
            disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Already calculated and stored.....Proceed to copy.')
            %Copy from the stored file here.
            if (~save_data_only_netcdf)
              N = ncread(conc_file,'N');
              constrained_dust_PSD(:,ilon,ilat,:,m)  = ncread(conc_file,'constrained_dust_PSD');
              param_min_latlon(:,ilon,ilat,:,m)  = ncread(conc_file,'param_min_latlon');
              range_warning_latlon(:,ilon,ilat,:,m)  = ncread(conc_file,'range_warning_latlon');
              biCorr_dust_fraction_norm(1:N,ilon,ilat,:,m)  = ncread(conc_file,'biCorr_dust_fraction_norm');
            else
              constrained_dust_PSD = 0;
              param_min_latlon = 0;
              range_warning_latlon = 0;
              biCorr_dust_fraction_norm = 0;
            end%save_data_only_netcdf
          end

        end %ilat
      end %ilat

      % delete(gcp('nocreate')); % First close any existing pool

    end % m

    clear model_column_loading_orig D_lower D_upper Nbin no_model_data

    disp(' ');
  elseif (annu_or_seas == char('seas'))

    disp(' ')
    disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Running....For seasonally-averaged calculation...')
    disp(' ')

%     where to store the data...
    nparam = 6; % number of parameters

    if (~save_data_only_netcdf)
      constrained_dust_PSD = NaN(nD,nlon,nlat,nlev,nseas,nmodel);
      param_min_latlon = NaN(nparam,nlon,nlat,nlev,nseas,nmodel);
      range_warning_latlon = NaN(nparam,nlon,nlat,nlev,nseas,nmodel);
      biCorr_dust_fraction_norm = NaN(nD_model,nlon,nlat,nlev,nseas,nmodel);
    end%save_data_only_netcdf

  % =================
    if (do_calc)
  %  Load the uncorrected dust mass fraction
      disp(' ');
      disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Collecting the model dust loading...');

      % first collect the annually-averaged
%       evalc('[model_column_loading_clim,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_all_3d(''annu'',model_name,false)');
      evalc('[model_column_loading_clim,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_data(''annu'',model_name)');

    % collect the seasonally-averaged
%       evalc('[model_column_loading_seas,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_all_3d(annu_or_seas,model_name,false)');
      evalc('[model_column_loading_seas,llon,llat,llev,D_lower,D_upper,Nbin,no_model_data] = get_model_dust_concentration_data(annu_or_seas,model_name)');

      clear llon llat llev;
    end % do_calc

    if ((nThreads > 1) && ~use_desktop_computer)
      pool = parpool(nThreads,'IdleTimeout', Inf); %Start the parallel pool using the 'local' profile.
    %           The pool is not defined here for desktop computer. When Parfor
    %           command is called later, the code defaults to the standard
    %           number of workers allowed on the computer.
    end

%     loop over the models
    for m=1:nmodel
      disp(['constrain_biasCorr_dust_fraction_clim_seas_3d: Continue only with...',model_name(m)]);
      disp('');

      if (do_calc && isnan(exist_models(m)))
        % convert to normalized fraction
        xdum = squeeze(model_column_loading_clim(m,1:Nbin(m),:,:,:));
        sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),Nbin(m),1,1,1,1));
        column_loading_ann  = xdum./sum_xdum;

        % convert to normalized fraction
        column_loading_seas = [];
        for ix=1:nseas
          xdum = squeeze(model_column_loading_seas(m,1:Nbin(m),:,:,:,ix));
          sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),Nbin(m),1,1,1,1));
          column_loading_seas(:,:,:,:,ix) = xdum./sum_xdum;
        end

        clear xdum sum_xdum ans;

        % get the diameter limits
        Dlower = D_lower(m,1:Nbin(m));
        Dupper = D_upper(m,1:Nbin(m));

        % apply the bias correction here...
        [biCorr_dust_fraction] = bias_correction_all_3d_v1(column_loading_seas,column_loading_ann,Dlower,Dupper,model_name(m),'seas');


%         do the same for annually averaged
        evalc('[biCorr_dust_fraction_ann] = bias_correction_all_3d_v1(column_loading_ann,0,Dlower,Dupper,model_name(m),''annu'')');
        evalc('[biCorr_dust_fraction_scale_ann,~] = scale_bicor_dustfraction_model_3d(biCorr_dust_fraction_ann,Dlower,Dupper,model_name(m),''annu'')');

        %  scale the bias corrected version and extend to 20 micron
        [biCorr_dust_fraction_scale,Dlower,Dupper,model_Nbin] = scale_bicor_dustfraction_model_3d(biCorr_dust_fraction,Dlower,Dupper,model_name(m),'seas');


        % apply the bias correction here again using the bias corrected and scaled annually averaged vales...
        [biCorr_dust_fraction_scale] = bias_correction_all_3d_v1(biCorr_dust_fraction_scale,biCorr_dust_fraction_scale_ann,Dlower,Dupper,model_name(m),'seas');

        ind_model_Nbin = 1:model_Nbin;
        %normalize
        sum_xdum = squeeze(repelem(sum(biCorr_dust_fraction_scale(ind_model_Nbin,:,:,:,:),1,'omitnan'),model_Nbin,1,1,1,1));
        biCorr_dust_fraction_scale_norm  = biCorr_dust_fraction_scale./sum_xdum;

        clear biCorr_dust_fraction biCorr_dust_fraction_scale column_loading_seas column_loading_ann sum_xdum
        clear biCorr_dust_fraction_scale_ann biCorr_dust_fraction_ann

      end % do_calc


      % now loop over the seasons and locations

      % begin loop.
      for iseas=1:nseas
        disp([seas_str(iseas)]);

        for ilat=1:nlat
          for ilon=1:nlon
  %           Check if the location data is already saved
            xlon_str = sprintf('%0.2f',xlon(ilon));
            xlat_str = sprintf('%0.2f',xlat(ilat));
            conc_file = char(strcat(out_dir,model_name(m),'_',fitPSD_filename,'_',xlon_str,'_',xlat_str,'_',seas_str(iseas),'.nc'));

          %  Check if the output file alreadt exist.
            if (~exist(conc_file, 'file'))
  %             data to work with
              loc_biCorr  = squeeze(biCorr_dust_fraction_scale_norm(:,indlon(ilon),indlat(ilat),:,iseas));

  %           Check if atleast one level has all the bias corrected values
              if (any(all(loc_biCorr)))

                if (~save_data_only_netcdf)

                  parfor ilev=1:nlev
  %                   for ilev=1:nlev
                      % disp([m nmodel ilon nlon ilat nlat ilev nlev iseas nseas])
                    if all((~isnan(loc_biCorr(:,ilev))))
                      [constrained_dust_PSD(:,ilon,ilat,ilev,iseas,m),param_min_latlon(:,ilon,ilat,ilev,iseas,m),range_warning_latlon(:,ilon,ilat,ilev,iseas,m),biCorr_dust_fraction_norm(ind_model_Nbin,ilon,ilat,ilev,iseas,m)] =  generalized_PSD_nonlinear_function(loc_biCorr(:,ilev), Dlower, Dupper, xD_OBS,0,0,xfit_error_threshold,false);
                    end
                  end % ilev

                  clear param_min range_warning chi_square C1 C2 C3 C4 C5;
    % now save the file here...
    %STore to a file
                  x = delete_file_ifexist(conc_file);
                  f = write_netcdf (conc_file, 2, {'nD','nlev'}, char('constrained_dust_PSD'), squeeze(constrained_dust_PSD(:,ilon,ilat,:,iseas,m)));
                  f = write_netcdf (conc_file, 2, {'nD_model','nlev'}, char('biCorr_dust_fraction_norm'), squeeze(biCorr_dust_fraction_norm(ind_model_Nbin,ilon,ilat,:,iseas,m)));
                  f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('param_min_latlon'), squeeze(param_min_latlon(:,ilon,ilat,:,iseas,m)));
                  f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('range_warning_latlon'), squeeze(range_warning_latlon(:,ilon,ilat,:,iseas,m)));
                else

                  constrained_dust_PSD = NaN(nD,nlev);
                  param_min_latlon = NaN(nparam,nlev);
                  range_warning_latlon = NaN(nparam,nlev);
                  biCorr_dust_fraction_norm = NaN(model_Nbin,nlev);

% disp([m nmodel ilon nlon ilat nlat iseas nseas])
                  parfor ilev=1:nlev
  %                   for ilev=1:nlev
                      % disp([m nmodel ilon nlon ilat nlat ilev nlev iseas nseas])
                    if all((~isnan(loc_biCorr(:,ilev))))
                      [constrained_dust_PSD(:,ilev),param_min_latlon(:,ilev),range_warning_latlon(:,ilev),biCorr_dust_fraction_norm(ind_model_Nbin,ilev)] =  generalized_PSD_nonlinear_function(loc_biCorr(:,ilev), Dlower, Dupper, xD_OBS,0,0,xfit_error_threshold,false);
                    end
                  end % ilev

                  clear param_min range_warning chi_square C1 C2 C3 C4 C5;
% now save the file here...
%STore to a file
                  x = delete_file_ifexist(conc_file);
                  f = write_netcdf (conc_file, 2, {'nD','nlev'}, char('constrained_dust_PSD'), constrained_dust_PSD);
                  f = write_netcdf (conc_file, 2, {'nD_model','nlev'}, char('biCorr_dust_fraction_norm'), biCorr_dust_fraction_norm(ind_model_Nbin,:));
                  f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('param_min_latlon'), param_min_latlon);
                  f = write_netcdf (conc_file, 2, {'nparam','nlev'}, char('range_warning_latlon'), range_warning_latlon);

                end %save_data_only_netcdf

                f = write_netcdf (conc_file, 1, {'nD_model'}, char('Dupper'), Dupper(ind_model_Nbin));
                f = write_netcdf (conc_file, 1, {'nD_model'}, char('Dlower'), Dlower(ind_model_Nbin));
                f = write_netcdf (conc_file, 1, {'nbins'}, char('D'), D_OBS);
                f = write_netcdf (conc_file, 1, {'nlev'}, char('lev'), clev(1:nlev));
                f = write_netcdf (conc_file, 1, {'no_of_model_bins'}, char('N'), model_Nbin);

              end %loc_biCorr
              clear loc_biCorr;

            else
              % disp([m nmodel ilon nlon ilat nlat iseas nseas])
              disp('constrain_biasCorr_dust_fraction_clim_seas_3d: Already calculated and stored.....Proceed to copy.')

              %Copy from the stored file here.
              %     e.g.
              if (~save_data_only_netcdf)
                N = ncread(conc_file,'N');
                constrained_dust_PSD(:,ilon,ilat,:,iseas,m)  = ncread(conc_file,'constrained_dust_PSD');
                biCorr_dust_fraction_norm(1:N,ilon,ilat,:,iseas,m)  = ncread(conc_file,'biCorr_dust_fraction_norm');
                param_min_latlon(:,ilon,ilat,:,iseas,m)  = ncread(conc_file,'param_min_latlon');
                range_warning_latlon(:,ilon,ilat,:,iseas,m)  = ncread(conc_file,'range_warning_latlon');
              else
                constrained_dust_PSD = 0;
                param_min_latlon = 0;
                range_warning_latlon = 0;
                biCorr_dust_fraction_norm = 0;
              end%save_data_only_netcdf
            end

          end %ilat
        end %ilat


      end % iseas


    end % m

%   delete(gcp('nocreate')); % First close any existing pool
    clear model_column_loading_orig D_lower D_upper Nbin no_model_data

    disp(' ');

  end

  % disp('Done...')
  % Done..

%   end
