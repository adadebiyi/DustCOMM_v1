function calc_annu_constrained_PSD_bootstrap (model_name,nboot,recalculate)
% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles
% ==========
% This function uses the constrained data calculated for each model and for each location using the
% mean globally-averaged observationally-constrained PSD to create the
% ==========
% Steps:
% 1. Collect the constrained PSD data for all the models
% 2. Randomly scale the mean distribution based on the realizations of the
% globally-averaged PSD from Kok et al, 2017
% 3. Use that with the constrained PSD dat randonly selecting from
% one of the models

% ==========
% Variables:
% model_name --  models to use
% annu_or_seas -- state whether annual or seasonal

% Begin here....
% =================
% define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          gPSD_OBS gPSD_OBS_median D_OBS ...
          dustcomm_dir_psd_annu fitPSD_filename_annu ...
          use_parallel_pool_in_function nThreads_Max use_desktop_computer

  disp('calc_annu_constrained_PSD_bootstrap: Calculate most-likely annually-averaged constrained PSD for each location...')
% =================
% Create the directory to store...if not available.
  system(char(strcat('/bin/mkdir -p',{' '},dustcomm_dir_psd_annu))); % make the directory if not available

%   ;; Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

% =================
  %dimension
  nlev_all = max(size(clev)); %
  nlev = max(size(indlev)); %
  nlat = max(size(clat)); %
  nlon = max(size(clon)); %
  nD = size(D_OBS,2);
  nmodel = max(size(model_name)); %

%   copy out the lon/lat not to use in par loop
  xclon = clon;
  xclat = clat;

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


  if ((nThreads > 1) && ~use_desktop_computer)
          pool = parpool(nThreads,'IdleTimeout', Inf); %Start the parallel pool using the 'local' profile.
%           The pool is not defined here for desktop computer. When Parfor
%           command is called later, the code defaults to the standard
%           number of workers allowed on the computer.
  end

% -----------
% %now get all the fitted PSD for all the model
  disp('calc_annu_constrained_PSD_bootstrap: Getting the constrained size distribution for all the models');
  mfit_PSD = NaN(nD,nlat,nlon,nlev_all,nmodel);
  for ilon=1:nlon
    disp([ilon nlon])
    parfor ilat=1:nlat
      [mfit_PSD(:,ilat,ilon,:,:),~] = constrain_biasCorr_dust_fraction_clim_seas_3d (model_name,[xclon(ilon), xclat(ilat)],'annu',false);
    end
  end

  delete(gcp('nocreate')); % First close any existing pool
% -----------

 % ; correct all the KokFIt distribution by the globally averaged distribution
  disp(' ');
  disp('calc_annu_constrained_PSD_bootstrap: Collecting the model dust loading...');
%   evalc('[model_column_loading_orig,~] = get_model_dust_concentration_all_3d(''annu'',model_name,false)');
  evalc('[model_column_loading_orig,~] = get_model_dust_concentration_data(''annu'',model_name)');
  model_load = squeeze(mean(mean(sum(model_column_loading_orig,2,'omitnan'),3,'omitnan'),4,'omitnan'));
  model_load = model_load./repelem(sum(model_load,2,'omitnan'),1,nlev_all);
  clear model_column_loading_orig;

  disp('calc_annu_constrained_PSD_bootstrap: Correcting the globally-averaged fit to get rid of random errors...');
  glob_mfit_PSD = squeeze(mean(mean(mfit_PSD,2,'omitnan'),3,'omitnan'));
  for m=1:nmodel
    Kokfit = sum(glob_mfit_PSD(:,:,m).*repelem(model_load(m,:),nD,1),2,'omitnan');
    Kokfit = reshape(Kokfit,1,nD);
%     scale_f = gPSD_OBS_median./(D_OBS.*Kokfit);
    scale_f = (gPSD_OBS_median./D_OBS)./Kokfit;
    rep_scale_f = repelem(scale_f',1,nlat,nlon,nlev_all);
    mfit_PSD(:,:,:,:,m) = mfit_PSD(:,:,:,:,m).*rep_scale_f;
    clear Kokfit scale_f rep_scale_f
  end

  clear glob_mfit_PSD model_load

   disp('calc_annu_constrained_PSD_bootstrap: Renormalize to unity...')
  % first re-normalize all the distribution for each location
  xsum = trapz(D_OBS,mfit_PSD,1);
  xsum = repelem(xsum,nD,1,1,1,1);
  mfit_PSD = mfit_PSD./xsum;
  clear xsum;

%       Now calculate the the bootstrap
% define filename
  check_dim = (size(mfit_PSD,1) == size(gPSD_OBS,2)) && ...
            (size(mfit_PSD,1) == max(size(gPSD_OBS_median)));

  if (check_dim)
    for j=1:nboot
      % for j=457:nboot
%  Check if the output file exist for every location of the globe
      j_str = sprintf('%05d',j);
      out_dir_j = strcat(dustcomm_dir_psd_annu,j_str,'/');
      out_dir_j_basefile = char(strcat(out_dir_j,fitPSD_filename_annu));
      system(char(strcat('/bin/mkdir -p',{' '},out_dir_j))); % make the directory if not available
      conc_file = char(strcat(out_dir_j_basefile,'_*_',j_str,'.nc'));
      listing = dir(conc_file);
      if (numel(listing) ~= nlat*nlon || recalculate)
        if (mod(j,10)==0)
            disp(['Bootstrap realization for j=',sprintf('%05d',j)])
        end
%         [exitflag] = create_biCorr_bootstrap (j_str,nboot,gPSD_OBS,gPSD_OBS_median,mfit_PSD,clon,clat,out_dir_j_basefile);

% Random numbers to use
        rand_j = ceil(nboot*rand(1));
        rand_m = ceil(nmodel*rand(1));

% Scale the constrained distribution based on ramdonly selected globally-averaged PSD
        scale_constant = gPSD_OBS(rand_j,:)./gPSD_OBS_median;
        scale_constant = reshape(scale_constant,max(size(scale_constant)),1);
        xscale_constant = repelem(scale_constant,1,nlat,nlon,nlev);
        fit_PSD_bootstrap = mfit_PSD(:,:,:,indlev,rand_m).*xscale_constant;

%             store the data here
        for ilat=1:nlat
          for ilon=1:nlon
            clon_str = sprintf('%0.2f',clon(ilon));
            clat_str = sprintf('%0.2f',clat(ilat));
            conc_file = char(strcat(out_dir_j,fitPSD_filename_annu,'_',clon_str,'_',clat_str,'_',j_str,'.nc'));
            x = delete_file_ifexist(conc_file);
            f = write_netcdf (conc_file, 2, {'nD','nlev'}, char(fitPSD_filename_annu), squeeze(fit_PSD_bootstrap(:,ilat,ilon,:)));
          end %ilon
        end % ilat
      else
        disp(['Files already exist for j=',j_str])
      end % conc_file

    end  % j
    clear xscale_constant scale_constant
  else
    warning('The dimensions are not in required order')
  end %check_dim

  clear mfit_PSD fit_PSD_bootstrap
  % -----------

  % disp('Done...')
  % Done..
% end
