function [exitflag] = calc_annu_constrained_MEE_Load_bootstrap (model_name,nboot,recalculate)
% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This code constraints the annually-averaged mass extinction efficiency
% and dust load,using the constrained annually-averaged PSD and dust extinction effiency.
% ==========

% ==========
% Steps:
% 1. Collect the constrained annually-averaged PSD over each location.
% 2. Use that to calculate the constrained MEE and dust load

% =================
% Begin here....
% =================

%  define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          mean_rho_d error_rho_d DOAD_data_annu ...
          gQe_OBS D_OBS dustcomm_dir_mee_annu mee_filename_annu...
          use_desktop_computer use_parallel_pool_in_function nThreads_Max

%  some usefule variables.
  out_dir = dustcomm_dir_mee_annu;
  base_filename = mee_filename_annu;
  system(char(strcat('/bin/mkdir -p',{' '},out_dir))); % make the directory if not available

% zerise Exit flag...
  exitflag = 0;

% Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

% convert the diameter to meters
  D_m = 10^(-6)*D_OBS; %in meters

% =================
  %dimension
  nlev = max(size(indlev)); %
  nlat = max(size(clat)); %
  nlon = max(size(clon)); %

% ; for the many matlab workers
  xclat = clat;
  xclon = clon;
  Qe = gQe_OBS;
  mean_rho = mean_rho_d;
  error_rho = error_rho_d;

% =================
  disp('calc_annu_constrained_MEE_Load_bootstrap: Getting the dust vertical weights ...');

  [model_norm_profile] = get_model_dust_Vweight_data('annu',model_name,false);
  model_norm_profile = squeeze(mean(model_norm_profile(:,:,:,indlev),1,'omitnan'));

  % =================
  disp('calc_annu_constrained_MEE_Load_bootstrap: Getting the dust AOD...')
  DAOD_m  = ncread(DOAD_data_annu,'DAOD_m');
  DAOD_sd = ncread(DOAD_data_annu,'DAOD_sd');

% =================
  disp('calc_annu_constrained_MEE_Load_bootstrap: Loop over each location...')

  for ilon=1:nlon

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

    parfor ilat=1:nlat
%     for ilat=1:nlat

      [ilon nlon ilat nlat]
%Some parameter that will be needed in all workstation
      xlat = xclat;
      xlon = xclon
      clon_str = sprintf('%0.2f',xlon(ilon));
      clat_str = sprintf('%0.2f',xlat(ilat));
      conc_file = char(strcat(out_dir,base_filename,'_',clon_str,'_',clat_str,'.nc'));

      if (~exist(conc_file, 'file') || recalculate)
% Now First get the bootstrap PSD
        [constrained_PSD] = collect_annu_constrained_PSD_bootstrap ([xlon(ilon),xlat(ilat)],nboot,false);

        MEE = NaN(nboot,nlev);
        Latm = NaN(1,nboot);
        Qext = Qe;

%calculate MEE
        for j=1:nboot
          rand_j = ceil(nboot*rand(1));
          if (any(all(~isnan(constrained_PSD(:,:,rand_j)))))
            %Radomly choose from the Qe
            Qe_rand = Qext(rand_j,:);
            rho_d = Gaussian(mean_rho,error_rho);
            dV_dD = 10^6*squeeze(constrained_PSD(:,:,rand_j));
            Qe_rand_lev = repelem(Qe_rand',1,nlev);
            D_m_lev = repelem(D_m',1,nlev);
            dV_dD_QD = (3/2)*(1/rho_d)*dV_dD.*(Qe_rand_lev./D_m_lev);
            MEE(j,:) = 10^-3*trapz(D_m,dV_dD_QD,1); %in m2/g

          % verically averged value
            MEE_weight = MEE(j,:).*squeeze(model_norm_profile(ilon,ilat,:))';
            mass_ext = squeeze(sum(MEE_weight,2,'omitnan'));
            DAOD = Gaussian(DAOD_m(ilon,ilat), DAOD_sd(ilon,ilat));
            Latm(j) = (DAOD./mass_ext); % in g/m2
          end
        end % j

        x = delete_file_ifexist(conc_file);
        f = write_netcdf (conc_file, 2, {'nboot','nlev'}, char('MEE'), MEE);
        f = write_netcdf (conc_file, 1, {'nlev'}, char('dust_weight'), squeeze(model_norm_profile(ilon,ilat,:)));
        f = write_netcdf (conc_file, 1, {'nboot'}, char('Latm'), Latm);
        exitflag = 1;
      else
        disp('calc_annu_constrained_MEE_Load_bootstrap: Data already exist...')
        exitflag = 0;
      end % is the data already available
    end  % ilat

    delete(gcp('nocreate')); % First close any existing pool
  end    % ilon

%   disp('Done...')

% end
