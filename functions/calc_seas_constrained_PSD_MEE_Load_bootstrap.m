function [PSD_constrained,MEE_lev_constrained,MEE_constrained,Load_constrained] = calc_seas_constrained_PSD_MEE_Load_bootstrap (indlon,indlat,model_name,nboot,recalculate)

% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This code calculates the seasonally-averaged bootstrap PSD, MEE and
% dust loading for DustCOMM.
% ==========

% ==========
% Steps:
% 1. Calculate or collect the boostrap PSD.
% 2. Use that to calculate the MEE
% 3. Use the MEE to calculate the dust loading

% =================
% Begin here....
% =================

%  define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          mean_rho_d error_rho_d...
          gQe_OBS D_OBS gPSD_OBS gPSD_OBS_median ...
          DOAD_data_seas dustcomm_dir_seas seas_filename ....


 % make the folder if it is not available
  system(char(strcat('/bin/mkdir -p',{' '},dustcomm_dir_seas))); % make the directory if not available

% general names of seasons
  seas_str = {'DJF', 'MAM', 'JJA','SON'}; % seasons

% =================
% get specific longitude/latitude based on the specified index
% =================
% get specific longitude and latitude
  xlat = clat(indlat);
  xlon = clon(indlon);
  clon_str = sprintf('%0.2f',xlon); % the strings
  clat_str = sprintf('%0.2f',xlat);

% Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

  % convert the diameter to meters
  D_m = 10^(-6)*D_OBS; %in meters

  % =================
% dimension
% =================
    nlev = max(size(indlev)); %
    nD = max(size(D_OBS));
    nmodel = max(size(model_name)); %
    nseas = max(size(seas_str));
    nparam = 6;

% =================
% check if the data already exist
% =================

  conc_file = char(strcat(dustcomm_dir_seas,seas_filename,'_',clon_str,'_',clat_str,'.nc'));

  if (~exist(conc_file, 'file') | recalculate)

% =================
% Fit the analytical equation for model simlation of the seasons
% and calculate the boostrap for that location
% =================
    disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Constraining the size distribution for all the models');
    [fitPSD_constrained,~] = constrain_biasCorr_dust_fraction_clim_seas_3d (model_name,[xlon, xlat],'seas',false);
    fitPSD_constrained = squeeze(fitPSD_constrained);

    if (any(~isnan(fitPSD_constrained(:))))

      disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Get the constrained annually-averaged size distribution to scale the seasonal values...');
      % ; get the median of the annually-averaged  for the bias-corrected boostratp for each location
      [annu_bootstrap] = collect_annu_constrained_PSD_bootstrap([xlon, xlat],nboot,true);
      xannu_bootstrap = repelem(annu_bootstrap,1,1,nmodel); % repeat for all the models

      % averaged this fit of each season to create annually-average for each
      % model and the location
      fitPSD_constrained_ann = squeeze(mean(fitPSD_constrained(:,indlev,:,:),3,'omitnan'));

      % scale the new fit back to the bootstrap mean
      scale_f = xannu_bootstrap./fitPSD_constrained_ann;
      scale_f = reshape(scale_f,nD,nlev,1,nmodel);
      scale_f = repelem(scale_f,1,1,nseas,1);

      % correct here
      xfitPSD_constrained = squeeze(fitPSD_constrained(:,indlev,:,:)).*scale_f;
      fitPSD_constrained = xfitPSD_constrained;

      clear scale_f annu_bootstrap xannu_bootstrap fitPSD_constrained_ann xfitPSD_constrained

      disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Renormalize to unity...')
    % first re-normalize all the distribution for each location
      xsum = trapz(D_OBS,fitPSD_constrained,1);
      xsum = repelem(xsum,nD,1,1,1,1);
      fitPSD_constrained = fitPSD_constrained./xsum;
      clear xsum;

  % =================
  % Get the modified dust AOD for all the season
  % =================
      disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Getting the dust AOD...')
      DAOD_m  = ncread(DOAD_data_seas,'DAOD_m');
      DAOD_sd = ncread(DOAD_data_seas,'DAOD_sd');
      DAOD_m = DAOD_m(:,indlon,indlat);
      DAOD_sd = DAOD_sd(:,indlon,indlat);

  %  =================
      disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Getting the dust vertical weight.....');
      evalc('[dust_weight] = get_model_dust_Vweight_data(''seas'',model_name,false)');
      dust_weight = squeeze(dust_weight(:,indlon,indlat,indlev,:));


  % =================
  % Now calculate the the bootstrap for all
  % Also calculate the the extinction and loading using the boostrap
  % =================
      biCorr_bootstrap = NaN(nD,nlev,nseas,nboot);
      MEE = NaN(nboot,nlev,nseas);
      mass_ext = NaN(nboot,nseas);
      Latm = NaN(nboot,nseas);

  % =================
      disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Calculating the bootstrap dust PSD, MEE and Load....')
  % =================

      for j=1:nboot
        if (mod(j,10)==0)
          disp(['Bootstrap realization for j=',sprintf('%05d',j)])
        end

        rand_j = ceil(nboot*rand(1));
        rand_m = ceil(nmodel*rand(1));

        % Scale the bias-corrected distribution based on ramdonly selected PSD
        scale_constant = gPSD_OBS(rand_j,:)./gPSD_OBS_median;
        scale_constant = reshape(scale_constant,max(size(scale_constant)),1);
        xscale_constant = repelem(scale_constant,1,nlev,nseas);
        boot_PSD = fitPSD_constrained(:,indlev,:,rand_m).*xscale_constant;
        biCorr_bootstrap(:,:,:,j) = boot_PSD;

        clear xscale_constant scale_constant

        % get other randonly selected variable
        Qe_rand = gQe_OBS(rand_j,:);
        Qe_rand_lev = repelem(Qe_rand',1,nlev);
        D_m_lev = repelem(D_m',1,nlev);
        rho_d = Gaussian(mean_rho_d,error_rho_d);

        for iseas=1:nseas

          if (any(all(~isnan(boot_PSD(:,:,iseas)))))

            dV_dD = 10^6*squeeze(boot_PSD(:,:,iseas));
            dV_dD_QD = (3/2)*(1/rho_d)*dV_dD.*(Qe_rand_lev./D_m_lev);
            MEE(j,:,iseas) = 10^-3*trapz(D_m,dV_dD_QD,1); %in m2/g

          % verically averged value
            MEE_weight = MEE(j,:,iseas).*squeeze(dust_weight(rand_m,:,iseas));
            mass_ext(j,iseas) = squeeze(sum(MEE_weight,2,'omitnan'));
            DAOD = Gaussian(DAOD_m(iseas), DAOD_sd(iseas));
            Latm(j,iseas) = (DAOD./mass_ext(j,iseas)); % in g/m2

            clear dV_dD dV_dD_QD MEE_weight DAOD

          end
        end

        clear Qe_rand Qe_rand_lev D_m_lev

      end  % j

      clear fitPSD_constrained
%       clear gPSD_OBS gPSD_OBS_median gQe_OBS

  % =================
    disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Now calculate the uncertainty...');
  % =================

    PSD_constrained = NaN(nparam,nD,nlev,nseas);
    MEE_lev_constrained = NaN(nparam,nlev,nseas);
    MEE_constrained = NaN(nparam,nseas);
    Load_constrained = NaN(nparam,nseas);

    for iseas=1:nseas
      for ilev=1:nlev
        for i=1:nD
  %PSD
          A=squeeze(sort(biCorr_bootstrap(i,ilev,iseas,:)));
          if (~all(isnan(A)))
            A = A(~isnan(A));
            xnboot = max(size(A));
            PSD_constrained(1,i,ilev,iseas) = mean(A,1,'omitnan'); % mean vale
            PSD_constrained(2,i,ilev,iseas) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
            PSD_constrained(3,i,ilev,iseas) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
            PSD_constrained(4,i,ilev,iseas) = A(max(1,(floor(0.5*xnboot)))); % median  value
            PSD_constrained(5,i,ilev,iseas) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
            PSD_constrained(6,i,ilev,iseas) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
          end

        end % nD
  %MEE with levels
        A=squeeze(sort(MEE(:,ilev,iseas)));
        if (~all(isnan(A)))
          A = A(~isnan(A));
          xnboot = max(size(A));
          MEE_lev_constrained(1,ilev,iseas) = mean(A,1,'omitnan'); % mean vale
          MEE_lev_constrained(2,ilev,iseas) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
          MEE_lev_constrained(3,ilev,iseas) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
          MEE_lev_constrained(4,ilev,iseas) = A(max(1,(floor(0.5*xnboot)))); % median  value
          MEE_lev_constrained(5,ilev,iseas) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
          MEE_lev_constrained(6,ilev,iseas) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
        end

      end % nlev

  %MEE without levels

      A=squeeze(sort(mass_ext(:,iseas)));
      if (~all(isnan(A)))
        A = A(~isnan(A));
        xnboot = max(size(A));
        MEE_constrained(1,iseas) = mean(A,1,'omitnan'); % mean vale
        MEE_constrained(2,iseas) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
        MEE_constrained(3,iseas) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
        MEE_constrained(4,iseas) = A(max(1,(floor(0.5*xnboot)))); % median  value
        MEE_constrained(5,iseas) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
        MEE_constrained(6,iseas) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
      end

  % Load
      A=squeeze(sort(Latm(:,iseas)));
      if (~all(isnan(A)))
        A = A(~isnan(A));
        xnboot = max(size(A));
        Load_constrained(1,iseas) = mean(A,1,'omitnan'); % mean vale
        Load_constrained(2,iseas) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
        Load_constrained(3,iseas) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
        Load_constrained(4,iseas) = A(max(1,(floor(0.5*xnboot)))); % median  value
        Load_constrained(5,iseas) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
        Load_constrained(6,iseas) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
      end

    end %  nseas

      clear biCorr_bootstrap MEE mass_ext Latm

    else

  % =================
    disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: All values are missing...NaN will be recorded');
  % =================
      PSD_constrained = NaN(nparam,nD,nlev,nseas);
      MEE_lev_constrained = NaN(nparam,nlev,nseas);
      MEE_constrained = NaN(nparam,nseas);
      Load_constrained = NaN(nparam,nseas);
    end

% =================
    disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Now saving the data....');
% =================

    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 4, {'param','D','lev','seas'}, char('PSD_constrained'), PSD_constrained);
    f = write_netcdf (conc_file, 3, {'param','lev','seas'}, char('MEE_lev_constrained'), MEE_lev_constrained);
    f = write_netcdf (conc_file, 2, {'param','seas'}, char('MEE_constrained'), MEE_constrained);
    f = write_netcdf (conc_file, 2, {'param','seas'}, char('Load_constrained'), Load_constrained);

    f = write_netcdf (conc_file, 1, {'param'}, char('param'), [0,-2,-1,0,1,2]);
    f = write_netcdf (conc_file, 1, {'D'}, char('D'), D_OBS);
    f = write_netcdf (conc_file, 1, {'lev'}, char('lev'), clev(indlev));
    f = write_netcdf (conc_file, 1, {'seas'}, char('seas'), [12,3,6,9]);

  else
    % =================
    disp('calc_seas_constrained_PSD_MEE_Load_bootstrap: Data already saved...');
    disp('Now get it...');
    % =================

    PSD_constrained = ncread(conc_file,'PSD_constrained');
    MEE_lev_constrained = ncread(conc_file,'MEE_lev_constrained');
    MEE_constrained = ncread(conc_file,'MEE_constrained');
    Load_constrained = ncread(conc_file,'Load_constrained');

  end  % recalculate

%   disp('Done...')

end
