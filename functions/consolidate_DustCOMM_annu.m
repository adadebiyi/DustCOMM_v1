function consolidate_DustCOMM_annu (nboot,recalculate)

% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This code consolidates the DustCOMM annually-averaged data.
% ==========

% =================
% Begin here....
% =================

%  define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          use_desktop_computer use_parallel_pool_in_function nThreads_Max ...
          D_OBS ...
          dustcomm_dir_all dustcomm_dir_mee_annu mee_filename_annu

% Create the folder if it is not available
  system(char(strcat('/bin/mkdir -p',{' '},dustcomm_dir_all))); % make the directory if not available

% Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

  % =================
  %dimension
  nlev = max(size(indlev)); %
  nlat = max(size(clat)); %
  nlon = max(size(clon)); %
  nD = max(size(D_OBS));

% For many matalab workers
  xclat = clat;
  xclon = clon;

%% =================
%  PSD
% ----------------
  base_filename = 'Dust_Size_Distr_dVdD_annual';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
      %   where to save it
    dVdD_mean = NaN(nD,nlev,nlat,nlon);
    dVdD_median = dVdD_mean;
    dVdD_Pos1sig = dVdD_mean;
    dVdD_Neg1sig = dVdD_mean;
    dVdD_Pos2sig = dVdD_mean;
    dVdD_Neg2sig = dVdD_mean;

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


    for ilon=1:nlon


      if ((nThreads > 1) && ~use_desktop_computer)
        pool = parpool(nThreads,'IdleTimeout', Inf); %Start the parallel pool using the 'local' profile.
      %           The pool is not defined here for desktop computer. When Parfor
      %           command is called later, the code defaults to the standard
      %           number of workers allowed on the computer.
      end

      parfor ilat=1:nlat

        disp(['PSD ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])
  %Some parameter that will be needed in all workstation
        xlat = xclat;
        xlon = xclon

      % Now First get the bootstrap PSD
      [fitPSD_bootstrap] = collect_annu_constrained_PSD_bootstrap ([xlon(ilon),xlat(ilat)],nboot,false);

      if (~all(isnan(fitPSD_bootstrap(:))))
        for i=1:nD
          for ilev=1:nlev
          A=squeeze(sort(fitPSD_bootstrap(i,ilev,:)));
          if (~all(isnan(A)))
            A = A(~isnan(A));
            xnboot = max(size(A));
            dVdD_mean(i,ilev,ilat,ilon) = mean(A,1,'omitnan'); % mean vale
            dVdD_Neg2sig(i,ilev,ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
            dVdD_Neg1sig(i,ilev,ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
            dVdD_median(i,ilev,ilat,ilon) = A(max(1,(floor(0.5*xnboot)))); % median  value
            dVdD_Pos1sig(i,ilev,ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
            dVdD_Pos2sig(i,ilev,ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value
          end
          end % nlev
        end % nD
      end

      end % ilat

      delete(gcp('nocreate')); % First close any existing pool

    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_mean'), dVdD_mean);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_median'), dVdD_median);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_Pos1sig'), dVdD_Pos1sig);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_Neg1sig'), dVdD_Neg1sig);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_Pos2sig'), dVdD_Pos2sig);
    f = write_netcdf (conc_file, 4, {'D','lev','lat','lon'}, char('dVdD_Neg2sig'), dVdD_Neg2sig);
    f = write_netcdf (conc_file, 1, {'D'}, char('D'), D_OBS);
    f = write_netcdf (conc_file, 1, {'lev'}, char('lev'), clev(indlev));
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);

  % Add the attributes
    ncwriteatt(conc_file,char('dVdD_mean'),char('description'),'Mean Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_mean'),char('unit'),'?');
    ncwriteatt(conc_file,char('dVdD_median'),char('description'),'Median Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_median'),char('unit'),'?');
    ncwriteatt(conc_file,char('dVdD_Pos1sig'),char('description'),'+1 Sigma Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_Pos1sig'),char('unit'),'?');
    ncwriteatt(conc_file,char('dVdD_Neg1sig'),char('description'),'-1 Sigma Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_Neg1sig'),char('unit'),'?');
    ncwriteatt(conc_file,char('dVdD_Pos2sig'),char('description'),'+2 Sigma Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_Pos2sig'),char('unit'),'?');
    ncwriteatt(conc_file,char('dVdD_Neg2sig'),char('description'),'-2 Sigma Normalized Dust Volume Size Distribution');
    ncwriteatt(conc_file,char('dVdD_Neg2sig'),char('unit'),'?');
    ncwriteatt(conc_file,char('D'),char('description'),'Dust Geometric Diameter');
    ncwriteatt(conc_file,char('D'),char('unit'),'microns');
    ncwriteatt(conc_file,char('lev'),char('description'),'Pressure levels');
    ncwriteatt(conc_file,char('lev'),char('unit'),'hPa');
    ncwriteatt(conc_file,char('lat'),char('description'),'Latitude');
    ncwriteatt(conc_file,char('lat'),char('unit'),'degrees_north');
    ncwriteatt(conc_file,char('lon'),char('description'),'Longitude');
    ncwriteatt(conc_file,char('lon'),char('unit'),'degrees_east');
    ncwriteatt(conc_file,'/','creation_date',datestr(now));

    clear dVdD_mean dVdD_median dVdD_Pos1sig dVdD_Neg1sig dVdD_Pos2sig dVdD_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist

    disp('  ');




%% =================
%  3D MEE
% ----------------
  base_filename = 'Dust_3D_MEE_annual';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
%   where to save it
    MEE_mean = NaN(nlev,nlat,nlon);
    MEE_median = MEE_mean;
    MEE_Pos1sig = MEE_mean;
    MEE_Neg1sig = MEE_mean;
    MEE_Pos2sig = MEE_mean;
    MEE_Neg2sig = MEE_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['3-D MEE ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])
%         get the file
        filename = char(strcat(dustcomm_dir_mee_annu,mee_filename_annu,...
          '_',sprintf('%0.2f',clon(ilon)),'_',sprintf('%0.2f',clat(ilat)),'.nc'));
        MEE = ncread(filename,'MEE');

      if (~all(isnan(MEE(:))))
        for ilev=1:nlev
          A=squeeze(sort(MEE(:,ilev)));
          if (~all(isnan(A)))
            A = A(~isnan(A));
            A = A(A~=0);
            xnboot = max(size(A));
            MEE_mean(ilev,ilat,ilon) = mean(A,1,'omitnan'); % mean vale
            MEE_Neg2sig(ilev,ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
            MEE_Neg1sig(ilev,ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
            MEE_median(ilev,ilat,ilon) = A(max(1,(floor(0.5*xnboot)))); % median  value
            MEE_Pos1sig(ilev,ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
            MEE_Pos2sig(ilev,ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value
          end
        end % nlev
      end

      end % ilat
    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_mean'), MEE_mean);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_median'), MEE_median);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_Pos1sig'), MEE_Pos1sig);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_Neg1sig'), MEE_Neg1sig);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_Pos2sig'), MEE_Pos2sig);
    f = write_netcdf (conc_file, 3, {'lev','lat','lon'}, char('MEE_Neg2sig'), MEE_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lev'}, char('lev'), clev(indlev));
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);

  % Add the attributes
    ncwriteatt(conc_file,char('MEE_mean'),char('description'),'Mean Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_mean'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_median'),char('description'),'Median Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_median'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Pos1sig'),char('description'),'+1 Sigma Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Pos1sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Neg1sig'),char('description'),'-1 Sigma Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Neg1sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Pos2sig'),char('description'),'+2 Sigma Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Pos2sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Neg2sig'),char('description'),'-2 Sigma Dust 3D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Neg2sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('lev'),char('description'),'Pressure levels');
    ncwriteatt(conc_file,char('lev'),char('unit'),'hPa');
    ncwriteatt(conc_file,char('lat'),char('description'),'Latitude');
    ncwriteatt(conc_file,char('lat'),char('unit'),'degrees_north');
    ncwriteatt(conc_file,char('lon'),char('description'),'Longitude');
    ncwriteatt(conc_file,char('lon'),char('unit'),'degrees_east');
    ncwriteatt(conc_file,'/','creation_date',datestr(now));

    clear MEE_mean MEE_median MEE_Pos1sig MEE_Neg1sig MEE_Pos2sig MEE_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist
  disp('  ');




%% =================
%  2D MEE
% ----------------
  base_filename = 'Dust_2D_MEE_annual';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
      %   where to save it
    MEE_mean = NaN(nlat,nlon);
    MEE_median = MEE_mean;
    MEE_Pos1sig = MEE_mean;
    MEE_Neg1sig = MEE_mean;
    MEE_Pos2sig = MEE_mean;
    MEE_Neg2sig = MEE_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['2-D MEE ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])
%         get the file
        filename = char(strcat(dustcomm_dir_mee_annu,mee_filename_annu,...
          '_',sprintf('%0.2f',clon(ilon)),'_',sprintf('%0.2f',clat(ilat)),'.nc'));
        MEE = ncread(filename,'MEE');
        dust_weight = ncread(filename,'dust_weight');
        dust_weight = repelem(dust_weight',nboot,1);
        MEE = sum(MEE.*dust_weight,2,'omitnan');
        MEE = MEE(MEE~=0);

        if (~all(isnan(MEE(:))))
          A=squeeze(sort(MEE(~isnan(MEE))));
          A = A(A~=0);
          xnboot = max(size(A));
          MEE_mean(ilat,ilon) = mean(A,1,'omitnan'); % mean vale
          MEE_Neg2sig(ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
          MEE_Neg1sig(ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
          MEE_median(ilat,ilon) = A(max(1,(floor(0.5*xnboot)))); % median  value
          MEE_Pos1sig(ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
          MEE_Pos2sig(ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value
        end
      end % ilat
    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_mean'), MEE_mean);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_median'), MEE_median);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_Pos1sig'), MEE_Pos1sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_Neg1sig'), MEE_Neg1sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_Pos2sig'), MEE_Pos2sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('MEE_Neg2sig'), MEE_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);

  % Add the attributes
    ncwriteatt(conc_file,char('MEE_mean'),char('description'),'Mean Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_mean'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_median'),char('description'),'Median Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_median'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Pos1sig'),char('description'),'+1 Sigma Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Pos1sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Neg1sig'),char('description'),'-1 Sigma Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Neg1sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Pos2sig'),char('description'),'+2 Sigma Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Pos2sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('MEE_Neg2sig'),char('description'),'-2 Sigma Dust 2D Mass Extinction Efficiency');
    ncwriteatt(conc_file,char('MEE_Neg2sig'),char('unit'),'m2/g');
    ncwriteatt(conc_file,char('lat'),char('description'),'Latitude');
    ncwriteatt(conc_file,char('lat'),char('unit'),'degrees_north');
    ncwriteatt(conc_file,char('lon'),char('description'),'Longitude');
    ncwriteatt(conc_file,char('lon'),char('unit'),'degrees_east');
    ncwriteatt(conc_file,'/','creation_date',datestr(now));

    clear MEE_mean MEE_median MEE_Pos1sig MEE_Neg1sig MEE_Pos2sig MEE_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist
  disp('  ');


  %% =================
%  Load
% ----------------
  base_filename = 'Dust_Load_annual';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
      %   where to save it
    Latm_mean = NaN(nlat,nlon);
    Latm_median = Latm_mean;
    Latm_Pos1sig = Latm_mean;
    Latm_Neg1sig = Latm_mean;
    Latm_Pos2sig = Latm_mean;
    Latm_Neg2sig = Latm_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['2-D Load ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])
%         get the file
        filename = char(strcat(dustcomm_dir_mee_annu,mee_filename_annu,...
          '_',sprintf('%0.2f',clon(ilon)),'_',sprintf('%0.2f',clat(ilat)),'.nc'));
        Latm = ncread(filename,'Latm');
        Latm = Latm(Latm~=0);

        if (~all(isnan(Latm(:))))
          A=squeeze(sort(Latm(~isnan(Latm))));
          A = A(A~=0);
          xnboot = max(size(A));
          Latm_mean(ilat,ilon) = mean(A,1,'omitnan'); % mean vale
          Latm_Neg2sig(ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*xnboot)))); %lower uncertainty value at -2 sigma
          Latm_Neg1sig(ilat,ilon) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*xnboot)))); %lower uncertainty value at -1 sigma
          Latm_median(ilat,ilon) = A(max(1,(floor(0.5*xnboot)))); % median  value
          Latm_Pos1sig(ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*xnboot)))); %upper uncertainty value at +1 sigma
          Latm_Pos2sig(ilat,ilon) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*xnboot)))); %upper uncertainty value
        end
      end % ilat
    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_mean'), Latm_mean);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_median'), Latm_median);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_Pos1sig'), Latm_Pos1sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_Neg1sig'), Latm_Neg1sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_Pos2sig'), Latm_Pos2sig);
    f = write_netcdf (conc_file, 2, {'lat','lon'}, char('Latm_Neg2sig'), Latm_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);

  % Add the attributes
    ncwriteatt(conc_file,char('Latm_mean'),char('description'),'Mean Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_mean'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Latm_median'),char('description'),'Median Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_median'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Latm_Pos1sig'),char('description'),'+1 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_Pos1sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Latm_Neg1sig'),char('description'),'-1 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_Neg1sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Latm_Pos2sig'),char('description'),'+2 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_Pos2sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Latm_Neg2sig'),char('description'),'-2 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Latm_Neg2sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('lat'),char('description'),'Latitude');
    ncwriteatt(conc_file,char('lat'),char('unit'),'degrees_north');
    ncwriteatt(conc_file,char('lon'),char('description'),'Longitude');
    ncwriteatt(conc_file,char('lon'),char('unit'),'degrees_east');
    ncwriteatt(conc_file,'/','creation_date',datestr(now));

    clear Latm_mean Latm_median Latm_Pos1sig Latm_Neg1sig Latm_Pos2sig Latm_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist
  disp('  ');

  % disp('Done...');


%   end
