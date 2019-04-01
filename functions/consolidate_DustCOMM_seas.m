function consolidate_DustCOMM_seas (recalculate)

% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This code consolidates the DustCOMM seasonally-averaged data.
% ==========

% =================
% Begin here....
% =================

%  define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          D_OBS dustcomm_dir_seas seas_filename...
          dustcomm_dir_all


  % Create the folder if it is not available
  system(char(strcat('/bin/mkdir -p',{' '},dustcomm_dir_all))); % make the directory if not available

% Define the standard season names
  seas_str = {'DJF', 'MAM', 'JJA','SON'}; % seasons
  seas_num = [12, 1, 2; 3, 4, 5; 6, 7, 8; 9, 10, 11]; % The corresponding index

% Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

  % =================
  %dimensions
  nlev = max(size(indlev)); %
  nlat = max(size(clat)); %
  nlon = max(size(clon)); %
  nD = max(size(D_OBS));
  nseas = max(size(seas_str));

  %% =================
%  PSD
% ----------------
  base_filename = 'Dust_Size_Distr_dVdD_seasonal';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
%   where to save it
    dVdD_mean = NaN(nD,nlev,nlat,nlon,nseas);
    dVdD_median = dVdD_mean;
    dVdD_Pos1sig = dVdD_mean;
    dVdD_Neg1sig = dVdD_mean;
    dVdD_Pos2sig = dVdD_mean;
    dVdD_Neg2sig = dVdD_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['PSD ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])

        clon_str = sprintf('%0.2f',clon(ilon)); % the strings
        clat_str = sprintf('%0.2f',clat(ilat));

        filename = char(strcat(dustcomm_dir_seas,seas_filename,'_',clon_str,'_',clat_str,'.nc'));

        pdata = ncread(filename,'PSD_constrained');
        dVdD_mean(:,:,ilat,ilon,:) = pdata(1,:,:,:);
        dVdD_Neg2sig(:,:,ilat,ilon,:) = pdata(2,:,:,:);
        dVdD_Neg1sig(:,:,ilat,ilon,:) = pdata(3,:,:,:);
        dVdD_median(:,:,ilat,ilon,:) = pdata(4,:,:,:);
        dVdD_Pos1sig(:,:,ilat,ilon,:) = pdata(5,:,:,:);
        dVdD_Pos2sig(:,:,ilat,ilon,:) = pdata(6,:,:,:);
        clear pdata
      end % ilat

    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_mean'), dVdD_mean);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_median'), dVdD_median);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_Pos1sig'), dVdD_Pos1sig);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_Neg1sig'), dVdD_Neg1sig);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_Pos2sig'), dVdD_Pos2sig);
    f = write_netcdf (conc_file, 5, {'D','lev','lat','lon','nseas'}, char('dVdD_Neg2sig'), dVdD_Neg2sig);
    f = write_netcdf (conc_file, 1, {'D'}, char('D'), D_OBS);
    f = write_netcdf (conc_file, 1, {'lev'}, char('lev'), clev(indlev));
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);
    f = write_netcdf (conc_file, 2, {'nseas','nmonth'}, char('month_season'), seas_num);

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

    ncwriteatt(conc_file,char('month_season'),char('description'),'Months in each season');

    clear dVdD_mean dVdD_median dVdD_Pos1sig dVdD_Neg1sig dVdD_Pos2sig dVdD_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist

  disp('  ');

%% =================
%  3D MEE
% ----------------
  base_filename = 'Dust_3D_MEE_seasonal';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
      %   where to save it
    MEE_mean = NaN(nlev,nlat,nlon,nseas);
    MEE_median = MEE_mean;
    MEE_Pos1sig = MEE_mean;
    MEE_Neg1sig = MEE_mean;
    MEE_Pos2sig = MEE_mean;
    MEE_Neg2sig = MEE_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['3-D MEE ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])

        clon_str = sprintf('%0.2f',clon(ilon)); % the strings
        clat_str = sprintf('%0.2f',clat(ilat));

        filename = char(strcat(dustcomm_dir_seas,seas_filename,'_',clon_str,'_',clat_str,'.nc'));

        pdata = ncread(filename,'MEE_lev_constrained');
        MEE_mean(:,ilat,ilon,:) = pdata(1,:,:);
        MEE_Neg2sig(:,ilat,ilon,:) = pdata(2,:,:);
        MEE_Neg1sig(:,ilat,ilon,:) = pdata(3,:,:);
        MEE_median(:,ilat,ilon,:) = pdata(4,:,:);
        MEE_Pos1sig(:,ilat,ilon,:) = pdata(5,:,:);
        MEE_Pos2sig(:,ilat,ilon,:) = pdata(6,:,:);
        clear pdata
      end % ilat

    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_mean'), MEE_mean);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_median'), MEE_median);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_Pos1sig'), MEE_Pos1sig);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_Neg1sig'), MEE_Neg1sig);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_Pos2sig'), MEE_Pos2sig);
    f = write_netcdf (conc_file, 4, {'lev','lat','lon','nseas'}, char('MEE_Neg2sig'), MEE_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lev'}, char('lev'), clev(indlev));
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);
    f = write_netcdf (conc_file, 2, {'nseas','nmonth'}, char('month_season'), seas_num);

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

    ncwriteatt(conc_file,char('month_season'),char('description'),'Months in each season');

    clear MEE_mean MEE_median MEE_Pos1sig MEE_Neg1sig MEE_Pos2sig MEE_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist

  disp('  ');

  %% =================
%  2D MEE
% ----------------
  base_filename = 'Dust_2D_MEE_seasonal';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
%   where to save it
    MEE_mean = NaN(nlat,nlon,nseas);
    MEE_median = MEE_mean;
    MEE_Pos1sig = MEE_mean;
    MEE_Neg1sig = MEE_mean;
    MEE_Pos2sig = MEE_mean;
    MEE_Neg2sig = MEE_mean;

    for ilon=1:nlon
      for ilat=1:nlat

        disp(['2-D MEE ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])

        clon_str = sprintf('%0.2f',clon(ilon)); % the strings
        clat_str = sprintf('%0.2f',clat(ilat));

        filename = char(strcat(dustcomm_dir_seas,seas_filename,'_',clon_str,'_',clat_str,'.nc'));

        pdata = ncread(filename,'MEE_constrained');
        MEE_mean(ilat,ilon,:) = pdata(1,:);
        MEE_Neg2sig(ilat,ilon,:) = pdata(2,:);
        MEE_Neg1sig(ilat,ilon,:) = pdata(3,:);
        MEE_median(ilat,ilon,:) = pdata(4,:);
        MEE_Pos1sig(ilat,ilon,:) = pdata(5,:);
        MEE_Pos2sig(ilat,ilon,:) = pdata(6,:);
        clear pdata
      end % ilat

    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_mean'), MEE_mean);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_median'), MEE_median);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_Pos1sig'), MEE_Pos1sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_Neg1sig'), MEE_Neg1sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_Pos2sig'), MEE_Pos2sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('MEE_Neg2sig'), MEE_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);
    f = write_netcdf (conc_file, 2, {'nseas','nmonth'}, char('month_season'), seas_num);

  % Add the attributes
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

    ncwriteatt(conc_file,char('month_season'),char('description'),'Months in each season');

    clear MEE_mean MEE_median MEE_Pos1sig MEE_Neg1sig MEE_Pos2sig MEE_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist

  disp('  ');

%% =================
%  2D Load
% ----------------
  base_filename = 'Dust_Load_seasonal';
  conc_file = char(strcat(dustcomm_dir_all,base_filename,'.nc'));

  if (~exist(conc_file, 'file') || recalculate)
    disp(['Working on...',base_filename]);
      %   where to save it
    Load_mean = NaN(nlat,nlon,nseas);
    Load_median = Load_mean;
    Load_Pos1sig = Load_mean;
    Load_Neg1sig = Load_mean;
    Load_Pos2sig = Load_mean;
    Load_Neg2sig = Load_mean;

    for ilon=1:nlon
      for ilat=1:nlat
        disp(['2-D Load ', num2str(ilon), '  ',num2str(nlon), '  ',num2str(ilat), '  ',num2str(nlat) ])

        clon_str = sprintf('%0.2f',clon(ilon)); % the strings
        clat_str = sprintf('%0.2f',clat(ilat));

        filename = char(strcat(dustcomm_dir_seas,seas_filename,'_',clon_str,'_',clat_str,'.nc'));

        pdata = ncread(filename,'Load_constrained');
        Load_mean(ilat,ilon,:) = pdata(1,:);
        Load_Neg2sig(ilat,ilon,:) = pdata(2,:);
        Load_Neg1sig(ilat,ilon,:) = pdata(3,:);
        Load_median(ilat,ilon,:) = pdata(4,:);
        Load_Pos1sig(ilat,ilon,:) = pdata(5,:);
        Load_Pos2sig(ilat,ilon,:) = pdata(6,:);
        clear pdata
      end % ilat
    end % ilon

    %  save data here
    x = delete_file_ifexist(conc_file);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_mean'), Load_mean);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_median'), Load_median);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_Pos1sig'), Load_Pos1sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_Neg1sig'), Load_Neg1sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_Pos2sig'), Load_Pos2sig);
    f = write_netcdf (conc_file, 3, {'lat','lon','nseas'}, char('Load_Neg2sig'), Load_Neg2sig);
    f = write_netcdf (conc_file, 1, {'lat'}, char('lat'), clat);
    f = write_netcdf (conc_file, 1, {'lon'}, char('lon'), clon);
    f = write_netcdf (conc_file, 2, {'nseas','nmonth'}, char('month_season'), seas_num);

  % Add the attributes
    ncwriteatt(conc_file,char('Load_mean'),char('description'),'Mean Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_mean'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Load_median'),char('description'),'Median Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_median'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Load_Pos1sig'),char('description'),'+1 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_Pos1sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Load_Neg1sig'),char('description'),'-1 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_Neg1sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Load_Pos2sig'),char('description'),'+2 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_Pos2sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('Load_Neg2sig'),char('description'),'-2 Sigma Atmospheric Dust Loading');
    ncwriteatt(conc_file,char('Load_Neg2sig'),char('unit'),'g/m2');
    ncwriteatt(conc_file,char('lat'),char('description'),'Latitude');
    ncwriteatt(conc_file,char('lat'),char('unit'),'degrees_north');
    ncwriteatt(conc_file,char('lon'),char('description'),'Longitude');
    ncwriteatt(conc_file,char('lon'),char('unit'),'degrees_east');
    ncwriteatt(conc_file,'/','creation_date',datestr(now));

    ncwriteatt(conc_file,char('month_season'),char('description'),'Months in each season');

    clear Load_mean Load_median Load_Pos1sig Load_Neg1sig Load_Pos2sig Load_Neg2sig
  else
    disp([base_filename,'....already exist']);
  end % check if the file exist

  disp('  ');
