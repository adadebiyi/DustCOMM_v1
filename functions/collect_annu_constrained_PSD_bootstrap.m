function [fitPSD_bootstrap] = collect_annu_constrained_PSD_bootstrap (loc_lonlat,nboot,median_true)
% ==========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% ==========
% This function collects annually-averaged climatology of the constrained PSD
% for a particular location

% ==========
% Variables:
% fitPSD_bootstrap  -- the output
% loc_lonlat -- regions of interest. It could also be just one point
% nboot -- number of bootstrap to take
% median_true -- Specify if you want to take only the median

% =================
%  define some global variables
  define_global_variables (true);
  global  clat clon clev max_level_in_hPa ...
          dustcomm_dir_psd_annu fitPSD_filename_annu ...
          D_OBS
  
% Select only heights up to a threshold
  indlev = find(clev >= max_level_in_hPa);

% Iinformation for the spatial range
  if (all(loc_lonlat == 0))
    xlat = clat;
    xlon = clon;
  else
    if (size(loc_lonlat,2) == 1 || size(loc_lonlat,1) == 1)
      [d, indlat] = min(abs( clat-loc_lonlat(2)));
      [d, indlon] = min(abs( clon-loc_lonlat(1)));
      xlon = clon(indlon);
      xlat = clat(indlat);
    else
      % ; get the indices for the logitude and latitude  and height
      indlat = NaN(1,max(size(loc_lonlat)));
      indlon = indlat;
      for i=1:max(size(loc_lonlat))
        [d, indlat(i)] = min(abs( clat-loc_lonlat(i,1)));
        [d, indlon(i)] = min(abs( clon-loc_lonlat(i,2)));
        xlon = clon(indlon);
        xlat = clat(indlat);
      end
    end
  end %(all(loc_lonlat == 0) || recalculate)

  clear clat clon
% =================
%dimension
  nlev  = max(size(indlev)); %
  nlat  = max(size(xlat)); %
  nlon  = max(size(xlon)); %
  nD    = max(size(D_OBS)); %

% Store the data here...
  fitPSD_bootstrap = NaN(nD,nlev,nlon,nlat,nboot);

  for j=1:nboot
    for ilat=1:nlat
      for ilon=1:nlon
        clon_str = sprintf('%0.2f',xlon(ilon));
        clat_str = sprintf('%0.2f',xlat(ilat));
        j_str = sprintf('%05d',j);
        out_dir_j = strcat(dustcomm_dir_psd_annu,j_str,'/');
        conc_file = char(strcat(out_dir_j,fitPSD_filename_annu,'_',clon_str,'_',clat_str,'_',j_str,'.nc'));
        if (exist(conc_file, 'file'))
          fitPSD_bootstrap(:,:,ilon,ilat,j)  = ncread(conc_file,char(fitPSD_filename_annu));
        end
      end %ilon
    end % ilat
  end % j

%   Check if all data is there...
  if (any(~isnan(fitPSD_bootstrap(:))))
    if (median_true)
      for i = 1:nD
        for ilat=1:nlat
          for ilon=1:nlon
            for ilev=1:nlev
              A=sort(squeeze(fitPSD_bootstrap(i,ilev,ilon,ilat,:)));
              if (~all(isnan(A)))
                biCorr_Kokfit(i,ilev,ilon,ilat)= A(round(0.5*nboot));
              end
            end
          end
        end
      end
      fitPSD_bootstrap = biCorr_Kokfit;
      clear biCorr_Kokfit
    end
  else
    fitPSD_bootstrap = NaN(nD,nlev,nlon,nlat);
  end

    fitPSD_bootstrap = squeeze(fitPSD_bootstrap);


%   disp('Done...')
% end
