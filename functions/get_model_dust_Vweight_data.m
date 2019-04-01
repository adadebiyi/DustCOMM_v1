function [model_norm_profile] = get_model_dust_Vweight_data(annu_or_seas,model_names,savedata)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function gets the vertical weight of model dust

%% ========================
%  define some global variables
	define_global_variables (false);
  global out_dir_vwgt
%
% general names of seasons
  seas_str = {'DJF', 'MAM', 'JJA','SON'};
%
  nseas = max(size(seas_str));
  nmodel = max(size(model_names));
	tepo = cellstr(model_names);
	if all([max(size((num2str(cell2mat(tepo(1)))))) max(size((num2str(cell2mat(tepo(end))))))] == max(size('GISS')))
    nmodel = 1;
    model_names = {char(model_names)};
  end
%  convert char if not...
  annu_or_seas = char(annu_or_seas);

% ========================

  if (annu_or_seas == char('annu'))
   disp('get_model_dust_Vweight_all_3d: Annually-averaged dust vertical weights...');

%  Check if any of the output files already exist.
    exist_models = zeros(1,nmodel);
    for ix=1:nmodel
      conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_annu.nc'));
      if (exist(conc_file, 'file'))
        exist_models(ix) = exist_models(ix) + ix;
      else
        exist_models(ix) = exist_models(ix) + NaN;
      end
    end % ix
   check_all_file_exist = any(isnan(exist_models));

    if (check_all_file_exist || savedata)
      disp(' ');
      disp('get_model_dust_Vweight_all_3d: Getting the 3D dust concentration...');
%       evalc('[model_concentration,llon,llat,llev,~] = get_model_dust_concentration_all_3d(annu_or_seas,model_names,false)');
      evalc('[model_concentration,llon,llat,llev,~] = get_model_dust_concentration_data(annu_or_seas,model_names)');

      nlon = size(model_concentration,3);
      nlat = size(model_concentration,4);
      nlev = size(model_concentration,5);

      xmodel_concentration = squeeze(sum(model_concentration,2,'omitnan'));
      sum_hgt = repelem(sum(xmodel_concentration,4),1,1,1,nlev);
      model_norm_profile = xmodel_concentration./sum_hgt;

      clear model_concentration xmodel_concentration sum_hgt

      if (savedata)
        for ix=1:nmodel
          conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_annu.nc'));
          x = delete_file_ifexist(conc_file);
          f = write_netcdf (conc_file, 3, {'nlon','nlat','nlev'}, char('model_norm_profile'), squeeze(model_norm_profile(ix,:,:,:)));
          f = write_netcdf (conc_file, 1, {'nlon'}, char('lon'), llon);
          f = write_netcdf (conc_file, 1, {'nlat'}, char('lat'), llat);
          f = write_netcdf (conc_file, 1, {'nlev'}, char('lev'), llev);
        end
      end
    else
%       all data is available; Collect here
      disp('Data are available...');
      model_norm_profile = [];
      for ix=1:nmodel
%         disp(['Collecting data for= ',model_names(ix)]);
        conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_annu.nc'));
        model_norm_profile(ix,:,:,:) = ncread(conc_file,'model_norm_profile');
      end
    end

% ========================
  elseif (annu_or_seas == char('seas'))
    disp('get_model_dust_Vweight_all_3d: Seasonally-averaged dust vertical weights...');

    %  Check if any of the output files already exist.
    exist_models = zeros(nmodel,max(size(seas_str)));
    for ix=1:nmodel
      for k=1:max(size(seas_str))
        conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_',seas_str(k),'.nc'));
        if (exist(conc_file, 'file'))
          exist_models(ix,k) = exist_models(ix) + ix;
        else
          exist_models(ix,k) = exist_models(ix) + NaN;
        end
      end
    end

    check_all_file_exist = any(isnan(exist_models(:)));

    if (check_all_file_exist || savedata)
      disp(' ');
      disp('get_model_dust_Vweight_all_3d: Getting the 3D dust concentration...');
%       evalc('[model_concentration,llon,llat,llev,~] = get_model_dust_concentration_all_3d(annu_or_seas,model_names,false)');
      evalc('[model_concentration,llon,llat,llev,~] = get_model_dust_concentration_data(annu_or_seas,model_names)');
      nlon = size(model_concentration,3);
      nlat = size(model_concentration,4);
      nlev = size(model_concentration,5);

      xmodel_concentration = squeeze(sum(model_concentration,2,'omitnan'));
      sum_hgt = repelem(sum(xmodel_concentration,4),1,1,1,nlev,1);
      model_norm_profile = xmodel_concentration./sum_hgt;

      clear model_concentration xmodel_concentration sum_hgt

      if (savedata)
        for ix=1:nmodel
          for k=1:nseas
            conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_',seas_str(k),'.nc'));
            x = delete_file_ifexist(conc_file);
            f = write_netcdf (conc_file, 3, {'nlon','nlat','nlev'}, char('model_norm_profile'), squeeze(model_norm_profile(ix,:,:,:,k)));
            f = write_netcdf (conc_file, 1, {'nlon'}, char('lon'), llon);
            f = write_netcdf (conc_file, 1, {'nlat'}, char('lat'), llat);
            f = write_netcdf (conc_file, 1, {'nlev'}, char('lev'), llev);
          end
        end
      end
    else
%       all data is available; Collect here
      disp('Data are available...');
      model_norm_profile = [];
      for ix=1:nmodel
        for k=1:nseas
          disp(['Collecting data for= ',model_names(ix),' for ',seas_str(k)]);
          conc_file = char(strcat(out_dir_vwgt,model_names(ix),'_3D_model_dust_Vweight_',seas_str(k),'.nc'));
          model_norm_profile(ix,:,:,:,k) = ncread(conc_file,'model_norm_profile');
        end
      end
    end

  end

%   model_norm_profile = squeeze(model_norm_profile);
%   Done...

% end
