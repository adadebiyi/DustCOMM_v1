function [model_concentration,llon,llat,llev,model_D_lower,model_D_upper,N,no_data_sets] = get_model_dust_concentration_data(annu_or_seas,model_names)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function loads the 3D model dust concentration (ug/m3) for
% annually-averaged and seasonally average values.

% %% ========================

%  define some global variables
  define_global_variables (false);
  global out_dir_conc

%
% general names of seasons
  seas_str = {'DJF', 'MAM', 'JJA','SON'};
%
% Dimensions
  nseas = max(size(seas_str));
  nmodel = max(size(model_names));
  tepo = cellstr(model_names);
	if all([max(size((num2str(cell2mat(tepo(1)))))) max(size((num2str(cell2mat(tepo(end))))))] == max(size('GISS')))
    nmodel = 1;
    model_names = {char(model_names)};
  end

%  convert char if not...
  annu_or_seas = char(annu_or_seas);

%% ========================
  if (annu_or_seas == char('annu'))
   disp('get_model_dust_concentration_all_3d: annual Concentration...');

%  Check if any of the output files already exist.
    exist_models = zeros(1,nmodel);
    for ix=1:nmodel
      conc_file = char(strcat(out_dir_conc,model_names(ix),'_3D_model_dust_concentration_annu.nc'));
      if (exist(conc_file, 'file'))
        exist_models(ix) = exist_models(ix) + ix;
      else
        exist_models(ix) = exist_models(ix) + NaN;
      end
    end % ix
   check_all_file_exist = any(isnan(exist_models));

% If the data already exist..dont calculate again
    if (check_all_file_exist)
      disp('get_model_dust_concentration_all_3d: Provide annually-averaged dust concentration data...');
    else

      disp('get_model_dust_concentration_all_3d: File data already exist....');
      disp('get_model_dust_concentration_all_3d: Result will be copied from the data...');
      disp(' ');

      %Copy from the stored file here.
      for m=1:nmodel
        %     store the data
        conc_file = char(strcat(out_dir_conc,model_names(m),'_3D_model_dust_concentration_annu.nc'));
%         disp(['get_model_dust_concentration_all_3d: Collecting data for= ',model_names(m)]);
        if (m == 1)
          no_data_sets = ncread(conc_file,'no_data_sets');
          llon = ncread(conc_file,'lon');
          llat = ncread(conc_file,'lat');
          llev = ncread(conc_file,'lev');
        end

        N(m) = ncread(conc_file,'N');
        model_concentration(m,1:N(m),:,:,:) = ncread(conc_file,'dust_concentration');
        model_D_lower(m,1:N(m)) = ncread(conc_file,'D_lower');
        model_D_upper(m,1:N(m)) = ncread(conc_file,'D_upper');
      end % end do

    end  % check_all_file_exist

  elseif (annu_or_seas == char('seas'))

    disp('get_model_dust_concentration_all_3d: Seasonal Concentration');


    %  Check if any of the output files alreadt exist.
    exist_models = zeros(nmodel,max(size(seas_str)));
    for ix=1:nmodel
      for k=1:max(size(seas_str))
        conc_file = char(strcat(out_dir_conc,model_names(ix),'_3D_model_dust_concentration_',seas_str(k),'.nc'));
        if (exist(conc_file, 'file'))
%                   disp([model_names(m), seas_str(k)]);
          exist_models(ix,k) = exist_models(ix) + ix;
        else
          exist_models(ix,k) = exist_models(ix) + NaN;
        end
      end
    end

    check_all_file_exist = any(isnan(exist_models(:)));

% If the data already exist..dont calculate again
    if (check_all_file_exist)
      disp('get_model_dust_concentration_all_3d: Provide seasonally-averaged dust concentration data...');
    else

      disp('get_model_dust_concentration_all_3d: File data already exist....');
      disp('get_model_dust_concentration_all_3d: Result will be copied from the data...');
      disp(' ');

      for m=1:nmodel
        for k=1:max(size(seas_str))

          conc_file = char(strcat(out_dir_conc,model_names(m),'_3D_model_dust_concentration_',seas_str(k),'.nc'));
%           disp(['get_model_dust_concentration_all_3d: Collecting data for= ',model_names(m),' for ',seas_str(k)]);

          if (k == 1)
            no_data_sets = ncread(conc_file,'no_data_sets');
            llon = ncread(conc_file,'lon');
            llat = ncread(conc_file,'lat');
            llev = ncread(conc_file,'lev');
          end
          N(m) = ncread(conc_file,'N');
          model_D_lower(m,1:N(m)) = ncread(conc_file,'D_lower');
          model_D_upper(m,1:N(m)) = ncread(conc_file,'D_upper');
          model_concentration(m,1:N(m),:,:,:,k) = ncread(conc_file,'dust_concentration');
        end
      end

    end  %check_all_file_exist

  end %annu_or_seas

% end
