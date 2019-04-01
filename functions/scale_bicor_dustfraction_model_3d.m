function [model_column_loading,model_D_lower,model_D_upper,model_N] = scale_bicor_dustfraction_model_3d(load_column_mass,D_lower,D_upper,model_name,annu_or_seas)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function corrects the model dust mass fraction and put the size between 0.2 and 20 micron.

 % =====================
  %  define some global variables
  define_global_variables (true);
  global D_max D_min

 disp(' ');

  %  convert to character if not...
  annu_or_seas = char(annu_or_seas);

  % Get the reference model. This can be GISS or CNRM
%   evalc('[model_column_loading_orig,llon,llat,llev,ref_Dlower,ref_Dupper,N,no_data_sets] = get_model_dust_concentration_all_3d(annu_or_seas,''GISS'',false)');
     evalc('[model_column_loading_orig,llon,llat,llev,ref_Dlower,ref_Dupper,N,no_data_sets] = get_model_dust_concentration_data(annu_or_seas,''GISS'')');
     clear llon llat no_data_sets

  if (annu_or_seas == char('annu'))
      disp('Annual: Calculating additional dust mass fraction for the last bin..')
    % =====================
    % get all the model mass to calculate the one for GISS as the reference model.
      % convert to normalized fraction
      m=1;
      xdum = squeeze(model_column_loading_orig(m,1:N(m),:,:,:));
      sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),N(m),1,1,1));
      ref_fraction  = xdum./sum_xdum;

      clear xdum sum_xdum

%       bias correct the refeence model here
      evalc('[biCorr_ref_fraction] = bias_correction_all_3d_v1(ref_fraction,0,ref_Dlower,ref_Dupper,''GISS'',annu_or_seas);');

      % scale here...
      [model_column_loading,model_D_lower,model_D_upper,model_N] = scale_model_no_3d(load_column_mass,D_lower,D_upper,model_name,biCorr_ref_fraction,D_min,D_max);

  elseif (annu_or_seas == char('seas'))

    disp('Seasonal: Calculating additional dust mass fraction for las bin..')
    nseas = size(model_column_loading_orig,6);
    m=1;
%      obtain the annual average first...
      xdum = squeeze(mean(model_column_loading_orig(m,1:N(m),:,:,:,:),6,'omitnan'));
      sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),N(m),1,1,1));
      ref_fraction_ann = xdum./sum_xdum;

    for ix=1:nseas
      xdum = squeeze(model_column_loading_orig(m,1:N(m),:,:,:,ix));
      sum_xdum = squeeze(repelem(sum(xdum,1,'omitnan'),N(m),1,1,1));
      ref_fraction  = xdum./sum_xdum;

      clear sum_xdum xdum

      %       bias correct it here
      evalc('[biCorr_ref_fraction] = bias_correction_all_3d_v1(ref_fraction,ref_fraction_ann,ref_Dlower,ref_Dupper,''GISS'',annu_or_seas);');

      % scale here...
      [model_column_loading(:,:,:,:,ix),model_D_lower,model_D_upper,model_N] = scale_model_no_3d(load_column_mass(:,:,:,:,ix),D_lower,D_upper,model_name,biCorr_ref_fraction,D_min,D_max);
    end

  else
    error('Specify calculation for annually or seaonally-averaged variables')
  end

  clear model_column_loading_orig

  disp(' ');

end
