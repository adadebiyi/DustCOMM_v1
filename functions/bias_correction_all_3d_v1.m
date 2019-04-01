function [biCorr_dust_fraction] = bias_correction_all_3d_v1(unCorr_dust_fraction,unCorr_dust_fraction_ann,D_lower,D_upper,model_name,annu_or_seas)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This code bias correct both annually-averaged and seasonally-averaged
% model dust distribtion.
%It uses the constrained globally-averaged PSD to bias correct the
%spatially-varying annully-averaged distribution

% % =====================
%  define some global variables
  define_global_variables (true);
  global gPSD_OBS_median D_OBS

%  convert char if not...
  annu_or_seas = char(annu_or_seas);

if (unCorr_dust_fraction_ann == 0 & annu_or_seas == char('annu'))
  %then the calculation is for anually-averaged

  disp('bias_correction_all_3d: Applying bias correction to the Annually-averaged dust mass fraction...');
%Dimensions sizes
    nbin = size(unCorr_dust_fraction,1);
    nlon = size(unCorr_dust_fraction,2);
    nlat = size(unCorr_dust_fraction,3);
    nlev = size(unCorr_dust_fraction,4);

  %Calculate the globally-averaged model PSD here
%     evalc('[model_norm_profile] = get_model_dust_Vweight_all_3d(annu_or_seas,model_name,false)');
    evalc('[model_norm_profile] = get_model_dust_Vweight_data(annu_or_seas,model_name,false)');
    model_norm_profile = repelem(model_norm_profile,nbin,1,1,1);

    % Globally-averaged annually-averaged model PSD
    PSD_clim = squeeze(mean(mean(sum(unCorr_dust_fraction.*model_norm_profile,4,'omitnan'),2,'omitnan'),3,'omitnan'));
    PSD_ann(1,:) = PSD_clim;
    PSD_ann = PSD_ann./sum(PSD_ann);

%     integrate the constrained globally-averaged PSD to the model sizes
    dV_dD_integral = zeros(1,nbin);
    for k=1:nbin
      i_bin_start = find(D_OBS>=D_lower(k),1,'first');
      i_bin_end = find(D_OBS<D_upper(k),1,'last');
      dV_DdD = (1./(10^(-6)*D_OBS(i_bin_start:i_bin_end))).*gPSD_OBS_median(i_bin_start:i_bin_end);
      dV_dD_integral(k) = (10^(-6))*trapz(D_OBS(i_bin_start:i_bin_end),dV_DdD,2);
      clear dV_DdD
    end %k

    % normalize the integral.
    dV_dD_integral = dV_dD_integral/sum(dV_dD_integral);

% calculate the the correcton factor for the specific model
    alpha_constant = dV_dD_integral./PSD_ann;
    alpha = repelem(alpha_constant',1,nlon,nlat,nlev); % repeat for each location
%correct the model dust mass faction
    biCorr_dust_fraction = squeeze(unCorr_dust_fraction(1:nbin,:,:,:)).*alpha;

% %normalize
% sum_column_dust_dist = repelem(sum(biCorr_dust_fraction,1),nbin,1,1);
% biCorr_dust_fraction = biCorr_dust_fraction./sum_column_dust_dist;

elseif (any(unCorr_dust_fraction_ann(:) ~= 0) & annu_or_seas == char('seas'))
  %then the calculation is for seasonally-averaged
    disp('bias_correction_all_3d: Applying bias correction to the Seasonally-averaged dust mass fraction...');
  %Dimensions sizes
    nbin = size(unCorr_dust_fraction,1);
    nlon = size(unCorr_dust_fraction,2);
    nlat = size(unCorr_dust_fraction,3);
    nlev = size(unCorr_dust_fraction,4);
    if (length(size(unCorr_dust_fraction)) >= 5)
      nseas= size(unCorr_dust_fraction,5);
    else
      nseas = 1;
    end

    %Calculate the equivalent bias-corrected annually-averaged distribution
%     [biCorr_dust_fraction_ann] = bias_correction_all_3d(unCorr_dust_fraction_ann,0,D_lower,D_upper,'annu');
    evalc('[biCorr_dust_fraction_ann] = bias_correction_all_3d_v1(unCorr_dust_fraction_ann,0,D_lower,D_upper,model_name,''annu'');');

% where to store the seasonally-averaged bias-correction
%     Corr_dust_fraction = [];
% now calculate for the season
    if (nseas == 1)
        biasCorr_column_loading_seas = unCorr_dust_fraction.*(biCorr_dust_fraction_ann./unCorr_dust_fraction_ann);
        biCorr_dust_fraction(1:nbin,:,:,:) = biasCorr_column_loading_seas./repelem(nansum(nanmean(nanmean(biasCorr_column_loading_seas,2),3),1),nbin,nlon,nlat,1);
    else
      for k=1:nseas
        unCorr_column_loading_seas = squeeze(unCorr_dust_fraction(1:nbin,:,:,:,k));
        biasCorr_column_loading_seas = unCorr_column_loading_seas.*(biCorr_dust_fraction_ann./unCorr_dust_fraction_ann);
        % biCorr_dust_fraction(1:nbin,:,:,:,k) = biasCorr_column_loading_seas./repelem(nansum(nanmean(nanmean(biasCorr_column_loading_seas,2),3),1),nbin,nlon,nlat,1);
        biCorr_dust_fraction(1:nbin,:,:,:,k) = biasCorr_column_loading_seas;
      end
    end


else
  error('Seasonal Calculation is requested by but the annually-averaged values are not available.');
  biCorr_dust_fraction = 0;
end

end
