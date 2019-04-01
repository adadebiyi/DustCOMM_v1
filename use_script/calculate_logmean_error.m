function  [fit_dust,bias] = calculate_logmean_error(fit_PSD,D_OBS,f_dust_bicor,Dlower,Dupper)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function calculates the log mean error..

% First check if the size matches between variables.
% make it all row data..
% convert to data row array
  fit_PSD = reshape(fit_PSD,1,max(size(fit_PSD)));
  D_OBS = reshape(D_OBS,1,max(size(D_OBS)));
  f_dust_bicor = reshape(f_dust_bicor,1,max(size(f_dust_bicor)));
  Dlower = reshape(Dlower,1,max(size(Dlower)));
  Dupper = reshape(Dupper,1,max(size(Dupper)));

  if (max(size(fit_PSD)) ~= max(size(D_OBS)))
    erro('fit_PSD and D_OBS must be the same size')
  end
  % if (max(size(f_dust_bicor)) ~= max(size(Dlower)) || max(size(f_dust_bicor)) ~= max(size(Dupper)) )
  %   erro('f_dust_bicor Dlower and Dupper are not the same size')
  % end

  nbin = max(size(f_dust_bicor));
  fit_dust = zeros(1,nbin);

  for p=1:nbin
    ibin = find(D_OBS>=Dlower(p) & D_OBS<Dupper(p));
    fit_dust(p)=trapz(D_OBS(ibin),fit_PSD(ibin),2);
  end %for, cycling over data

  fit_dust = fit_dust/sum(fit_dust);
  bias = sum(abs(log10(fit_dust./f_dust_bicor)))*(1./nbin);
end
