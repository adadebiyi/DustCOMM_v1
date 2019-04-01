function [chi_square] =  chi_square_function_subbin (C,f_dust_model, Dlower, Dupper, D)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function minimize the chi parameter used to fit the generalized function in to the bias corrected dust mass fraction

% ========
% number of bins in the model
  nbin = max(size(f_dust_model)); % number of bins

  f_dust_model = reshape(f_dust_model,1,max(size(f_dust_model)));
  Dlower = reshape(Dlower,1,max(size(Dlower)));
  Dupper = reshape(Dupper,1,max(size(Dupper)));
  D = reshape(D,1,max(size(D)));

% the parameter for the function
  C1 = C(1);
  C2 = C(2);
  C3 = C(3);
  C4 = C(4);
  C5 = C(5);
  Cv_def = C(6);


% define the generalized function
  PSD = @(D) (1/Cv_def).*(1+erf(log(D./C1)/(sqrt(2)*log(C2)))).*(D.^C3).*exp(-(D/C4).^C5); %
  diff = [];

  % loop over each model bin
  for p=1:nbin %size(x,2) %cycling over data
    ibin = find(D>=Dlower(p) & D<Dupper(p));
    diff(p)=(log(trapz(D(ibin),PSD(D(ibin)),2))-log(f_dust_model(p)))^2;
  end %for, cycling over data
  chi_square = sqrt(sum(diff(:)));

end
