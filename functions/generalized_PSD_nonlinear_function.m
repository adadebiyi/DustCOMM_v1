function [fit_PSD,param,range_warning,f_dust_model_out] =  generalized_PSD_nonlinear_function(f_dust_model, Dlower, Dupper, D,init_param,param_bounds,fit_error_threshold,refit)

% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% This function fits the the generalized function in the model dust
% faction based on a paramter bound.


% =====================
% convert to data row array
  f_dust_model = reshape(f_dust_model,1,max(size(f_dust_model)));
  Dlower = reshape(Dlower,1,max(size(Dlower)));
  Dupper = reshape(Dupper,1,max(size(Dupper)));
  D = reshape(D,1,max(size(D)));

  if (param_bounds == 0)
% define the range of parameters allowed
    param_bounds(1,1) = 0.25; %min_Ds
    param_bounds(1,2) = 6; %max_Ds
    param_bounds(2,1) = 1.6; %min_sigma_s
    param_bounds(2,2) = 4; %max_sigma_s
    param_bounds(3,1) = -10; %min_b
    param_bounds(3,2) = 4; %max_b
    param_bounds(4,1) = 1; %min_lamda
    param_bounds(4,2) = 30; %max_lamda
    param_bounds(5,1) = 1; %min_alpha
    param_bounds(5,2) = 6; %max_alpha
    param_bounds(6,1) = 0; %min_CV
    param_bounds(6,2) = 200; %max_CV
  end

% define the initial range for the parameter to fit the data...
  if (init_param == 0)
%     set the initial parameters to the globally-averaged values.
    init_param(1) = 3.7;
    init_param(2) = 2.2;
    init_param(3) = -2;
    init_param(4) = 13;
    init_param(5) = 3;
    init_param(6) = 1.;
  end

% % =====================
%     First Calculate the fit using the 5-param function with alpha
%   Fit the generalized function to the model fraction to calculate the subbin distribution
% % =====================

% The function
% PSD = @(D) (1/CV)*(1+erf(log(D./D_s)/(sqrt(2)*log(sigma_s)))).*(D.^b).*exp(-(D/lamda).^alpha);
%   define the function handle for the chi- square
  gen_kok_fun = @(init_param) chi_square_function_subbin (init_param,f_dust_model, Dlower, Dupper, D);

% Options to consider
  options = optimset('MaxFunEvals',2000,'TolFun',1e-12);
  [param,chi_square] = fminsearchbnd(gen_kok_fun,init_param,param_bounds(:,1),param_bounds(:,2),options);

  PSD = @(D) (1/param(6))*(1+erf(log(D./param(1))/(sqrt(2)*log(param(2))))).*(D.^param(3)).*exp(-(D/param(4)).^param(5));
  param(6) = param(6)*trapz(D,PSD(D)); PSD = @(D) (1/param(6))*(1+erf(log(D./param(1))/(sqrt(2)*log(param(2))))).*(D.^param(3)).*exp(-(D/param(4)).^param(5));
  fit_PSD = PSD(D);

  %   calculate the error
  [C,PM20error] = calculate_logmean_error(fit_PSD,D,f_dust_model,Dlower,Dupper);
  
  f_dust_model_out = f_dust_model;

  if (refit)
  % % =====================
  % 	Check if the error is more than the thresholds
  % % =====================
  % Check if the error is within the allowed limit.
   if (PM20error > fit_error_threshold)
     disp('generalized_PSD_nonlinear_function: Error is too high...Trying to refit the generalized equation...' )
      %    remove the last bin...This is the bin that was scaled in to
      %    account for the full distribution
      xf_dust_model = f_dust_model(1:end-1);
      xDlower = Dlower(1:end-1);
      xDupper = Dupper(1:end-1);

      nbin = max(size(xf_dust_model));

  %     check if there is enough point for power law
  %     This signifies that the maximum bin in the distrubtion has few
  %     fraction that the maximum
      if ((nbin - find(xf_dust_model == max(xf_dust_model)) ~= 0))
        % There are atleast two points for power law dustribution
        % choose the last two points
        Dl = Dlower(end-2:end-1);
        Du = Dupper(end-2:end-1);
        Dg = sqrt(Dl.*Du);
        twobin = xf_dust_model(end-1:end);

  % ; get the subbin diameter between the two points and calculate the Plaw
        D_new = D(D >= Dg(1) & D <= Dg(2));
        [PSD_new,Dx ] = PSD_Plaw(twobin,Dg, D_new);
  % extend the Plaw to 20microns
        p = polyfit(log(D_new),log(PSD_new),1);
        D_new = D(D >= Dg(1) & D <= Dupper(end));
        Ap = exp(polyval(p,log(D_new)));

  % identify the best value for the last bin.
        Dlast = sqrt(Du(2)*Dupper(end));
        [C, zbin] = min(abs(D_new-Dlast));
        lastbin = Ap(zbin);

% add the data on the dV/dlnD space and convert back to dV/dD
        f_dust_model_new = [xf_dust_model./(log(xDupper) - log(xDlower)),lastbin].*(log(Dupper) - log(Dlower));
% add the data on the dV/dD space
%         f_dust_model_new = [xf_dust_model,lastbin];
% Normalize
        f_dust_model_new = f_dust_model_new/sum(f_dust_model_new);

        %   define the function handle for the chi- square
        gen_kok_fun = @(init_param) chi_square_function_subbin (init_param,f_dust_model_new, Dlower, Dupper, D);
        [aparam,chi_square] = fminsearchbnd(gen_kok_fun,init_param,param_bounds(:,1),param_bounds(:,2),options);

        PSD = @(D) (1/aparam(6))*(1+erf(log(D./aparam(1))/(sqrt(2)*log(aparam(2))))).*(D.^aparam(3)).*exp(-(D/aparam(4)).^aparam(5));
        aparam(6) = aparam(6)*trapz(D,PSD(D)); PSD = @(D) (1/aparam(6))*(1+erf(log(D./aparam(1))/(sqrt(2)*log(aparam(2))))).*(D.^aparam(3)).*exp(-(D/aparam(4)).^aparam(5));
        afit_PSD_kok = PSD(D);
        %   calculate the error
        [C,xerrorPM20] = calculate_logmean_error(afit_PSD_kok,D,f_dust_model_new,Dlower,Dupper);

        if (xerrorPM20 <= PM20error)
%         if (xerrorPM20 <= fit_error_threshold)
          fit_PSD = afit_PSD_kok;
          param = aparam;
          f_dust_model_out = f_dust_model_new;
        end % check if xerrorPM20 is good
      end  % if ther is enough points

   else
  %    if less...do nothing
  %       accept the original as final.
      disp('error is within limit... No need for additional calculation.')
   end  % Check error threshold

  end % refit

  %   ; now check if the parameter is within the limit, and flag if it is at the edge.
  range_warning = zeros(1,max(size(param)));
  if (single(param(1)) == param_bounds(1,1)); range_warning(1) = 1; elseif (single(param(1)) == param_bounds(1,2)); range_warning(1) = 2; end
  if (single(param(2)) == param_bounds(2,1)); range_warning(2) =  1; elseif (single(param(2)) == param_bounds(2,2)); range_warning(2) = 2; end
  if (single(param(3)) == param_bounds(3,1)); range_warning(3) = 1; elseif (single(param(3)) == param_bounds(3,2)); range_warning(3) = 2; end
  if (single(param(4)) == param_bounds(4,1)); range_warning(4) =  1; elseif (single(param(4)) == param_bounds(4,2)); range_warning(4) = 2; end
  if (single(param(5)) == param_bounds(5,1)); range_warning(5) =  1; elseif (single(param(5)) == param_bounds(5,2)); range_warning(5) = 2; end
  if (single(param(6)) == param_bounds(6,1)); range_warning(6) =  1; elseif (single(param(6)) == param_bounds(6,2)); range_warning(6) = 2; end

% end
