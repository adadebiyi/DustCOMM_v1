function [dV_DdD_integral] = PSD_V_integral(PSD_lnV,D,D1,D2)

% Simple function to calculate the intergral between two diameter for dV/dD, Using the dV/dlnD.

  if (~isempty(find(D >=D1 & D<=D2))) % First check oif D1 and D2 is between the D_old
    i_bin_start = find(D>=D1,1,'first'); % greater than previous bin
    i_bin_end = find(D<=D2,1,'last');
    dV_DdD = (1./(10^(-6)*D(i_bin_start:i_bin_end))).*PSD_lnV(i_bin_start:i_bin_end);
    dV_DdD_integral = (10^(-6))*trapz(D(i_bin_start:i_bin_end),dV_DdD,2);
  else
    disp('D1 or D2 is outside the range of Ds provided');
    dV_DdD_integral = NaN;
  end

end
