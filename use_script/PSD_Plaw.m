
function [ PSD_new,D_med_new ] = PSD_Plaw(PSD_model,D_median_model, D_new)
% function: It takes in a pointwise values of PSD and estimate the straight line joining the two points, in log space and equal to the number of desired sub-bin values.

PSD_new = null(1);
nD = size(D_median_model,2); % number of particle bins

%Loop over the radius of the median
  for i=2:nD
    tepo = D_new(find(D_new>=D_median_model(i-1) & D_new<=D_median_model(i)));
    a = 1+size(PSD_new,2);
    b = size(tepo,2)+size(PSD_new,2);
    A = exp(linspace(log(PSD_model(i-1)),log(PSD_model(i)),size(tepo,2)));
    PSD_new(a:b) = A;
    D_med_new(a:b) = tepo;
    clear tepo a b;
  end

end  % function
