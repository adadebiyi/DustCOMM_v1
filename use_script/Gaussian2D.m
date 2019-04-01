function [value] = Gaussian2D(mu, sigma)

  nmu = size(squeeze(mu));
  nsigma = size(squeeze(sigma));

  % In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1).
  % r = a + (b-a).*rand(N,1).
  % randomly generate numbers between -1:1

  % check if they are both the same
    if (isequal(nmu,nsigma))
      % randomly generate numbers between -1:1
      rand_no = -1 + 2.*rand(nmu(1),nmu(2));
      value = mu + sigma.*rand_no;
    else
      disp('The dimension of mu is not the same as sigma')
    end
