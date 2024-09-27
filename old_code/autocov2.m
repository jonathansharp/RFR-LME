function [acov,acor,lag,dof] = autocov2 (t,y,max_lag)
  
%   Computes Autocovariance of a residual time-series
%   accounting for missing data values. Note: no
%   trend is removed. This is so that this function
%   can be used with colored noise models that 
%   may have small trends. Removing the small
%   trend will alter the covariance.
%   
%   t is time, and assumed to be full. max_lag is the 
%   the maximum lag to compute, in same units as t (i.e., years, months, days, etc)
%   Missing data values must be defined by %nan.
%   
%   This version does not interpolate to get full matrices, but
%   uses loops to edit out missing data
%   
%   Bin size is determined by the time difference between first 2 points.
%   
%   DOF is computed by decorrelation time where correlation goes below 0.5
%  
  n    = size(t);
  nmax = max(n);
  del  = mode(diff(t));
  
  nlags = round(max_lag/del);
  
  epoch = t(1);
  dt = t - epoch;
  
  %check for %nan values 
  
  a = isnan(y);
  
  j = 1;
  for i = 1:nmax
    if (~a(i))
      yt(j,:) = y(i);
      tt(j,:) = dt(i);
      j = j+1;
    end
  end   
  
  n = size(yt);
  nmax2 = max(n);  
    
  acov = zeros(nlags,1);
  ncov = zeros(nlags,1);
  nl = nlags-1;
  
  for i = 1:nmax2
   for j = i-nl:i+nl
     if (j <= nmax2 & j >= 1)
       k = round(abs((tt(j)-tt(i))/del))+1;
       if (k >=1 & k <= nlags)
       acov(k,1) = acov(k,1) + yt(i)*yt(j);
       ncov(k,1) = ncov(k,1) + 1.0;
       end
     end
   end
 end
 
 acov(1) = var(yt);
 for i = 2:nlags
   acov(i,:) = acov(i)/(ncov(i)-1);
 end
 lag = 0:nlags-1;
 lag = lag*del;
 lag = lag';
 acor = acov/acov(1);
   
 % compute effective degrees of freedom based on 2x time for acor to first drop 
 % below 1/e
 
 test = 0.0;
 efold = 1/exp(1);
 decor = lag(5);
 for i=2:nlags
     if acor(i) <= efold && test ~= 123.45
         decor = lag(i)*2.0;
         test = 123.45;
     end
 end
 
 ncor = round(decor/del);
 dof = round(nmax2/ncor);
