% Find critical T when the third virial coefficient change sign.
function result = criticalT_thirdB(Trange, coeffss)
addpath('..');
if (nargin<2)
  close all;
  xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5].';
  num_of_points = length(xis);
  coeffss=ones(num_of_points,1) * [1,2.5,4,1,1];
  coeffss(:, 5) = xis;
  Trange = [0.01, 0.5];
end
ncoeff = size(coeffss, 1);
Tcs = zeros(ncoeff, 1);
options = optimset('TolX',1e-4);
for rp=1:ncoeff
  coeffs = coeffss(rp,:);
  Tlow = Trange(1);
  Thigh = Trange(2);
  B3low = thirdB(Tlow, coeffs);
  B3high = thirdB(Thigh, coeffs);
  Tol = Thigh-Tlow;
  if(~ B3low<0 && B3high>0)
      warning('Zero point not found.');
      Tcs(rp) = nan;
      continue;
  end
  while(Tol>1e-4)
    Tcen = (Thigh+Tlow)/2;
    B3tmp = thirdB(Tcen, coeffs);
    if(B3tmp>0)
      Thigh = Tcen;
    else
      Tlow=Tcen;
    end
    Tol = Thigh-Tlow;
  end
  Tcs(rp) = Tcen;
end
if(nargin<2)
  plot(xis, Tcs.');
  %set(gca,'yscale','log');
end
if(nargout==1)
    result = Tcs;
end
end
