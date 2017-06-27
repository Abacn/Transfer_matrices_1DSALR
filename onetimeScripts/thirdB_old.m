% Calculate third virial coefficient
% Did not divide the integration region manually
function result = thirdB_old(Ts, coeffs)
addpath('..');
if (nargin<2)
  close all;
  coeffs=[1,2.5,4,1,1];
  Ts=linspace(0.25, 0.35, 10);
end
betas = 1./Ts;

sigma = coeffs(1);
lambda = coeffs(2);
kappa = coeffs(3);
epsilon = coeffs(4);
xi = coeffs(5);
upperr = kappa*sigma;

nbeta = size(betas, 2);
B3s = zeros(size(betas));
Intfun = @(r1, r2, bcoeffs)((boltzmanweight(r1,bcoeffs)-1).* ...
  (boltzmanweight(r2,bcoeffs)-1).*(boltzmanweight(r2-r1,bcoeffs)-1));
for rp=1:nbeta
  bcoeffs = coeffs;
  bcoeffs(4) = bcoeffs(4)*betas(rp);
  tmp = integral2(@(r1, r2)(Intfun(r1, r2, bcoeffs)), 0, upperr, 0, upperr);
  B3s(rp) = -tmp;
end
if(nargin<2)
  plot(Ts, alignvalue(B3s));
  %set(gca,'yscale','log');
end
if(nargout==1)
    result = B3s;
end
end

function out = alignvalue(in)
    out = zeros(size(in));
    bools1 = in>10;
    out(bools1) = log10(in(bools1));
    bools2 = in<-10;
    out(bools2) = -log10(-in(bools2));
    bools3 = ~(bools1 | bools2);
    out(bools3) = in(bools3)/10;
end
