% Calculate fourth virial coefficient
% Did not divide the integration region manually
% Unsuccessful attempt
function result = fourthB(Ts, coeffs)
addpath('..');
if (nargin<2)
  close all;
  coeffs=[1,2.5,2.5,1,0];
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
B4s = zeros(size(betas));

mfun = @(r, bcoeffs)(boltzmanweight(r,bcoeffs)-1); % mayer function

for rp=1:nbeta
  bcoeffs = coeffs;
  bcoeffs(4) = bcoeffs(4)*betas(rp);

  B4A = integral3(@(r1, r2, r3)(mfun(r1, bcoeffs).*mfun(r2, bcoeffs).*mfun(r3, bcoeffs) ...
      .*mfun(r2-r1, bcoeffs).*mfun(r3-r2, bcoeffs).*mfun(r3-r1, bcoeffs)), ...
      0, upperr, 0, upperr, 0, upperr, ...
      'AbsTol',1e-2,'RelTol',1e-5,'method','iterated');
  B4B = integral3(@(r1, r2, r3)(mfun(r1, bcoeffs).*mfun(r2, bcoeffs).*mfun(r3, bcoeffs) ...
      .*mfun(r3-r1, bcoeffs).*mfun(r3-r2, bcoeffs)), ...
      0, upperr, 0, upperr, 0, upperr, ...
      'AbsTol',1e-2,'RelTol',1e-5,'method','iterated');
  B4C = integral3(@(r1, r2, r3)(mfun(r1, bcoeffs).*mfun(r2, bcoeffs) ...
      .*mfun(r3-r1, bcoeffs).*mfun(r3-r2, bcoeffs)), ...
      0, upperr, 0, upperr, 0, 2*upperr, ...
      'AbsTol',1e-2,'RelTol',1e-5,'method','iterated');
  
  B4s(rp) = -0.5*(B4A+6*B4B+3*B4C);   % 4*(-1/8)*(B4A+6*B4B+3*B4C)
end
if(nargin<2)
  plot(Ts, alignvalue(B4s));
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
