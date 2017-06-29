%Calculate third virial coefficient

function result = thirdB(Ts, coeffs, Tol)
addpath('..');
if (nargin<3)
    Tol = 1e-5;
end
if (nargin<2)
  close all;
  coeffs=[1,2.5,4,1,1];
  Ts=linspace(0.29, 0.31, 10);
end
betas = 1./Ts;

sigma = coeffs(1);
lambda = coeffs(2);
kappa = coeffs(3);
epsilon = coeffs(4);
xi = coeffs(5);
centrr = lambda*sigma;
upperr = kappa*sigma;

nbeta = size(betas, 2);
B3s = zeros(size(betas));
Intfun = @(r1, r2, bcoeffs)((boltzmanweight(r1,bcoeffs)-1).* ...
  (boltzmanweight(r2,bcoeffs)-1).*(boltzmanweight(r2-r1,bcoeffs)-1));
Intfun_special = @(r1, r2, bcoeffs) ...
  ((1-boltzmanweight(r2,bcoeffs)).*(boltzmanweight(r2-r1,bcoeffs)-1));
for rp=1:nbeta
  beta = betas(rp);
  w1 = -1;
  w2 = exp(beta*epsilon)-1;
  w3 = @(r)(exp(-beta*potential(r, coeffs))-1);
  i11 = -sigma^2;
  i12 = integral2(@(r1, r2)(w1*w2.*w3(r2-r1)), 0, sigma, sigma, centrr,'AbsTol',1e-2,'RelTol',Tol,'method','iterated');
  i13 = integral2(@(r1, r2)(w1*w3(r2).*w3(r2-r1)), 0, sigma, centrr, upperr,'AbsTol',1e-2,'RelTol',Tol,'method','iterated');
  i22 = integral2(@(r1, r2)(w2*w2.*w3(r2-r1)), sigma, centrr, sigma, centrr,'AbsTol',1e-2,'RelTol',Tol,'method','iterated');
  i23 = integral2(@(r1, r2)(w2.*w3(r2).*w3(r2-r1)), sigma, centrr, centrr, upperr,'AbsTol',1e-2,'RelTol',Tol,'method','iterated');
  i33 = integral2(@(r1, r2)(w3(r1).*w3(r2).*w3(r2-r1)), centrr, upperr, centrr, upperr,'AbsTol',1e-2,'RelTol',Tol,'method','iterated');
  
  B3s(rp) = -(i11+i22+i33+2*(i12+i13+i23));
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