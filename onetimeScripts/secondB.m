%Calculate second virial coefficient

function result = secondB(Ts, coeffs)
if (nargin<2)
  close all;
  coeffs=[1,2.5,4,1,100];
  Ts=linspace(0.1, 1);
end
betas = 1./Ts;

sigma = coeffs(1);
lambda = coeffs(2);
kappa = coeffs(3);
epsilon = coeffs(4);
xi = coeffs(5);

nbeta = size(betas, 2);
B2s = zeros(size(betas));
for rp=1:nbeta
    beta = betas(rp);
    if(xi==0)
      B2s(rp) = -sigma*(exp(beta*epsilon)*(lambda-1)-lambda);
    else
      B2s(rp) = -sigma*((exp(beta*epsilon)*(lambda-1)-lambda)+ ...
        (1-exp(-beta*xi*epsilon*(kappa-lambda)))/(beta*xi*epsilon)+lambda-kappa);
    end
end
if(nargin<2)
  plot(Ts, B2s);
  set(gca,'yscale','log');
end
if(nargout==1)
    result = B2s;
end
end