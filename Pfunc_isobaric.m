% Calculate the partition function (\Zeta/N)
% P - pressure, beta - temperature,
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
% which defines SALR potential u = \begin{case}
% \infty, (r < \sigma) \\
% -\epsilon, (\sigma <= r < \lambda*\sigma) \\
% \zeta\epsilon(\kappa-r/\sigma), (\lambda*\sigma <= r < \kappa*\sigma) \\
% 0, (r > \kappa*\sigma)
% set 2*pi*m/h^2 = 1
% consider next nearest neighbor
function zeta = Pfunc_isobaric(P, beta, coeffs, divide)
    % set 2*pi*m/h^2 = 1
    if ( nargin<4)
        divide = 100;
    end
    sigma = coeffs(1);
    kappasigma = coeffs(3)*sigma;
    
    ds = (kappasigma-sigma)/(divide-1);
    rlist = linspace(sigma, kappasigma, divide).'+ds/2; % Magic. Increased the accurancy pretty much
    pot = potential(rlist, coeffs);
    pot = pot + P*(rlist-sigma);
    tmat = pot*ones(1,divide);
    for rp=1:divide
        % plus next nearest neighbor
        tmat(:,rp) = exp(-beta*(tmat(:,rp) + potential(rlist+rlist(rp), ...
            coeffs)));
    end
    bpks = beta*P*(coeffs(3)-1)*sigma;  % contribution of tail
    tail = exp(-bpks)/(beta*P);
    tmat(divide, 1:divide) = tail/ds; % tail
    lmax = eigs(tmat,1);
    zeta = sqrt(1/beta)*lmax*ds*exp(-beta*P*sigma);
end
