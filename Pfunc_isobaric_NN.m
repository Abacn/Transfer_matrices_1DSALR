% Calculate the partition function (\Zeta/N), Nearest neigbor
% P - pressure, beta - temperature,
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
% which defines SALR potential u = \begin{case}
% \infty, (r < \sigma) \\
% -\epsilon, (\sigma <= r < \lambda*\sigma) \\
% \zeta\epsilon(\kappa-r/\sigma), (\lambda*\sigma <= r < \kappa*\sigma) \\
% 0, (r > \kappa*\sigma)
% set 2*pi*m/h^2 = 1
function pfunc = Pfunc_isobaric_NN(P, beta, coeffs, divide)
    % set 2*pi*m/h^2 = 1
    sigma=coeffs(1);
    lambda=coeffs(2);
    kappa=coeffs(3);
    epsilon=coeffs(4);
    zeta = coeffs(5);
    fenmu = beta*(P-zeta*epsilon/sigma);
    if abs(fenmu)<1e-10
        % The analytical form suffers with 0/0
      pfunc = (exp(beta*epsilon)*(exp(-beta*P*sigma)-exp(-beta*P*lambda* ...
      sigma))+exp(-beta*P*kappa*sigma))/(beta*P) + ...
      (exp(-beta*kappa*P*sigma)*kappa-exp(-beta*lambda*P*sigma+beta*epsilon ...
      *(-kappa+lambda)*zeta)*lambda)*sigma;
    else
      pfunc = (exp(beta*epsilon)*(exp(-beta*P*sigma)-exp(-beta*P*lambda* ...
      sigma))+exp(-beta*P*kappa*sigma))/(beta*P) + ...
      (-exp(-beta*kappa*P*sigma)+exp(-beta*lambda*P*sigma+beta*epsilon ...
      *(-kappa+lambda)*zeta))/fenmu;
    end
    pfunc = sqrt(1/beta)*pfunc;
end