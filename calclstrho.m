% Calculate density the of total cluster vs. pressure
% Calculate the cumulative gap distribution function at \lambda \sigma.

% option: 1-total cluster; 2-free particle
option = 2;

coeffs = [1, 2.5, 4, 1, 1];
T = 0.2;
beta = 1/T;
bps=[10.^linspace(-5, 0, 20)];
ps = bps*T;
nps = length(ps);
rhos = zeros(size(ps));    % density
rhocls = zeros(size(ps));  % density of total cluster
rp = 1;
while rp<=nps
    rhos(rp) = findrho(ps(rp), beta, coeffs);
    if 1==option
        [rlist, pros] = gapdf(ps(rp), beta, coeffs);
        procl = sum(pros(rlist > coeffs(1)*coeffs(2)))*(rlist(2)-rlist(1));
        rhocls(rp) = procl*rhos(rp);
    else
        [ks, ~] = cdf(ps(rp), beta, coeffs, 200, 1);
        rhocls(rp) = rhos(rp)*ks(1);
    end
    rp = rp + 1;
end

loglog(rhos, rhocls);