% draw (\beta p - rho)/rho^2 vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,2.5,4,1,2];
betas=[1];
suffix = '.mat';
%ps=[linspace(0.1,1,10), 1.5, linspace(2,10,9)];
ps=2.0;
%ps=10.^linspace(-5,-1.2,20);
rq=1;
Nbeta = size(betas,2);
rho = zeros(size(ps,1),Nbeta); % NNN
for beta=betas
  rp=1;
  for p=ps/beta
    rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 100))),p, 1e-5, 'central2');
    rp=rp+1;
  end
  rq=rq+1;
end

y = (ps./rho-1)./rho;

Nbeta = size(betas,2);
figure; hold on;
for rp=1:Nbeta
    plot(ps,rho(:,rp));
end
xlabel('\beta p');ylabel('\rho');
legend(string(betas),'Location','northwest');

figure;
plot(rho,y);
xlabel('\rho');ylabel('(\beta p - \rho)/\rho^2');
set(gca,'xscale','log');set(gca,'yscale','log');
legend(string(betas),'Location','northwest');

