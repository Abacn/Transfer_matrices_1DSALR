% draw p vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,2.5,4,1,2];
betas=[2];
Foldername = 'data_low';
suffix = '.mat';
% use pressure, not beta*pressure here
pr=[linspace(0.01,1,10), 1.5, linspace(2, 10, 9)];
% pr=1;
% pr=10.^linspace(-4,-1,20).';
rq=1;
Nbeta = size(betas,2);
rho = zeros(size(pr,1),Nbeta); % NNN
for beta=betas
  rp=1;
  for p=pr
    rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 150))),p);
    rp=rp+1;
  end
  rq=rq+1;
end

y=zeros(size(rho));
for i=1:size(betas,2)
    y(:,i) = (pr(:)*betas(i)./rho(:,i)-1)./rho(:,i);
end

Nbeta = size(betas,2);
figure; hold on;
for rp=1:Nbeta
    plot(pr,rho(:,rp));
end
xlabel('p');ylabel('\rho');
legend(string(betas),'Location','northwest');

figure;
plot(rho,y);
xlabel('\rho');ylabel('(\beta p - \rho)/\rho^2');
set(gca,'xscale','log');set(gca,'yscale','log');
legend(string(betas),'Location','northwest');
% original_color_order=get(0,'DefaultAxesColorOrder');
% set(0,'DefaultAxesColorOrder',original_color_order);
