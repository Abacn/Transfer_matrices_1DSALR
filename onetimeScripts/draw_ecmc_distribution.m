function draw_ecmc_distribution()
% draw rho-h(rho)
% coeffs = [1, 2.5, 4, 1, 1];
close all;
addpath('..');
choice=2;
lambda = 5;
if(lambda==2)
    fin = '../../simulationdata/ecmc_nvt/secondarydata/cdf/1,2.2,4,1,1,0.5.dat';
elseif(lambda==5)
    fin = '../../simulationdata/ecmc_nvt/secondarydata/cdf/1,2.5,4,1,1,0.5.dat';
end
cdata = load(fin);
rhos = [0.01 0.05 0.1 0.2 0.3];
nrho=length(rhos);
xx=[1:0.1:10].';
yys=zeros(length(xx), nrho);
figure;hold on;
for ind=1:nrho
    plot(nan, nan, 'x-');
end
set(gca,'ColorOrderIndex',1);
for ind=1:nrho
    filted = find(cdata(:,1) == rhos(ind));
    yys(:,ind)=spline(cdata(filted,2), cdata(filted,3), xx);
    plot(cdata(filted,2), cdata(filted,3), 'x');
end
set(gca,'ColorOrderIndex',1);
for ind=1:nrho
    plot(xx, yys(:,ind),'-');
end
xlabel('N');ylabel('p(N)');
legend('\rho=0.01','\rho=0.05','\rho=0.1','\rho=0.2','\rho=0.3');
legend boxoff;
xlim([0.9 10.1]);
ylim([0 0.6]);
set(gca,'XTick',linspace(1,9,5));
set(gca,'YTick',linspace(0,0.6,4));
set(gca,'fontsize',20);
end