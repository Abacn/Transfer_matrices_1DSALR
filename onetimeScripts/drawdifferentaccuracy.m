% Compare of different accuracy 
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
addpath('..');
choice=2;  % 1 normal, 2 small range
coeffs=[1,2.5,4,1,1];
beta=1;
% targetp = 2.0;ps = [-0.02 -0.01 0 0.01 0.02]+targetp;
rp=1;
if(~exist('done'))
    if(choice==1)
        ps=[0.01 linspace(0.1,5,20)];
    else
        ps=[0.001 linspace(0.01,0.1,20)];
    end
    rho = zeros(size(ps));    % next next nearest neighbor
    rho_NNN = zeros(size(ps)); % nearest neighbor
    rho_NN = zeros(size(ps)); % nearest neighbor
    rho_HS = zeros(size(ps)); % hard sphere
    for p=ps/beta
        rho(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p);
        rho_NNN(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric(bp, beta, coeffs, 200))),p);
        rho_NN(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_NN(bp, beta, coeffs))),p);
        rho_HS(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_HS(bp, beta, coeffs))),p);
        rp=rp+1;
    end
    done = 1;
end
figure;
plot(ps,rho_HS,ps,rho_NN,ps,rho_NNN,ps,rho);
if(choice==1)
  xlabel('\beta p');ylabel('\rho');
  legend('Hard sphere','Nearest neighbor','NN neighbor','NNN neighbor','Location','southeast');
  set(gca,'YTick',linspace(0,1,6))
else
  set(gca,'YTick',linspace(0,0.12,4));
  set(gca,'XTick',linspace(0,0.1,3));
  ax = gca;
  ax.Box = 'off';
  ax.XAxisLocation = 'origin';
end
legend boxoff;
set(gca,'fontsize',16);
% gibbs_HS = gibbs; ps_HS=ps; rho_HS=rho;
