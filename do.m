% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,2.5,4,1,1];
beta=1;
% targetp = 2.0;ps = [-0.02 -0.01 0 0.01 0.02]+targetp;
ps=[0.01,linspace(0.1,5,20)];
rp=1;
rst = zeros(size(ps));    % next next nearest neighbor
rst_NNN = zeros(size(ps)); % nearest neighbor
rst_NN = zeros(size(ps)); % nearest neighbor
rst_HS = zeros(size(ps)); % hard sphere
for p=ps/beta
    rst_NNN(rp) = Pfunc_isobaric(p, beta, coeffs);
    rst(rp) = Pfunc_isobaric_3NN(p, beta, coeffs);
    rst_NN(rp) = Pfunc_isobaric_NN(p, beta, coeffs);
    rst_HS(rp) = Pfunc_isobaric_HS(p, beta, coeffs);
    rp=rp+1;
end
% rst is partition function
% calculate Gibbs free energy (G/n*beta)
gibbs = -log(rst);
gibbs_NN = -log(rst_NN);
gibbs_NNN = -log(rst_NNN);
gibbs_HS = -log(rst_HS);
plot(ps,gibbs,ps,gibbs_NNN,ps,gibbs_NN,ps,gibbs_HS);
xlabel('\beta p');ylabel('\beta G/N');
legend('3NN','NNN','NN','HS','Location','northwest');

% pressure - density
rho = gradient(ps)./gradient(gibbs);
rho_HS = ps ./ (coeffs(1)*ps + 1);
rho_NN = gradient(ps)./gradient(gibbs_NN);
rho_NNN = gradient(ps)./gradient(gibbs_NNN);
figure;
plot(ps,rho,ps,rho_NNN,ps,rho_NN,ps,rho_HS);
xlabel('\beta p');ylabel('\rho');
legend('3NN','NNN','NN','HS','Location','northwest');
% gibbs_HS = gibbs; ps_HS=ps; rho_HS=rho;
