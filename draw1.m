% draw (\beta p - rho)/rho^2 vs. rho lines
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,1.5,4,1,1];
beta=1;
%ps=linspace(0.1, 10, 50);
ps=10.^linspace(-4,0,50);
rp=1;
rho = zeros(size(ps));    % next next nearest neighbor
rho_NNN = zeros(size(ps)); % nearest neighbor
rho_NN = zeros(size(ps)); % nearest neighbor
rho_HS = zeros(size(ps)); % hard sphere
for p=ps/beta
    rho_NNN(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric(bp, beta, coeffs))),p);
    %rho(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs))),p);
    rho_NN(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_NN(bp, beta, coeffs))),p);
    rp=rp+1;
end
rho_HS = ps ./ (coeffs(1)*ps + 1);

y_NN = (ps./rho_NN-1)./rho_NN; % (beta*p-rho)/rho^2
y_NNN = (ps./rho_NNN-1)./rho_NNN;
y = (ps./rho-1)./rho; 
y_HS = coeffs(1)*(ps*coeffs(1)+1);

figure;
plot(ps,rho,ps,rho_NNN,ps,rho_NN,ps,rho_HS);
xlabel('\beta p');ylabel('\rho');
legend('3NN','NNN','NN','HS','Location','northwest');

figure;
plot(rho,y,rho_NNN,y_NNN,rho_NN,y_NN,rho_HS,y_HS);
xlabel('\rho');ylabel('(\beta p - \rho)/\rho^2');
set(gca,'xscale','log');
legend('3NN','NNN','NN','HS','Location','northwest');
% gibbs_HS = gibbs; ps_HS=ps; rho_HS=rho;
