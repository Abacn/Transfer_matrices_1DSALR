% check convergence of m
close all;
addpath('..');
coeffs = [1,2.5,4,1,1];
P = 1e-2;
T = 0.2;
beta = 1/T;

xstart = 103;
xinterval = 12;
xs = xstart:xinterval:420;
ys = [];
x = 0;
for N=xs
    ys(x+1) = findrho(P, 1/T, coeffs, N);
    x = x+1;
end

plot(1./xs, ys);
xlim([0, 1/xs(1)]);

save('onetimeData/checkm_mod1.mat', 'xs', 'ys', 'coeffs', 'P', 'T');