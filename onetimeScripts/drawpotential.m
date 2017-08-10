addpath('..');
xright = 6;
a = linspace(0.999, xright, 1000);
coeffs1 = [1, 2.5, 4, 1, 1];
coeffs2 = [1, 2.5, 4, 1, 0.5];
coeffs3 = [1, 2.5, 4, 1, 0.1];
b1=potential(a, coeffs1);
b2=potential(a, coeffs2);
b3=potential(a, coeffs3);
b1(1) = 999;
b2(1) = 999;
b3(1) = 999;
plot(a,b1,a,b2,a,b3);
xlim([0 xright]);
ylim([-1.5 2]);
ax = gca;
ax.Box = 'off';
ax.XAxisLocation = 'origin';
xlabel('r');
ylabel('u(r)');
legend('\xi=1', '\xi=0.5', '\xi=0.1');
legend boxoff;
set(gca,'fontsize',16);

