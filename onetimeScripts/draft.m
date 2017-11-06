% draft.m
option = 1;
if 1==option
%% Rend the relaxation curve
close all;
thisdata = alldata{3};
upper = 8;
x=linspace(0,upper);
ndata = size(thisdata,1);
figure;
figure('rend','painters','pos',[10 10 400 300]);
hold on;
for rp=1:ndata
    N = thisdata(rp,1);
    tau = thisdata(rp,2);
    y = func(tau/(N^2*log(N)), x);
    plot(x, y);
end
set(gca, 'XTick', linspace(0,upper,3));
set(gca, 'YTick', linspace(0.2,1,3));
set(gca,'fontsize',20);
xlabel('t/N^2 log(N)'); ylabel('C_{bind}');
elseif 2==option
%% 
end