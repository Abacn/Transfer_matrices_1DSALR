% draw comparison between calculation and simulation
close all;
choice = 2;
si_data=load('../../simulationdata/ecmc_nvt/secondarydata/p-1,2.5,4,1,1.dat');
tm_data = dlmread('../data_rho/1.0_2.5_4.0_1.0_1.00.dat','', 1, 0);
tm_0_3 = find(tm_data(:,1)==0.3);
tm_0_5 = find(tm_data(:,1)==0.5);
tm_1 = find(tm_data(:,1)==1);
% plot calculation data
plot(tm_data(tm_0_3, 2),tm_data(tm_0_3, 3),'b',tm_data(tm_0_5, 2), tm_data(tm_0_5, 3),'r', tm_data(tm_1, 2), tm_data(tm_1, 3),'k');
hold on;
% plot simulation data

if(choice==1)  % big figure
    si_0_3 = and(si_data(:,1)==0.3, si_data(:,2)>0.05 | si_data(:,2)<=0.01);
    si_0_5 = and(si_data(:,1)==0.5,si_data(:,2)>0.05 | si_data(:,2)<=0.01);
    si_1 = and(si_data(:,1)==1,si_data(:,2)>0.05 | si_data(:,2)<=0.01);
    plot(si_data(si_0_3, 3)/0.3,si_data(si_0_3, 2),'ob',si_data(si_0_5, 3)/0.5, si_data(si_0_5, 2), 'or',si_data(si_1, 3), si_data(si_1, 2), 'ok');
    xlim([0 5]);
    legend('T=0.3','T=0.5','T=1.0', 'Location','southeast');
    legend boxoff;
    xlabel('\beta p'); ylabel('\rho');
    set(gca,'fontsize',16);
    set(gca,'YTick',linspace(0,1,6))
elseif (choice==2) % small figure
    si_0_3 = and(si_data(:,1)==0.3, si_data(:,3)<=0.1);
    si_0_5 = and(si_data(:,1)==0.5,si_data(:,3)<=0.1);
    si_1 = and(si_data(:,1)==1,si_data(:,3)<=0.1);
    plot(si_data(si_0_3, 3)/0.3,si_data(si_0_3, 2),'ob',si_data(si_0_5, 3)/0.5, si_data(si_0_5, 2), 'or',si_data(si_1, 3), si_data(si_1, 2), 'ok');
    xlim([0 0.1]);
    set(gca,'YTick',linspace(0,0.3,4))
    set(gca,'XTick',linspace(0,0.1,3))
    ax = gca;
    ax.Box = 'off';
    set(gca,'fontsize',16);
    ax.XAxisLocation = 'origin';
end
%set(gca,'XTick',linspace(0,5,6))
