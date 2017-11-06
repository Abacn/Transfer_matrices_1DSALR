% draw comparison between calculation and simulation
% choice:
% 1-plot simulation data, big figure
% 2-plot simulation data, small figure (deprecated)
close all;
choice = 1;
si_data=dlmread('../../simulationdata/ecmc_nvt/secondarydata/pressure/1,2.5,4,1,1.dat','',1,0);
cv_data=dlmread('../../simulationdata/mc_nvt/secondarydata/pressure/1,2.5,4,1,1.dat','',1,0);
vm_data=dlmread('../../simulationdata/vmmc/secondarydata/pressure/1,2.5,4,1,1.dat','',1,0);
tm_data = dlmread('../data_rho/1.0_2.5_4.0_1.0_1.00.dat','', 1, 0);
tm_0_3 = find(tm_data(:,1)==0.3);
tm_0_5 = find(tm_data(:,1)==0.5);
tm_0_2 = find(tm_data(:,1)==0.2);
% plot calculation data
psize=2;unifyfigure;hold on;
plot(tm_data(tm_0_5, 3), tm_data(tm_0_5, 2),'r', tm_data(tm_0_3, 3),tm_data(tm_0_3, 2),'b',tm_data(tm_0_2, 3), tm_data(tm_0_2, 2),'m');

% plot simulation data
box off;
if(choice==1)  % big figure
    si_0_2 = and(si_data(:,1)==0.2, si_data(:,2)>=0.05 | si_data(:,2)<=0.01);
    si_0_3 = and(si_data(:,1)==0.3, si_data(:,2)>=0.05 | si_data(:,2)<=0.01);
    si_0_5 = and(si_data(:,1)==0.5,si_data(:,2)>=0.05 | si_data(:,2)<=0.01);
    %si_1 = and(si_data(:,1)==1,si_data(:,2)>=0.05 | si_data(:,2)<=0.01);
    vm_0_2 = and(vm_data(:,1)==0.2, vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    vm_0_3 = and(vm_data(:,1)==0.3, vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    vm_0_5 = and(vm_data(:,1)==0.5,vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    %vm_1 = and(vm_data(:,1)==1,vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    plot(si_data(si_0_5, 2), si_data(si_0_5, 3)/0.5, 'or',si_data(si_0_3, 2), si_data(si_0_3, 3)/0.3, 'ob',si_data(si_0_2, 2),si_data(si_0_2, 3)/0.2, 'om', 'markers', 8);
    plot(vm_data(vm_0_5, 2), vm_data(vm_0_5, 3)/0.5, 'xr',vm_data(vm_0_3, 2), vm_data(vm_0_3, 3)/0.3, 'xb', vm_data(vm_0_2, 2),vm_data(vm_0_2, 3)/0.2,'xm', 'markers', 8);
    xlim([0 0.5]);
    legend('T=0.5','T=0.3','T=0.2', 'Location','northwest');
    legend boxoff;
    xlabel('\rho'); ylabel('$\beta p$', 'Interpreter', 'Latex');
    set(gca,'fontsize',16);
    set(gca,'XTick',linspace(0,0.4,3))
    set(gca,'YTick',linspace(0,0.4,3))
    print(fig,'render/algocomp.eps','-depsc');
elseif (choice==2) % small figure
    si_0_3 = and(si_data(:,1)==0.3, si_data(:,3)<=0.1);
    si_0_5 = and(si_data(:,1)==0.5,si_data(:,3)<=0.1);
    si_0_2 = and(si_data(:,1)==0.2,si_data(:,3)<=0.1);
    plot(si_data(si_0_5, 2), si_data(si_0_5, 3)/0.5, 'or',si_data(si_0_3, 2), si_data(si_0_3, 3)/0.3, 'ob',si_data(si_0_2, 2),si_data(si_0_2, 3)/0.2, 'om');
    xlim([0 0.1]);
    set(gca,'YTick',linspace(0,0.3,4))
    set(gca,'XTick',linspace(0,0.1,3))
    ax = gca;
    ax.Box = 'off';
    set(gca,'fontsize',20);
    ax.XAxisLocation = 'origin';
elseif (choice==3)  % test desk
    cv_0_2 = and(cv_data(:,1)==0.2, cv_data(:,2)>=0.05 | cv_data(:,2)<=0.01);
    cv_0_3 = and(cv_data(:,1)==0.3, cv_data(:,2)>=0.05 | cv_data(:,2)<=0.01);
    cv_0_5 = and(cv_data(:,1)==0.5,cv_data(:,2)>=0.05 | cv_data(:,2)<=0.01);
    cv_1 = and(cv_data(:,1)==1,cv_data(:,2)>=0.05 | cv_data(:,2)<=0.01);
    vm_0_2 = and(vm_data(:,1)==0.2, vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    vm_0_3 = and(vm_data(:,1)==0.3, vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    vm_0_5 = and(vm_data(:,1)==0.5,vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    vm_1 = and(vm_data(:,1)==1,vm_data(:,2)>=0.05 | vm_data(:,2)<=0.01);
    plot(cv_data(cv_0_5, 2), cv_data(cv_0_5, 3)/0.5, 'or',cv_data(cv_0_3, 2), cv_data(cv_0_3, 3)/0.3, 'ob',cv_data(cv_0_2, 2),cv_data(cv_0_2, 3)/0.2, 'om');
    plot(vm_data(vm_0_5, 2), vm_data(vm_0_5, 3)/0.5, 'xr',vm_data(vm_0_3, 2), vm_data(vm_0_3, 3)/0.3, 'xb', vm_data(vm_0_2, 2),vm_data(vm_0_2, 3)/0.2,'xm');
    xlim([0 0.5]);
    legend('T=0.5','T=0.3','T=0.2', 'Location','northwest');
    legend boxoff;
    xlabel('\rho'); ylabel('\beta p');
    set(gca,'fontsize',20);
    set(gca,'XTick',linspace(0,0.5,6))
    set(gca,'YTick',linspace(0,1,6))

end
%set(gca,'XTick',linspace(0,5,6))
