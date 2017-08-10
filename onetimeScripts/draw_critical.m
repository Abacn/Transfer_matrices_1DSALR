function draw_critical()
% draw rho-h(rho)
% coeffs = [1, 2.5, 4, 1, 1];
close all;
addpath('..');
choice=2;
lambda = 2;
if(choice==1)
    if(lambda==2)
        rhodata=load('../pyout/rhoc-1,2.2,4,1.dat');
    elseif(lambda==5)
        rhodata=load('../pyout/rhoc-1,2.5,4,1.dat');
    end
    pickedXi = [0.0, 0.1, 0.5, 1, 4];
    lt = length(pickedXi);
	figure; hold on;
	for ind=1:lt
        tdata = rhodata(rhodata(:,1) == pickedXi(ind), 2:3); % data of specific temperature
        plot2(tdata);
    end
    if(lambda==2)
        xlim([0.06 0.3])
        ylim([1e-8, 2e-2]);
        set(gca,'XTick',linspace(0.06,0.3,4));
    elseif(lambda==5)
        xlim([0.08 0.38]);
        ylim([1e-8, 2e-2]);
        set(gca,'XTick',linspace(0.08,0.38,4));
    end
    set(gca,'yscale','log');
    set(gca,'fontsize',16);
    set(gca,'YTick',10.^linspace(-8,-2,4));
    legend('\xi=0','\xi=0.1', '\xi=0.5','\xi=1','\xi=4' ,'Location','southeast');
    legend boxoff;
    xlabel('T'); ylabel('\rho_{cmc}');
% mark the minimum
elseif(choice==2)
 figure; hold on;
 rhodata=load('../pyout/Tc-1,2.5,4,1.dat');
 plot(rhodata(:,1), rhodata(:,2),'b');
  % Tc given by virial coefficents
 rhodata=dlmread('../data_Bs/Tc_1.0_2.5_4.0_1.0.dat','',1,0);
 plot(rhodata(:,1), rhodata(:,2),'--b');
 rhodata=load('../pyout/Tc-1,2.2,4,1.dat');
 plot(rhodata(:,1), rhodata(:,2),'r');
  rhodata=dlmread('../data_Bs/Tc_1.0_2.2_4.0_1.0.dat','',1,0);
 plot(rhodata(:,1), rhodata(:,2),'--r');
 legend('\lambda=2.5','','\lambda=2.2','');
 legend boxoff;
 xlim([0 4]);
 ylim([0.15 0.4]);
 set(gca,'XTick',linspace(0,4,5));
 set(gca,'YTick',linspace(0.15,0.4,6));
 set(gca,'fontsize',16);
 xlabel('\xi');
 ylabel('T_c');
elseif(choice==3)
 if(lambda == 5)
 fin = '../data_Bs/B3_1.0_2.5_4.0_1.0.dat';
 elseif(lambda == 2)
 fin = '../data_Bs/B3_1.0_2.2_4.0_1.0.dat';
 end
 xis = [0, 0.1, 0.5, 1, 4];
 bdata = dlmread(fin, '', 1, 0);
 figure; hold on;
 for rp=1:length(xis)
     filted = find(bdata(:,1) == xis(rp));
     plots(bdata(filted,2), bdata(filted,3), 2);
 end
 xlabel('T');ylabel('B_3(T)');
 if(lambda == 5)
    xlim([0.2 0.4]);
    set(gca,'YTick',linspace(-6,3,4));
    set(gca,'YTickLabel',['-10^6';'-10^3';'    0';' 10^3']);
 elseif(lambda == 2)
    xlim([0.1 0.3]);
    set(gca,'YTick',linspace(-12,6,4));
    set(gca,'YTickLabel',['-10^{12}';'   -10^6';'       0';'    10^6']);
 end
 set(gca,'fontsize',16);
 legend('\xi=0','\xi=0.1', '\xi=0.5','\xi=1','\xi=4', 'Location', 'southeast');
 legend boxoff;
end
end

function minpoint = plots(x, y, choice)  % log scale (-inf, -10) ~ linear (-10, 10) ~ log(10, inf)
    if(nargin < 3)
        choice = 1;
    end
    f1 = y >= -10 & y <=10;
    if(choice==2)
     y(f1) = y(f1)/10;
     y(y < -10) = -log10(-y(y < -10));
     y(y > 10) = log10(y(y > 10));
    end
    [~, ind] = min(y);
    plot(x, y);
    minpoint = [x(ind), y(ind)];
end