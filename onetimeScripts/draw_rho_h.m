% option
% 1: draw rho-h(rho) curve
% 2: draw rho-p curve

function draw_rho_h(option)
% draw rho-h(rho)
% coeffs = [1, 2.5, 4, 1, 1];
close all;
addpath('..');
if nargin<1
    option = 2;
end
choice=4;
if(choice==1)
 fname = '../data_min/1.0_2.5_4.0_1.0_1.0.dat';
elseif(choice==2)
 fname = '../data_min/1.0_2.2_4.0_1.0_1.0.dat';
elseif(choice==3)
 fname = '../data_min/1.0_2.0_4.0_1.0_1.00.dat';
elseif(choice==4)
 fname = '../data_min/1.0_2.5_4.0_1.0_0.05.dat';
else
 return;
end
if 1==option
dn = chopnum(fname); % get data line
rhodata=dlmread(fname,'',[dn(1) 0 dn(2) 3]);
pickedT = [0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.3];
lt = length(pickedT);
psize=3;unifyfigure; hold on;
minpoint = zeros(lt,2);
for ind=1:lt
    tdata = rhodata(rhodata(:,1) == pickedT(ind), 2:4); % data of specific temperature
    [~, inds] = sort(tdata(:,2));
    stdata = tdata(inds, :);     % sorted
    minpoint(ind,:) = plots(stdata(:,2), stdata(:,3));
end
set(gca,'XTick',[1e-8 1e-6 1e-4 1e-2 1e0]);
set(gca,'YTick',[-1e6 -1e4 -1e2 -1e0]);
xlabel('\rho'); ylabel('h(\rho)');
% mark the minimum
set(0,'defaultfigurecolor',[0 0 0]);
for ind=1:lt
    if(choice==1)
        tthres = 0.28;
    elseif(choice==2)
        tthres = 0.18;
    else
        tthres = 0;
    end
    if(pickedT(ind) <=tthres)
        plot(minpoint(ind,1), minpoint(ind,2), '*');
    end
end
legend('T=0.10','T=0.12','T=0.14','T=0.16','T=0.20','T=0.24', 'T=0.30', 'Location','southeast');
legend boxoff;
xlim([1e-8 1e0]);
ylim([-1e6 -1]);
set(gca,'xscale','log');
set(gca,'yscale','log');
elseif 2==option
dn = chopnum(fname); % get data line
rhodata=dlmread(fname,'',[dn(1) 0 dn(2) 3]);
pickedT = [0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.3];
lt = length(pickedT);
figure; hold on;  % figure
minpoint = zeros(lt,2);
for ind=1:lt
    tdata = rhodata(rhodata(:,1) == pickedT(ind), 2:4); % data of specific temperature
    [~, inds] = sort(tdata(:,2));
    stdata = tdata(inds, :);     % sorted
    minpoint(ind,:) = plots(stdata(:,2), stdata(:,1));
end
set(gca,'XTick',[1e-8 1e-6 1e-4 1e-2 1e0]);
set(gca,'YTick',[-1e6 -1e4 -1e2 -1e0]);
xlabel('\rho'); ylabel('$\beta p$','Interpreter', 'LaTeX');
% mark the minimum
set(0,'defaultfigurecolor',[0 0 0]);
for ind=1:lt
    if(choice==1)
        tthres = 0.28;
    elseif(choice==2)
        tthres = 0.18;
    else
        tthres = 0;
    end
    if(pickedT(ind) <=tthres)
        plot(minpoint(ind,1), minpoint(ind,2), '*');
    end
end
legend('T=0.10','T=0.12','T=0.14','T=0.16','T=0.20','T=0.24', 'T=0.30', 'Location','southeast');
legend boxoff;
%xlim([1e-8 1e0]);
%ylim([-1e6 -1]);
set(gca,'xscale','log');
set(gca,'yscale','log');
end
end

function minpoint = plots(x, y)  % log scale (-inf, -10) ~ linear (-10, 10) ~ log(10, inf)
    f1 = y >= -10 & y <=10;
    %y(f1) = y(f1)/10;
    %y(y < -10) = -log10(-y(y < -10));
    %y(y > 10) = log10(y(y > 10));
    [~, ind] = min(y);
    plot(x, y);
    minpoint = [x(ind), y(ind)];
end