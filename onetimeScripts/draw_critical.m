% choice
% 1: draw T-rho_c, cmc with T
% 2: draw T-rho_c where B_3 changes sign.
% 3: draw terminal clustering temperature (x)
% 4: draw cmc from cdf
% 5: draw terminal clustering temperature by CDF (x)
% 6: draw terminal clustering temperature from h and CDF

function draw_critical(choice, lambda)
% draw rho-h(rho)
% coeffs = [1, 2.5, 4, 1, 1];
close all;
addpath('..');
if nargin<2
    lambda = 5;
end
if nargin<1
    choice = 6;
end
if(choice==1)
    if(lambda==2)
        rhodata=load('../pyout/rhoc-1,2.2,4,1.dat');
    elseif(lambda==5)
        rhodata=load('../pyout/rhoc-1,2.5,4,1.dat');
    end
    pickedXi = [0.0, 1, 2];
    lt = length(pickedXi);
	psize=2;unifyfigure;
	for ind=1:lt
        tdata = rhodata(rhodata(:,1) == pickedXi(ind), [3,2]); % data of specific temperature
        plot2(tdata);
        if 1==ind
             hold on;
        end
    end
    set(gca, 'ColorOrderIndex', 1);
    	for ind=1:lt
        tdata = rhodata(rhodata(:,1) == pickedXi(ind), [3,2]); % data of specific temperature
        plot(tdata(end,1),tdata(end,2), '*', 'Markersize', 8);
    end
    if(lambda==2)
        ylim([0.06 0.3])
        xlim([1e-8, 2e-2]);
        set(gca,'XTick',10.^linspace(-8,-2,4));
        set(gca,'YTick',linspace(0.06,0.3,4));
    elseif(lambda==5)
        ylim([0.0 0.4]);
        xlim([1e-7, 1e-1]);
        set(gca,'XTick',10.^linspace(-7,-1,4));
        set(gca,'YTick',linspace(0.0,0.4,3));
    end
    set(gca,'xscale','log');
    %set(gca,'fontsize',20);
    %legend('\xi=0','\xi=1','\xi=2' ,'Location','southeast');
    ylabel('$T$', 'Interpreter', 'LaTeX');
    xlabel('$\rho_{ccd}[h]$', 'Interpreter', 'LaTeX');
    %legend boxoff;
    unifyfigureend;
    if 2==lambda
        print(fig,'render/rhoc-2.2.eps','-depsc');
    elseif 5==lambda
        print(fig,'render/rhoc-2.5.eps','-depsc');
    end
% draw B3 with T
elseif(choice==2)
 if(lambda == 5)
 fin = '../data_Bs/B3_1.0_2.5_4.0_1.0_1.dat';
 elseif(lambda == 2)
 fin = '../data_Bs/B3_1.0_2.2_4.0_1.0.dat';
 end
 xis = [0, 0.1, 0.5, 1, 4];
 bdata = dlmread(fin, '', 1, 0);
 psize=2;unifyfigure; hold on;
 for rp=1:length(xis)
     filted = find(bdata(:,1) == xis(rp));
     % plots(bdata(filted,2), bdata(filted,3), 2);
     plot(bdata(filted,2), bdata(filted,3));
 end
 xlabel('$T$','Interpreter', 'LaTeX');
 if(lambda == 5)
    ylabel('$B_3(T)$','Interpreter', 'LaTeX');  % two plot together
    xlim([0.25 0.45]);
    plot([0.25 0.45],[0 0],'-k');
    set(gca,'XTick',[0.25,0.35,0.45]);
    ylim([-2000 1000]);
    set(gca,'YTick',[-2000,-1000,0,1000]);
    %set(gca,'YTick',linspace(-8,4,4));
    %set(gca,'YTickLabel',['-10^8';'-10^4';'    0';' 10^4']);
 elseif(lambda == 2)
    xlim([0.15 0.3]);
    set(gca,'XTick',[0.2,0.3]);
    plot([0.15 0.3],[0 0],'-k');
    ylim([-20000 20000]);
    %ylim([-12 6]);
    set(gca,'YTick',[-20000,0, 20000]);
    %set(gca,'YTickLabel',['-10^{12}';'   -10^6';'       0';'    10^6']);
 end
 legend('\xi=0','\xi=0.1', '\xi=0.5','\xi=1','\xi=4', 'Location', 'southeast');
 legend boxoff;
 if 2==lambda
     print(fig,'render/B3-2.2-v2.eps','-depsc');
 elseif 5==lambda
     print(fig,'render/B3-2.5-v2.eps','-depsc');
 end
% draw terminal clustering temperature
elseif(choice==3)
 psize=2;unifyfigure; hold on;
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
 xlim([0 4]);
 ylim([0.15 0.4]);
 set(gca,'XTick',linspace(0,4,5));
 set(gca, 'YTick', [0.2,0.3,0.4]);
 legend boxoff;
 
 xlabel('\xi');
 ylabel('$T_t$','Interpreter', 'LaTeX');
 print(fig,'render/Tc.eps','-depsc');
% draw cmc from cdf
elseif(choice==4)
 pickedXi = {'0.00';'1.00'; '2.00'};
 sprefix  = strcat('../data_cdf/t_1.0_2.', num2str(lambda) ,'_4.0_1.0_');
 fin = cellfun(@(str)strcat(sprefix, str, '.dat'), pickedXi, 'UniformOutput', false);
 lenfin = length(fin);
 psize=2;unifyfigure;
 set(gca,'ColorOrderIndex',1);
 for rp=1:lenfin
    plot(nan, nan);
    if(1==rp)
        hold on;
    end
 end
 for rp=1:lenfin
    fname = fin{rp};
    data=dlmread(fname,'',1,0);
    [index, content] = changpoint(data(:,4));
    index = [1, index];
    ntypes = length(content);
    for rq=1:ntypes
        nowtype = content(rq);
        if 1==nowtype
            mark = ':';
        elseif 2==nowtype
            mark = '--';
        else % 3
            mark = '-';
        end
        nowrange = index(rq):index(rq+1);
        set(gca,'ColorOrderIndex',rp);
        plot(data(nowrange, 3), data(nowrange, 1), mark);
        if 3==nowtype
            set(gca,'ColorOrderIndex',rp);
            plot(data(nowrange(end), 3), data(nowrange(end), 1), '*', 'markers',10);
        end
    end
 end
 set(gca,'xscale','log');
 ylim([0 0.4]);
 xlim([1e-7, 1e-1]);
 set(gca,'YTick',linspace(0,0.8,5));
 set(gca,'XTick',10.^linspace(-7,-1,4));
 ylabel('$T$','Interpreter', 'LaTeX');
 %legend('\xi=0','\xi=1', '\xi=2','Location', 'southeast');
 %legend boxoff;
 unifyfigureend;
 if 2==lambda
        print(fig,'render/rhoccdf-2.2.eps','-depsc');
 elseif 5==lambda
        xlabel('$\rho_{ccd}$[CDF]','Interpreter', 'LaTeX');
        print(fig,'render/rhoccdf-2.5.eps','-depsc');
 end

% draw terminal clustering temperature by CDF
elseif 5==choice
 restat = false;
 if restat
     data25=tmtcdf('t_1.0_2.5_4.0_1.0');
     data22=tmtcdf('t_1.0_2.2_4.0_1.0');
 else
    data25=dlmread('../data_cdf/Tc_1.0_2.5_4.0_1.0.dat','',1,0);
    data22=dlmread('../data_cdf/Tc_1.0_2.2_4.0_1.0.dat','',1,0);
 end
 alldata = {data25; data22};
 lenfin = length(alldata);
 psize=2;unifyfigure; hold on;
 for rp=1:lenfin
    plot(nan, nan);
 end
 for rp=1:lenfin
     nowdata = alldata{rp};
     [index, content] = changpoint(nowdata(:,3));
     index = [1, index];
     ntypes = length(content);
     for rq=1:ntypes
        nowtype = content(rq);
        if 1==nowtype
            mark = ':';
        elseif 2==nowtype
            mark = '--';
        else % 3
            mark = '-';
        end
        nowrange = index(rq):index(rq+1);
        set(gca,'ColorOrderIndex',rp);
        plot(nowdata(nowrange, 1), nowdata(nowrange, 2), mark);
     end
 end
 ylim([0.15 0.6]);
 xlabel('\xi');ylabel('$T_{t-CDF}$', 'Interpreter', 'LaTeX');
 %legend('\lambda=2.5','\lambda=2.2');
 %legend boxoff;
 set(gca,'XTick',linspace(0,4,5));
 set(gca,'YTick',linspace(0.2,0.6,5));
  print(fig,'render/Tccdf.eps','-depsc');
 %end  % debug. should be commented
elseif 6==choice
 restat = false;
 if 5==lambda
     rhodata=load('../pyout/Tc-1,2.5,4,1.dat');
     Bdata = dlmread('../data_Bs/Tc_1.0_2.5_4.0_1.0.dat','',1,0);
     if restat
        data2 = tmtcdf('t_1.0_2.5_4.0_1.0');
     else
         data2 = dlmread('../data_cdf/Tc_1.0_2.5_4.0_1.0.dat','',1,0);
     end
 elseif 2==lambda
     rhodata=load('../pyout/Tc-1,2.2,4,1.dat');
     Bdata = dlmread('../data_Bs/Tc_1.0_2.2_4.0_1.0.dat','',1,0);
     if restat
        data2 = tmtcdf('t_1.0_2.2_4.0_1.0');
     else
         data2 = dlmread('../data_cdf/Tc_1.0_2.2_4.0_1.0.dat','',1,0);
     end
 end
 psize=2;unifyfigure;
 %plot(nan, nan, '-k'); hold on;  plot(nan, nan, '--k'); plot(nan, nan, ':k'); plot(nan, nan, '-.k');
 plot(rhodata(:,1), rhodata(:,2),'b');
 hold on;
 plot(Bdata(:,1), Bdata(:,2),'--b');
 colors = {'b';'r'};
 [index, content] = changpoint(data2(:,3));
 index = [1, index];
 ntypes = length(content);
 for rq=1:ntypes
     nowtype = content(rq);
     if 1==nowtype
         mark = ':';
     elseif 2==nowtype
         mark = '--';
     else % 3
         mark = '-.';
     end
     nowrange = index(rq):index(rq+1);
     plot(data2(nowrange, 1), data2(nowrange, 2), strcat(mark, colors{1}));
     if rq==ntypes
         plot(data2(nowrange(1), 1), data2(nowrange(1), 2), strcat('*', colors{1}), 'markers', 10);
     end
 end
 
 xlim([0 2]);
 ylim([0.15 0.6]);
 xlabel('\xi');ylabel('$T_{t}$', 'Interpreter', 'LaTeX');
 legend({'$T_t[h]$','$T_t[B_3]$','$T_t$[CDF-condensation]', '$T_t$[CDF-clustering]'}, 'Interpreter', 'LaTeX');
 legend boxoff;
 set(gca,'XTick',linspace(0,2,5));
 set(gca,'YTick',linspace(0.2,0.6,5));
 unifyfigureend;
  print(fig,strcat('render/Tccdf-l2.',num2str(lambda),'.eps'),'-depsc');
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

function [index, content] = changpoint(data)
    % cut index marks
    % Ex: [1,1,1,1,2,2,3,4,4] will return index=[4,6,7,9], content=[1,2,3,4]
    index = [];
    content = [];
    nownumber = data(1);
    for rp=2:length(data)
        if data(rp) ~= nownumber
            nownumber = data(rp);
            index = [index, rp-1];
            content = [content, data(rp-1)];
        end
    end
    index = [index, length(data)];
    content = [content, data(end)];
end

function dataout = tmtcdf(fileprefix)
    prefixes = '../data_cdf/';
    suffixes = '_*.dat';
    lsfiles = strsplit(ls(strcat(prefixes, fileprefix,suffixes)));
    nfnames = length(lsfiles);
    dataout = [];
    patt = 't_\d\.?\d*_\d\.?\d*_\d\.?\d*_\d\.?\d*_(\d\.?\d*)\.dat';
    for rp=1:nfnames
        fname = lsfiles{rp};
        tok = regexp(fname, patt,'tokens');
        if(isempty(tok))
            continue;
        end
        xi = str2double(tok{1});
        nowdata=sortrows(dlmread(strcat(prefixes,fname),'',1,0));
        ind = find(nowdata(:,4) == 2);
        Tc = nowdata(ind(1)-1, 1);
        Typec = nowdata(ind(1)-1, 4);
        dataout = [dataout; xi, Tc, Typec];
    end
    dataout = sortrows(dataout);
    fout = fopen(strcat(prefixes, 'Tc_', fileprefix(3:end), '.dat'), 'w');
    for rp=1:size(dataout, 1)
        fprintf(fout, '%.3f\t%.3f\t%d\n', dataout(rp,:));
    end
    fclose(fout);
end