% draw cdf
close all;

% option:
% 1 - compare results of different methods
% 2 - compare results at different density (transfer matrix), deprecated
% 3 - same as 2, log scale.
option = 1;
Tstr = '0.2';
parastr='1.0_2.5_4.0_1.0_1.00'; % 1.0_2.5_4.0_1.0_0.00
parastr1='1,2.5,4,1,1';
% clear 'done'
if 1==option  % compare results of different methods
    xs = 1:10;
    useT = str2double(Tstr);
    useRhos = [1e-4, 1e-3,1e-2, 1e-1]; % , 1e-1
    tmat = strcat('../data_cdf/',parastr,'.dat');
    ecmc = strcat('../../simulationdata/ecmc_nvt/secondarydata/cdf/',parastr1,'.dat');
    vmmc = strcat('../../simulationdata/vmmc/secondarydata/cdf/',parastr1,'.dat');
    cvmc = strcat('../../simulationdata/mc_nvt/secondarydata/cdf/',parastr1,'.dat');
    data_tmat = gettmat(tmat);
    data_ecmc = getsimu(ecmc);
    data_vmmc = getsimu(vmmc);
    data_cvmc = getsimu(cvmc);  % conventional monte carlo
    done = 1;
    tm_check = data_tmat(data_tmat(:,1)==useT & sum(data_tmat(:,2)==useRhos, 2), xs+2);
    ec_check = data_ecmc(data_ecmc(:,1)==useT & sum(data_ecmc(:,2)==useRhos, 2), xs+2);
    cv_check = data_cvmc(data_cvmc(:,1)==useT & sum(data_cvmc(:,2)==useRhos, 2), xs+2);
    vm_check = data_vmmc(data_vmmc(:,1)==useT & sum(data_vmmc(:,2)==useRhos, 2), xs+2);
    psize=2;unifyfigure;hold on;
    nrho = length(useRhos);
    set(gca,'ColorOrderIndex',1);
    plot(xs, tm_check, '-');
    plot(nan, nan, 'ko');
    plot(nan, nan, 'kx');
    set(gca,'ColorOrderIndex',1);
    for ind=1:nrho
        plot(xs, ec_check(ind,:), 'o','markers',8);
    end
    set(gca,'ColorOrderIndex',1);
    for ind=1:nrho
        plot(xs, vm_check(ind,:), 'X','markers',8);
    end
    legend('\rho=10^{-4}','\rho=10^{-3}','\rho=10^{-2}','\rho=10^{-1}', 'ECMC', 'CMMC');
    legend boxoff;
    xlim([1 10]); ylim([0 1]);
    xlabel('$n$', 'Interpreter', 'Latex');ylabel('$K(n)$', 'Interpreter', 'Latex');
    set(gca,'XTick',linspace(2,10,5));
    set(gca,'YTick',linspace(0,1,6));
    set(gca,'fontsize',16);
    
    % small figure
    curpos = get(gcf, 'Position');
    axnow=axes('Position', [0.26, 0.5, 0.35, 0.35]);
    for ind=1:nrho
        filted = tm_check(ind,:);
        plot(xs, filted, '-');
        if 1==ind
            hold on;
        end
    end
    set(gca,'ColorOrderIndex',1);
    for ind=1:nrho
        plot(xs, ec_check(ind,:), 'o','markers',6);
    end
    set(gca,'ColorOrderIndex',1);
    for ind=1:nrho
        plot(xs, vm_check(ind,:), 'X','markers',6);
    end
    box(axnow,'off');
    set(axnow, 'YScale' ,'log');
    xlim([1 10]);
    set(axnow, 'XTick', [1 5 10],'FontName','Times New Roman','fontsize',14);
    ylim([1e-4 1]);
    set(axnow, 'YTick', [1e-4 1e-2 1e0]);
        
    print(fig,strcat('render/cdf-T',Tstr,'-compare.eps'),'-depsc');
elseif 2==option || 3==option % different density
    tmat2 = strcat('../data_cdf/',parastr,'.dat');
    data_tmat2 = gettmat(tmat2);
    xs = 1:20;
    useT = str2double(Tstr);
    useRho = [1e-4, 1e-3, 1e-2, 1e-1];
    %useRho = [1e-4, 2e-4, 5e-4, 5e-3];
    nrho = length(useRho);
    tm_check = data_tmat2(data_tmat2(:,1)==useT & sum(data_tmat2(:,2)==useRho,2), xs+2);
    ticklabel = 'x-';
    psize=2;unifyfigure;hold on;
    for ind=1:nrho
        plot(nan, nan, ticklabel);
    end
    set(gca,'ColorOrderIndex',1);
    xx=[1:0.1:10].';
    for ind=1:nrho
        filted = tm_check(ind,:);
    %    yys(:,ind)=pchip(xs, filted, xx);
        plot(xs, filted, ticklabel);
    end
    set(gca,'ColorOrderIndex',1);
    %legend('10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}');
    %legend boxoff;
    xlim([1 10]);ylim([0 1]);
    set(gca,'XTick',linspace(2,10,5));
    set(gca,'YTick',linspace(0,1,6));
    set(gca,'fontsize',16);
    xlabel('$n$', 'Interpreter', 'Latex');ylabel('$K(n)$', 'Interpreter', 'Latex');
    
    if 3==option
        curpos = get(gcf, 'Position');
        axnow=axes('Position', [0.4, 0.5, 0.4, 0.32]);
        for ind=1:nrho
            filted = tm_check(ind,:);
            plot(xs, filted, '-');
            if 1==ind
                hold on;
            end
        end
        box(axnow,'off');
        set(axnow, 'YScale' ,'log');
        xlim([1 20]);
        set(axnow, 'XTick', [1 10 20],'FontName','Times New Roman','fontsize',14);
        ylim([1e-10 1]);
        set(axnow, 'YTick', [1e-10 1e-5 1e0]);
    end
    print(fig,strcat('render/cdf-T',Tstr,'.eps'),'-depsc');
end

%% Get cdf results from data file of transfer-matrix calculation
function result=gettmat(filename)
    fin = fopen(filename, 'r');
    str = fgetl(fin);  % discard the first line
    % k = strfind(str, '=');
    status = 0;
    result = [];
    while(~feof(fin))
        str = fgetl(fin);
        if(status==0)
            Trhop = strread(str);  % line of T, rho, p
            str = fgetl(fin);
            cds = strread(str);% line of cluster distribution function
            result = [Trhop([1 2]), cds];
            status=1;
        elseif(status==1)
            Trhop = strread(str);
            status=2;
        else
            cds = strread(str);
            result = [result; Trhop([1 2]), cds];
            status=1;
        end
    end
    fclose(fin);
end

%% Get cdf results from data file of simulations
function result=getsimu(filename)
    fin = fopen(filename, 'r');
    status = 0;
    result = [];
    while(~feof(fin))
        str = fgetl(fin);
        if(status==0)
            Trhop = strread(str);  % line of T, rho, p
            str = fgetl(fin);
            cds = strread(str);% line of cluster distribution function
            result = [Trhop, cds];
            status=1;
        elseif(status==1)
            Trhop = strread(str);
            status=2;
        else
            cds = strread(str);
            result = [result; Trhop, cds];
            status=1;
        end
    end
    fclose(fin);
end
