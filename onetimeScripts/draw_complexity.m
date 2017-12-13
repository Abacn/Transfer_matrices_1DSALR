% draw_complexity
% calculate the relaxation time (step) of different algorithms

% 1 - with respect to N, both eqdis and compact
% 2 - with respect to T
% 3 - draw raw data
% 4 - same as 1 but drow compact and eqdis in same plot
close all;
addpath('..');
coeffs = [1,2.5,4,1,1];
marks = {'*-b';'*:b';'o-r';'o:r';'x-k'};
algonames = {'MMC'; 'CMMC'; 'HMC'; 'CHMC'; 'ECMC'};
option = 4;
mode = 'compact'; % 'eqdis' or 'compact'

lsqopts = optimset('Display','off');

if 1==option    % with respect to N
    Tstr = '0.22';
    rhostr = '0.001';
    calcp = false;  % switch for debug. p calucation is time-consuming
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/T',Tstr,',D',rhostr,'/gdf/');
    Tdb = str2double(Tstr);
    Temperature = Tdb;
    rhodb = str2double(rhostr);
    savestr = strcat('onetimeData/draw_complexity-T',Tstr,',D',rhostr,'.mat');
    if or(calcp, ~exist(savestr, 'file'))
        [Pdb,  PHI] = getPPHI(Tdb, rhodb, coeffs);
        calcp = true;
    else
        load(savestr);
    end
    if strcmp('compact',mode)
        str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 's_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 's_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 's_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 's_chmc', str_suffix);
    elseif strcmp('eqdis',mode)
        str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 'eq_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 'eq_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 'eq_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 'eq_chmc', str_suffix);
    else
        error('Error: Mode not valid. (compact or eqdis)');
    end
    Cpath = {str_mmc, str_cmmc, str_hmc, str_chmc, str_ecmc};
    lenCpath = length(Cpath);
    alldata = {};

    fprintf('N\ttau/sampleLength\tPHI_end/PHI_eq\n');
    for rq=1:lenCpath
        foldername = Cpath{rq};
        fnames = strsplit(ls(foldername));
        dlst = sort(fnames(~cellfun(@isempty,regexp(fnames,'^[0-9]+\.dat$'))));
        record = zeros(size(dlst,2),5);  % N tau err_tau
        rp = 1;
        
        disp(strcat(algonames{rq},':'));
        
        for dl= dlst
            dl = dl{1};
            sl = strsplit(dl,'.');
            Npart = str2double(sl{1});
            [record(rp,:), ~ ] = getfit(strcat(foldername,dl), Npart, PHI, mode, Npart);
            rp = rp + 1;
        end
        record = sortrows(record);
        alldata{rq} = record;
    end
    % draw
    psize=2;unifyfigure;hold on;
    for rp=1:lenCpath
        thisdata = alldata{rp};
        plot(thisdata(:,1), thisdata(:,2), marks{rp});
    end
    set(gca,'XSCale', 'log');
    set(gca,'YSCale', 'log');
    if strcmp('compact',mode)
        xlim([90 10000]);
        ylim([1e2, 1e10]);
        set(gca, 'YTick', [1e2 1e6 1e10]);
        %legend(algonames, 'Location','northwest', 'Box','off');
    elseif strcmp('eqdis',mode)
        if Temperature == 0.5
            xlim([90 10000]);
            ylim([1, 500]);
            set(gca, 'YTick', [1 10 100]);
        else
            xlim([90 10000]);
            ylim([1e2, 1e6]);
            set(gca, 'YTick', [2e2, 1e4, 1e6]);
        end
    end
    xlabel('N');ylabel('\tau_{mix}');
    
    if(calcp)
        save(savestr, 'alldata', 'Cpath', 'Tdb', 'rhodb', 'Pdb', 'PHI');
    end
    unifyfigureend;
    if strcmp('compact',mode)
        print(fig,strcat('render/relax-T',Tstr,'.eps'),'-depsc');
    elseif strcmp('eqdis',mode)
        print(fig,strcat('render/relax-eqdis-T',Tstr,'.eps'),'-depsc');
    end
elseif 2==option  % with respect to T
    Nstr = '1000';
    N2str = '10000';
    rhostr = '0.001';
    calcp = false;  % switch for debug. p calucation is time-consuming
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/N',Nstr,',D',rhostr,'/gdf/');
    savestr = strcat('onetimeData/draw_complexity-Tvariate,D',rhostr,'.mat');
    NpartA = str2double(Nstr);
    NpartB = str2double(N2str);
    
    rhodb = str2double(rhostr);
    if or(calcp, ~exist(savestr, 'file'))
        thermos = []; % T, P, PHI
    else
        load(savestr);
    end
    if strcmp('compact',mode)
        str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 's_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 's_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 's_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 's_chmc', str_suffix);
    elseif strcmp('eqdis',mode)
        str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 'eq_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 'eq_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 'eq_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 'eq_chmc', str_suffix);
    else
        error('Error: Mode not valid. (compact or eqdis)');
    end

    Cpath = {str_mmc, str_cmmc, str_hmc, str_chmc, str_ecmc};
    lenCpath = length(Cpath);
    alldata = {[],[],[]};
    fprintf('N\ttau/sampleLength\tPHI_end/PHI_eq\n');
    for rq=1:lenCpath
        foldername = Cpath{rq};
        fnames = strsplit(ls(foldername));
        dlst = sort(fnames(~cellfun(@isempty,regexp(fnames,'^[0-9\.]+\.dat$'))));
        record = zeros(size(dlst,2),5);  % N tau err_tau
        
        disp(strcat(algonames{rq},':'));
        
        rp = 1;
        for dl= dlst
            dl = dl{1};
            Temperature = str2double(dl(1:end-4));
            
            if isempty(thermos)
                [Pdb,  PHI] = getPPHI(Temperature, rhodb, coeffs);
                thermos= [thermos; Temperature, Pdb, PHI];
            else
                tindex = find(thermos(:,1)==Temperature);
                if tindex
                    Pdb = thermos(tindex, 2);
                    PHI = thermos(tindex, 3);
                else
                    [Pdb,  PHI] = getPPHI(Temperature, rhodb, coeffs);
                    thermos= [thermos; Temperature, Pdb, PHI];
                end
            end
            nowfilename = strcat(foldername, dl);
            upgfilename = strrep(nowfilename, Nstr, N2str);
            if exist(upgfilename, 'file')
                Npart = NpartB;
                nowfilename = upgfilename;
            else
                Npart = NpartA;
            end
            [record(rp,:), ~ ] = getfit(nowfilename, Npart, PHI, mode, Temperature);
            rp = rp + 1;
        end
        record = sortrows(record);
        alldata{rq} = record;
    end
    % draw
    psize=2;unifyfigure;hold on;
    for rp=1:lenCpath
        thisdata = alldata{rp};
        plot(thisdata(:,1), thisdata(:,2), marks{rp});
    end
    xlim([0.18 0.5]);
    set(gca,'YSCale', 'log');
    set(gca, 'XTick', linspace(0.2,0.5,4));
    if strcmp('compact',mode)
        tickpos = round( logspace(log10(1e0),log10(1e9), 4) );
        ylimarr = [1e3 1e9];
    elseif strcmp('eqdis',mode)
        tickpos = round( logspace(log10(1e0),log10(1e6), 3) );
        ylimarr = [1e0 1e7];
    end
    ylim(ylimarr);
    plot([0.26 0.26], ylimarr, ':k', 'LineWidth', 1);
    set(gca, 'YTick', tickpos);
    xlabel('$T$','Interpreter','Latex');ylabel('\tau_{mix}');
    
    % legend(algonames, 'Location','northeast', 'Box','off');
    save(savestr, 'alldata', 'thermos');
    unifyfigureend;
    if strcmp('compact',mode)
        print(fig,strcat('render/relax-T,D',rhostr,'.eps'),'-depsc');
    elseif strcmp('eqdis',mode)
        print(fig,strcat('render/relax-eqdis-T,D',rhostr,'.eps'),'-depsc');
    end
elseif 3==option  % draw example
    Nstr = '1000';
    N2str = '10000';
    rhostr = '0.001';
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/N',Nstr,',D',rhostr,'/gdf/');
    savestr = strcat('onetimeData/draw_complexity-Tvariate,D',rhostr,'.mat');
    NpartA = str2double(Nstr);
    NpartB = str2double(N2str);
    
    Tstr = '0.22';
    rhostr = '0.001';
    xlimup = 1e6;
    Temperature = str2double(Tstr);
    calcp = false;  % switch for debug. p calucation is time-consuming
    rhodb = str2double(rhostr);
    if or(calcp, ~exist(savestr, 'file'))
        thermos = []; % T, P, PHI
    else
        load(savestr);
    end
    if strcmp('compact',mode)
        str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 's_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 's_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 's_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 's_chmc', str_suffix);
    elseif strcmp('eqdis',mode)
        str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
        str_mmc = strcat(str_prefix, 'eq_mmc', str_suffix);
        str_cmmc = strcat(str_prefix, 'eq_cmmc', str_suffix);
        str_hmc = strcat(str_prefix, 'eq_hmc', str_suffix);
        str_chmc = strcat(str_prefix, 'eq_chmc', str_suffix);
    else
        error('Error: Mode not valid. (compact or eqdis)');
    end

    Cpath = {str_mmc, str_cmmc, str_hmc, str_chmc, str_ecmc};
    lenCpath = length(Cpath);
    alldata = {[],[],[]};
    fprintf('Algo\ttau/sampleLength\tPHI_end/PHI_eq\n');
    psize=2;unifyfigure;hold on;
    for rq=1:lenCpath
        plot(nan, nan);
    end
    for rq=lenCpath:-1:1
        foldername = Cpath{rq};
        fnames = strsplit(ls(foldername));
        dlst = fnames(~cellfun(@isempty,regexp(fnames,'^[0-9\.]+\.dat$')));
        record = zeros(size(dlst,2),5);  % N tau err_tau
        
        disp(strcat(algonames{rq},':'));
        
        findflag = false;
        for dl= dlst
            dl = dl{1};
            if Temperature == str2double(dl(1:end-4))
                findflag = true;
                break;
            end
        end  
        if ~findflag
            error(strcat('Error: ',algonames{rq} , ' Not found'));
        end
        if isempty(thermos)
            [Pdb,  PHI] = getPPHI(Temperature, rhodb, coeffs);
            thermos= [thermos; Temperature, Pdb, PHI];
        else
            tindex = find(thermos(:,1)==Temperature);
            if tindex
                Pdb = thermos(tindex, 2);
                PHI = thermos(tindex, 3);
            else
                [Pdb,  PHI] = getPPHI(Temperature, rhodb, coeffs);
                thermos= [thermos; Temperature, Pdb, PHI];
            end
        end
        nowfilename = strcat(foldername, dl);
        upgfilename = strrep(nowfilename, Nstr, N2str);
        if exist(upgfilename, 'file')
            Npart = NpartB;
            nowfilename = upgfilename;
        else
            Npart = NpartA;
        end
        if strcmp(mode, 'compact')
            T0 = (Npart-1)/Npart;
        end
        
        data=dlmread(nowfilename,'',1,0);
        data = sortrows(data);
            if strcmp('compact',mode)
                    %func=@(c,x)(PHI+(T0-PHI)*(c(3)*exp(-x/c(1))+(1-c(3))*exp(-x/c(2))));fittype=3;
                    %func=@(c,x)(PHI+(T0-PHI)*exp(-(x/c(1)).^c(2)));fittype = 2;
                    func=@(c,x)(PHI+(T0-PHI)*exp(-x/c));fittype = 1;
            elseif strcmp('eqdis',mode)
                    %func=@(c,x)(PHI*(1-c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype = 3;
                    %func=@(c,x)(PHI*(1-exp(-(x/c(1).^c(2)))));fittype = 2;
                    func=@(c,x)(PHI*(1-exp(-x/c)));fittype = 1;
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            iend = size(data, 1);
            xrescale = data(end,1)/10;
            cinitial = data(floor(iend/2)+1, 1)/xrescale;          
            if fittype == 1
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitial,data(:,1)/xrescale,data(:,2), [],[], lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = fitc*xrescale / Npart;
            elseif 2==fittype
                cinitials = [cinitial, 1];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2), [],[],lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = fitc(1)*xrescale / Npart;
                record(rp,:) = [Npart, fittau,nan, fitc(2), max(data(:,1))];
            elseif fittype == 3
                cinitials = [cinitial, cinitial/100, 0.9];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2), [],[], lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = max(fitc(1), fitc(2))*xrescale / Npart;
                % TODO, see how to estimate the error of fitted parameter
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            vldchk = sum(data(end-100:end , 2))/100/PHI;
            fprintf('%s\t%.4f\n', algonames{rq},  vldchk);
            set(gca,'ColorOrderIndex',rq);
            rawdata = sortrows(data);
            if strcmp(foldername, 's_ecmc')
                smthdata = smooth(rawdata(:,1), rawdata(:,2), 3);
            else
                smthdata = smooth(rawdata(:,1), rawdata(:,2), 20);
            end
            if strcmp('compact',mode) && ~strcmp('s_ecmc',foldername)
                ind = 1:100:size(data,1);
            else
                ind = 1:20:size(data,1);
            end
            plot(data(ind,1)/Npart, smthdata(ind), marks{rq}(2:end));
            if 1==rq
                hold on;
            end
            %{
            if strcmp('compact',mode)
                x=10.^linspace(0, 6+log10(Npart));
            else
                x=10.^linspace(0, 3+log10(Npart));
            end
            y=arrayfun(@(xx)(func(fitc, xx/xrescale)), x);
            set(gca,'ColorOrderIndex',rq);
            plot(x/Npart, y);
            %}
    end
    % draw
    plot([0.1, 1e7],[PHI, PHI], ':k');
    if strcmp('compact',mode)
        ylimlow = 0;

        xlim([1 2e6]);
        ylim([ylimlow 1]);

        %set(gca,'XSCale', 'log');
        set(gca, 'XTick', [0 1e6 2e6]);
        set(gca, 'YTick', [0 0.5 1]);
        %legend(algonames, 'Location', 'Northeast');
    else % eqdis
        ylimlow = 0;
        %set(gca,'XSCale', 'log');
        xlim([0 2e6]);
        ylim([0 0.2]);
        set(gca, 'XTick', [0 1e6 2e6]);
        set(gca, 'YTick', [0 0.1 0.2]);      
        %legend(algonames, 'Location', 'Southeast');
    end
    legend boxoff;
    xlabel('$t_{sweep}$','Interpreter','Latex');ylabel('$C_{bind}$','Interpreter','Latex');
    unifyfigureend;
    if strcmp('compact',mode)
        print(fig,strcat('render/draw-T',Tstr,',D',rhostr,'.eps'),'-depsc');
    elseif strcmp('eqdis',mode)
        print(fig,strcat('render/draw-eqdis-T',Tstr,',D',rhostr,'.eps'),'-depsc');
    end
elseif 4==option
    Tstr = '0.22';
    rhostr = '0.001';
    calcp = false;  % switch for debug. p calucation is time-consuming
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/T',Tstr,',D',rhostr,'/gdf/');
    Tdb = str2double(Tstr);
    Temperature = Tdb;
    rhodb = str2double(rhostr);
    savestr = strcat('onetimeData/draw_complexity-T',Tstr,',D',rhostr,'.mat');
    if or(calcp, ~exist(savestr, 'file'))
        [Pdb,  PHI] = getPPHI(Tdb, rhodb, coeffs);
        calcp = true;
    else
        load(savestr);
    end
    
    str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
    str_mmc = strcat(str_prefix, 's_mmc', str_suffix);
    str_cmmc = strcat(str_prefix, 's_cmmc', str_suffix);
    str_hmc = strcat(str_prefix, 's_hmc', str_suffix);
    str_chmc = strcat(str_prefix, 's_chmc', str_suffix);
    Cpathcp = {str_mmc, str_cmmc, str_hmc, str_chmc, str_ecmc};
    
    str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
    str_mmc = strcat(str_prefix, 'eq_mmc', str_suffix);
    str_cmmc = strcat(str_prefix, 'eq_cmmc', str_suffix);
    str_hmc = strcat(str_prefix, 'eq_hmc', str_suffix);
    str_chmc = strcat(str_prefix, 'eq_chmc', str_suffix);
    Cpatheq = {str_mmc, str_cmmc, str_hmc, str_chmc, str_ecmc};
    
    lenCpath = length(Cpathcp);

    fprintf('N\ttau/sampleLength\tPHI_end/PHI_eq\n');
    for rr=1:2
        if 1==rr
            Cpath = Cpathcp;
            mode = 'compact';
        else
            Cpath = Cpatheq;
            mode = 'eqdis';
        end
        alldata = {};
        for rq=1:lenCpath
            foldername = Cpath{rq};
            fnames = strsplit(ls(foldername));
            dlst = sort(fnames(~cellfun(@isempty,regexp(fnames,'^[0-9]+\.dat$'))));
            record = zeros(size(dlst,2),5);  % N tau err_tau
            rp = 1;

            disp(strcat(algonames{rq},':'));

            for dl= dlst
                dl = dl{1};
                sl = strsplit(dl,'.');
                Npart = str2double(sl{1});
                [record(rp,:), ~ ] = getfit(strcat(foldername,dl), Npart, PHI, mode, Npart);
                rp = rp + 1;
            end
            record = sortrows(record);
            alldata{rq} = record;
        end
        if 1==rr
            alldatacp = alldata;
        else
            alldataeq = alldata;
        end
    end
    % draw
    psize=21;unifyfigure;hold on;
    for rp=1:lenCpath
        thisdata = alldatacp{rp};
        plot(thisdata(:,1), thisdata(:,2), marks{rp});
    end
    for rp=1:lenCpath
        thisdata = alldataeq{rp};
        plot(thisdata(:,1), thisdata(:,2), marks{rp});
    end
    set(gca,'XSCale', 'log');
    set(gca,'YSCale', 'log');
    if Temperature == 0.5
        xlim([90 10000]);
        ylim([1e0, 1e9]);
        set(gca, 'YTick', [1 1e3 1e6 1e9]);
        %legend(algonames, 'Location','northwest', 'Box','off');
    else
        
        xlim([90 10000]);
        ylim([1e0, 1e9]);
        set(gca, 'YTick', [1, 1e3, 1e6, 1e9]);
    end
    xlabel('N');ylabel('\tau_{mix}');
    
    if(calcp)
        save(savestr, 'alldata', 'Cpath', 'Tdb', 'rhodb', 'Pdb', 'PHI');
    end
    unifyfigureend;
    print(fig,strcat('render/relax-T',Tstr,',D',rhostr,'.eps'),'-depsc');
end

function [P, PHI] = getPPHI(T, rho, coeffs)
    P = findp(rho, T, coeffs, 1e-3);
    [rlist, pros] = gapdf(P, 1/T, coeffs);
    PHI = sum(pros(rlist < coeffs(1)*coeffs(2)))*(rlist(2)-rlist(1));
end

% read from file and fit
% fname: readin filename; Npart: number of particle; 
% mode: 'compact' or 'eqdis'; record_value: identifier that return in ret_record 
function [ret_record, ret_data] = getfit(fname, Npart, PHI, mode, record_value)
            lsqopts = optimset('Display','off');
            T0 = (Npart-1)/Npart;
            data=dlmread(fname,'',1,0);
            data = sortrows(data);
            if strcmp('compact',mode)
                    %func=@(c,x)(PHI+(T0-PHI)*(c(3)*exp(-x/c(1))+(1-c(3))*exp(-x/c(2))));fittype=3;
                    func=@(c,x)(PHI+c(2)*exp(-x/c(1)));fittype = 2;
                    %func=@(c,x)(PHI+(T0-PHI)*exp(-x/c));fittype = 1;
            elseif strcmp('eqdis',mode)
                    %func=@(c,x)(PHI*(1-c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype = 3;
                    func=@(c,x)(PHI*(1-c(2)*exp(-x/c(1))));fittype = 2;
                    %func=@(c,x)(PHI*(1-exp(-x/c)));fittype = 1;
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            ret_data = [data(:,1)/Npart, data(:,2)];
            data = data(floor(size(data,1)/100):end, :);
            iend = size(data, 1);
            xrescale = data(end,1)/10;
            cinitial = data(floor(iend/2)+1, 1)/xrescale;          
            if fittype == 1
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitial,data(:,1)/xrescale,data(:,2), [],[], lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = fitc*xrescale / Npart;
                ret_record = [record_value, fittau, fittau-errc(1), errc(2)-fittau, max(data(:,1))];
            elseif 2==fittype
                cinitials = [cinitial, 1];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2), [],[],lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = fitc(1)*xrescale / Npart;
                ret_record = [record_value, fittau,nan, fitc(2), max(data(:,1))];
            elseif fittype == 3
                cinitials = [cinitial, cinitial/100, 0.9];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2), [],[], lsqopts);
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale / Npart;
                fittau = max(fitc(1), fitc(2))*xrescale / Npart;
                ret_record = [record_value, fittau, nan, nan, max(data(:,1))];
                % TODO, see how to estimate the error of fitted parameter
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            vldchk = sum(data(end-100:end , 2))/100/PHI;
            if record_value>10
                fprintf('%d\t%.3f\t%.4f\n', record_value, ret_record(2)*Npart/ret_record(5), vldchk);
            else
                fprintf('%.3f\t%.3f\t%.4f\n', record_value, ret_record(2)*Npart/ret_record(5), vldchk);
            end
end