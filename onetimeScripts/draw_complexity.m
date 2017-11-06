% draw_complexity
% calculate the relaxation time (step) of different algorithms
close all;
addpath('..');
coeffs = [1,2.5,4,1,1];
marks = ['*-b';'x-r';'o-k'];
option = 2;
mode = 'compact'; % 'eqdis' or 'compact'
if 1==option    % with respect to N
    Tstr = '0.25';
    rhostr = '0.1';
    calcp = false;  % switch for debug. p calucation is time-consuming
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/T',Tstr,',D',rhostr,'/gdf/');
    Tdb = str2double(Tstr);
    rhodb = str2double(rhostr);
    savestr = strcat('onetimeData/draw_complexity-T',Tstr,'.mat');
    if or(calcp, ~exist(savestr, 'file'))
        [Pdb,  PHI] = getPPHI(Tdb, rhodb, coeffs);
        calcp = true;
    else
        load(savestr);
    end
    if strcmp('compact',mode)
        str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
        str_vmmc = strcat(str_prefix, 's_vmmc', str_suffix);
        str_mcnvt = strcat(str_prefix, 's_mc_nvt', str_suffix);
    elseif strcmp('eqdis',mode)
        str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
        str_vmmc = strcat(str_prefix, 'eq_vmmc', str_suffix);
        str_mcnvt = strcat(str_prefix, 'eq_mc_nvt', str_suffix);
    else
        error('Error: Mode not valid. (compact or eqdis)');
    end
    Cpath = {str_mcnvt, str_vmmc, str_ecmc};
    lenCpath = length(Cpath);
    alldata = {[],[],[]};


    for rq=1:lenCpath
        foldername = Cpath{rq};
        fnames = strsplit(ls(foldername));
        dlst = fnames(~cellfun(@isempty,regexp(fnames,'^[0-9]+\.dat$')));
        record = zeros(size(dlst,2),4);  % N tau err_tau
        rp = 1;
        for dl= dlst
            dl = dl{1};
            sl = strsplit(dl,'.');
            Npart = str2double(sl{1});
            T0 = (Npart-1)/Npart;
            data=dlmread(strcat(foldername, dl),'',1,0);
            data = sortrows(data);
            fittype = 1;
            if strcmp('compact',mode)
                %func=@(c,x)(PHI+(T0-PHI)*exp(-x/c));fittype = 1;
                func=@(c,x)(PHI+(T0-PHI)*(c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype=3;
            elseif strcmp('eqdis',mode)
                %func=@(c,x)(PHI*(1-exp(-x/c)));fittype = 1;
                func=@(c,x)(PHI*(1-c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype = 3;
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end

            iend = size(data, 1);
            xrescale = data(end,1)/10;
            if 1==fittype
            cinitial = data(floor(iend/2)+1, 1)/xrescale;
            [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitial,data(:,1)/xrescale,data(:,2));
            errc = nlparci(fitc,resid,'jacobian',J) .* xrescale;
            fitc = fitc*xrescale;
            record(rp,:) = [Npart, fitc, fitc-errc(1), errc(2)-fitc];
            elseif 3==fittype
                cinitials = [cinitial, cinitial/100, 0.9];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2));
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale;
                fittau = max(fitc(1), fitc(2))*xrescale;
                record(rp,:) = [Npart, fittau, nan, nan];
                % TODO, see how to estimate the error of fitted parameter
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            rp = rp + 1;
        end
        record = sortrows(record);
        alldata{rq} = record;
    end
    % draw
    psize=2;unifyfigure;hold on;
    for rp=1:lenCpath
        thisdata = alldata{rp};
        plot(thisdata(:,1), thisdata(:,2), marks(rp,:));
    end
    set(gca,'XSCale', 'log');
    set(gca,'YSCale', 'log');
    if strcmp('compact',mode)
        xlim([40 10000]);
        tickpos = round( logspace(log10(1e2),log10(1e10), 5) );
        set(gca, 'YTick', tickpos);
    elseif strcmp('eqdis',mode)
        xlim([400 100000]);
        tickpos = round( logspace(log10(1e3),log10(1e10), 5) );
        %set(gca, 'YTick', tickpos);
    end
    xlabel('N');ylabel('\tau_{mix}');
    legend({'MMC', 'CMMC', 'ECMC'}, 'Location','northwest', 'Box','off');
    if(calcp)
        save(savestr, 'alldata', 'Cpath', 'Tdb', 'rhodb', 'Pdb', 'PHI');
    end
    if strcmp('compact',mode)
        print(fig,strcat('render/relax-T',Tstr,'.eps'),'-depsc');
    elseif strcmp('eqdis',mode)
        print(fig,strcat('render/relax-eqdis-T',Tstr,'.eps'),'-depsc');
    end
elseif 2==option  % with respect to T
    Nstr = '1000';
    rhostr = '0.1';
    calcp = false;  % switch for debug. p calucation is time-consuming
    str_prefix = '../../simulationdata/';
    str_suffix = strcat('/outputs/1,2.5,4,1,1/N',Nstr,',D',rhostr,'/gdf/');
    savestr = strcat('onetimeData/draw_complexity-Tvariate.mat');
    Npart = str2double(Nstr);
    T0 = (Npart-1)/Npart;
    rhodb = str2double(rhostr);
    if or(calcp, ~exist(savestr, 'file'))
        thermos = []; % T, P, PHI
    else
        load(savestr);
    end
    if strcmp('compact',mode)
        str_ecmc = strcat(str_prefix, 's_ecmc', str_suffix);
        str_vmmc = strcat(str_prefix, 's_vmmc', str_suffix);
        str_mcnvt = strcat(str_prefix, 's_mc_nvt', str_suffix);
    elseif strcmp('eqdis',mode)
        str_ecmc = strcat(str_prefix, 'eq_ecmc', str_suffix);
        str_vmmc = strcat(str_prefix, 'eq_vmmc', str_suffix);
        str_mcnvt = strcat(str_prefix, 'eq_mc_nvt', str_suffix);
    else
        error('Error: Mode not valid. (compact or eqdis)');
    end

    Cpath = {str_mcnvt, str_vmmc, str_ecmc};
    lenCpath = length(Cpath);
    alldata = {[],[],[]};

    for rq=1:lenCpath
        foldername = Cpath{rq};
        fnames = strsplit(ls(foldername));
        dlst = fnames(~cellfun(@isempty,regexp(fnames,'^[0-9\.]+\.dat$')));
        record = zeros(size(dlst,2),4);  % N tau err_tau
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
            data=dlmread(strcat(foldername, dl),'',1,0);
            data = sortrows(data);
            if strcmp('compact',mode)
                func=@(c,x)(PHI+(T0-PHI)*exp(-x/c));fittype = 1;
                %func=@(c,x)(PHI+(T0-PHI)*(c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype=3;
            elseif strcmp('eqdis',mode)
                %func=@(c,x)(PHI*(1-exp(-x/c)));fittype = 1;
                func=@(c,x)(PHI*(1-c(3)*exp(-x/c(1))-(1-c(3))*exp(-x/c(2))));fittype = 3;
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            iend = size(data, 1);
            xrescale = data(end,1)/10;
            cinitial = data(floor(iend/2)+1, 1)/xrescale;          
            if fittype == 1
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitial,data(:,1)/xrescale,data(:,2));
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale;
                fitc = fitc*xrescale;
                record(rp,:) = [Temperature, fitc, fitc-errc(1), errc(2)-fitc];
            elseif fittype == 3
                cinitials = [cinitial, cinitial/100, 0.9];
                [fitc,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(func,cinitials,data(:,1)/xrescale,data(:,2));
                errc = nlparci(fitc,resid,'jacobian',J) .* xrescale;
                fittau = max(fitc(1), fitc(2))*xrescale;
                record(rp,:) = [Temperature, fittau, nan, nan];
                % TODO, see how to estimate the error of fitted parameter
            else
                error('Error: Mode not valid. (compact or eqdis)');
            end
            
            rp = rp + 1;
        end
        record = sortrows(record);
        alldata{rq} = record;
    end
    % draw
    psize=2;unifyfigure;hold on;
    for rp=1:lenCpath
        thisdata = alldata{rp};
        plot(thisdata(:,1), thisdata(:,2), marks(rp,:));
    end
    %set(gca,'XSCale', 'log');
    set(gca,'YSCale', 'log');
    set(gca, 'XTick', linspace(0.2,1,5));
    if strcmp('compact',mode)
        tickpos = round( logspace(log10(1e4),log10(1e10), 4) );
    elseif strcmp('eqdis',mode)
        tickpos = round( logspace(log10(1e4),log10(1e10), 4) );
    end
    ylim([1e4 1e11]);
    set(gca, 'YTick', tickpos);
    xlabel('$T$','Interpreter','Latex');ylabel('\tau_{mix}');
    legend({'MMC', 'CMMC', 'ECMC'}, 'Location','northeast', 'Box','off');
    save(savestr, 'alldata', 'thermos');
    if strcmp('compact',mode)
        print(fig,'render/relax-T.eps','-depsc');
    elseif strcmp('eqdis',mode)
        print(fig,'render/relax-eqdis-T.eps','-depsc');
    end
end

set(gca,'fontsize',16);

function [P, PHI] = getPPHI(T, rho, coeffs)
    P = findp(rho, T, coeffs, 1e-2);
    [rlist, pros] = gapdf(P, 1/T, coeffs);
    PHI = sum(pros(rlist < coeffs(1)*coeffs(2)))*(rlist(2)-rlist(1));
end