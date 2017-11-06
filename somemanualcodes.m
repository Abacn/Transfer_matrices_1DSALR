part = 7;
if (1==part)
%% Set coeffs
coeffs = [1, 2.2, 4, 1, 1];
T = 2;
beta = 1/T;
ps = 1e-4;
p = ps/beta;

elseif (2==part)
%% Calculate density
rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 150))),p)

elseif (3==part)
%% plot, check
T = 0.16;
dcind = find(dc(:,1) == T);
dsind = find(ds(:,1) == T);
figure;hold on;xlabel('\beta p');ylabel('\rho');
plot(dc(dcind,2),dc(dcind,3));
plot(ds(dsind,3)/T,ds(dsind,2),'ob');

elseif(4==part)
%% Save variables
save('onetimeData/drawdiff.mat', 'rho', 'rho_NNN', 'rho_NN', 'rho_HS', 'ps', 'coeffs');

%% Plot entropy
elseif(5==part)
   
%% Compute cdf
elseif(6==part)
    coeffs = [1, 2.5, 4, 1, 1];
    T = 0.1;
    beta = 1/T;
    rho = 1e-3;
    p = findp(rho, T, coeffs)
    % p = 8.335e-06;
    N = 10;
    [ks, ~] = cdf(p, beta, coeffs, 200, N)
    % plot(data(:,2),data(:,3),data(:,2),ks);
    % legend('tmat','ecmc');xlabel('N');ylabel('k(n)');title('coeffs=[1,2.5,4,1,1],T=0.5,\rho=0.2');
    % 
%% See rho varying with transfer matrix size
elseif(7==part)
    close all;
    coeffs = [1, 2.5, 4, 1, 1];
    P = 1e-2;
    T = 0.2;
    xs = 0:5;
    start = 297;
    ys = [];
    Zs = [];
    for x=xs
        Zs(x+1) = Pfunc_isobaric_3NN(P, 1/T, coeffs, start+x);
        ys(x+1) = findrho(P, 1/T, coeffs, start+x);
    end
    figure;plot(xs+start, ys);
    figure;plot(xs+start, Zs);
%% draw S(k), data from g(r)
elseif(8==part)
    Tstr = '0.2';
    foldername = strcat('../simulationdata/vmmc/outputs/1,2.5,4,1,1/T',Tstr,'/gofr/');
    fnames = strsplit(ls(foldername));
    dlst = fnames(~cellfun(@isempty,regexp(fnames,'^\d.+\.dat')));
    record = zeros(size(dlst,1), 3);
    rp = 1;
    for dl= dlst
        dl = dl{1};
        grs=dlmread(strcat(foldername, dl),'',1,0);
        rho = str2double(dl(1:end-4));
        result=calcsq(rho, grs);
        [pks,locs] = findpeaks(result(:,2),'NPeaks',1,'SortStr','descend');
        record(rp,:) =  [rho, pks,result(locs,1)];
        rp = rp + 1;
    end
    record = sortrows(record);
    idx = abs(record(:,3)-5.1836) < 1e-4;
    plot2(record(idx,[1 2]), '-x');
    xlabel('\rho');ylabel('S(k_m)');
    title(strcat('S(k) at 5.2. T=',Tstr));
%% draw peak s(k) value
elseif(9==part)
    rhostr = '0.1';
    rho = str2double(rhostr);
    foldername = '../simulationdata/vmmc/outputs/1,2.5,4,1,1/';
    fnames = strsplit(ls(foldername));
    dlst = fnames(~cellfun(@isempty,regexp(fnames,'^T.+')));
    record = zeros(size(dlst,2),3);
    rp = 1;
    for dl= dlst
        dl = dl{1};
        T = str2double(dl(2:end));
        grs=dlmread(strcat(foldername, dl,'/gofr/', rhostr, '.dat'),'',1,0);
        result=calcsq(rho, grs);
        [pks,locs] = findpeaks(result(:,2),'NPeaks',1,'SortStr','descend');
        record(rp,:) =  [T, pks,result(locs,1)];
        rp = rp + 1;
    end
    record = sortrows(record);
    idx = abs(record(:,3)-5.2360) < 1e-4;
    plot(1 ./ (record(idx,1)), record(idx,2), '-x');
    xlabel('1/T');ylabel('S(k_m)');
%% Plot relaxation time
elseif(10==part)
    fname1 = '../simulationdata/d_mc_nvt/secondarydata/ssf/1,2.5,4,1,1.dat';
    fname2 = '../simulationdata/d_vmmc/secondarydata/ssf/1,2.5,4,1,1.dat';
    data1=dlmread(fname1,'',1,0);
    data2=dlmread(fname2,'',1,0);
    data1(:,1) = 1./data1(:,1);
    data2(:,1) = 1./data2(:,1);
    srho='0.1';
    %sT='0.25';
    plotraw(2,str2double(srho),1,3,data1,data2,'-x');
    %plotraw(1,str2double(sT),2,3,data1,data2,'-x');
    set(gca, 'YScale', 'Log');
    xlabel('1/T');ylabel('\tau');legend('mc-nvt', 'vmmc');
    title(strcat('\rho=', srho));
    %xlabel('\rho');ylabel('\tau');legend('mc-nvt', 'vmmc');
    %title(strcat('\rho=', srho));
end