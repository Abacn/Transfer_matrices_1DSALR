% minsearcher.m
% find minimum of (\beta p - rho)/rho^2 vs. rho / p
function minsearcher()
%% Script
close all;
coeffs=[1,2.5,4,1,0.05];
%Ts = [0.32];
Ts = linspace(0.1, 0.4, 31);
betas=1./Ts;
Lowlimit = -8;
NMag = 20;
NDivide = 10;
Accu = 1e-2;

% Pick up a filename
suffix='.dat';
nm = sprintf('%.1f_%.1f_%.1f_%.1f_%.2f',coeffs(1), coeffs(2), coeffs(3), ...
    coeffs(4), coeffs(5));
% Scan if file exist
flag = true;
count = 0;
if ~exist('data_min', 'dir')
    mkdir('data_min');
end
while flag
    if count==0
        fn = strcat('data_min/',nm,suffix);
    else
        fn = sprintf('data_min/%s_%d%s',nm,count,suffix);
    end
    if exist(fn, 'file')
        count = count + 1;
    else
        flag = false;
    end
end
fprintf('Saveto: %s\n', fn);
FilePtr = fopen(fn, 'w');
records = [];
rho = zeros(1, NMag);
for beta=betas
    flag=1;
    rp=1;
    ps=10.^linspace(Lowlimit, -1, NMag);
    for p=ps/beta
        rho(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p);
        rp = rp + 1;
    end
    y = (ps./rho-1)./rho;
    [val, ind] = min(y);
    if(ind==1 || ind==NMag)
        flag = false;   % minimum not found
        % fprintf(FilePtr, 'T=%.3f: not found\n', 1/beta);
        raw_record(FilePtr, 1/beta, ps, rho, y);
        continue;
    else
        % there must be a minimum within the data range
        raw_record(FilePtr, 1/beta, ps, rho, y);
        prange = [ps(ind-1), ps(ind+1)];
        ps = linspace(prange(1), prange(2), NDivide);
    end
    rho = zeros(1, NDivide);
    while(flag)
        rp=1;
        for p=ps/beta
            rho(rp) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p);
            rp = rp + 1;
        end
        y = (ps./rho-1)./rho;
        ind = min2(y);
        if(ind==1)
            prange = [ps(1), ps(2)];
        elseif(ind==NDivide)
            prange = [ps(NDivide-1), ps(NDivide)];
        else
            prange = [ps(ind-1), ps(ind+1)];
        end
        err = (prange(2)-prange(1))/2;
        if(err/rho(ind)<Accu)
            % accuracy demanding reached
            records(:,end+1) = [1/beta, ps(ind), rho(ind), y(ind), err];
            raw_record(FilePtr, 1/beta, ps, rho, y);
            flag = false;
        else
            % continue iteration
            raw_record(FilePtr, 1/beta, ps, rho, y);
            ps = linspace(prange(1), prange(2), NDivide);
        end
    end
end
min_record(FilePtr, records);
fclose(FilePtr);
end

%% Functions
% record the raw data
function raw_record(FilePtr, T, ps, rho, y)
rp = 1;
for pbeta=ps
    fprintf(FilePtr, '%.3f\t%.2e\t%.2e\t%.2e\n', T, ps(rp), rho(rp), y(rp));
    rp = rp + 1;
end
end

% record the minimuns
function min_record(FilePtr, records) % [Ts, ps, rho, y, err]
rp = 1;
fprintf(FilePtr, '%s\t%s\t%s\t%s\t%s\n', 'T', 'beta*p', 'rho', 'y', 'err of rho');
fprintf(FilePtr, '%.3f\t%.2e\t%.2e\t%.2e\t%.2e\n', records);
end

function ind = min2(vec)
% find local minimum from back
sz = size(vec(:),1);
while(sz>1)
    if(vec(sz-1)>vec(sz))
        break;
    end
    sz = sz - 1;
end
ind = sz;
return;
end
