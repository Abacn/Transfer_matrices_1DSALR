% cdfsearcher
% find the density where a minimum of k(2) happens
for xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 2, 4]
% for xis=[, 4]
fclose('all');
writefile = false;
coeffs = [1, 2.5, 4, 1, xis];
Ts = [ 0.15,0.2:0.02:0.58,0.6:0.1:1]; %
%Ts = [0.32:0.02:0.38];
tol = 1e-2;
bplow = 1e-5;
bphigh = 1;
k1lowlimit = 1e-2;
foldername = 'data_cdf';
nm = sprintf('t_%.1f_%.1f_%.1f_%.1f_%.2f',coeffs(1), coeffs(2), coeffs(3), ...
             coeffs(4), coeffs(5)); % t means transition
         
%Ts=linspace(0.1,1,19);
if(writefile)
% Pick up a filename
suffix='.dat';

% Scan if file exist
flag = true;
count = 0;
if ~exist(foldername, 'dir')
mkdir(foldername);
end
while flag
if count==0
fn = strcat(foldername,'/',nm,suffix);
else
fn = sprintf('%s/%s_%d%s',foldername,nm,count,suffix);
end
if exist(fn, 'file')
count = count + 1;
else
flag = false;
end
end
fprintf('Saveto: %s\n', fn);
FilePtr = fopen(fn, 'w');
else
   FilePtr = 1;
end

logfile = strcat(foldername,'/cdfsearcher.log');
LogPtr = fopen(logfile, 'a');

fprintf(LogPtr, '-- %s --\np\tk(1)\tk(2)\tk(3)\n', datestr(clock));
fprintf(FilePtr, 'T\tp\trho\ttype\n'); % type: 2: k(2)>k(3) or 3: k(2) ,= k(3)
for T=Ts
    beta = 1/T;
    plow = bplow*T;
    phigh = bphigh*T;
    % first test the low limit
    lowflag = 1;
    fprintf(LogPtr,'T=%.3f\n', T);
    while(lowflag)
        [klow, ~] = cdf(plow, beta, coeffs, 200, 4);
        fprintf(LogPtr, '%.3e\t', [plow, klow(:).']); fprintf(LogPtr, '\n');
        if not (klow(1) > klow(2) && klow(2) > klow(3) && klow(3) > klow(4))
            phigh = plow;
            khigh = klow;
            plow = plow / 10;
            lowflag = 2; % default pmin changed
        else
            break;
        end
    end
    if 1==lowflag % change phigh
        while(phigh > plow)
            phightmp = phigh/10;
            [khightmp, ~] = cdf(phightmp, beta, coeffs, 200, 4);
            fprintf(LogPtr, '%.3e\t', [phightmp, khightmp(:).']); fprintf(LogPtr, '\n');
            if not (khightmp(1) > khightmp(2) && khightmp(2) > khightmp(3) && khightmp(3) > khightmp(4))
                phigh = phightmp;
            else
                if phightmp<plow
                    plow = phightmp;
                end
                break;
            end
        end
    end
    %{
    if 2==lowflag
        highflag = 2;  % check if phigh is too large
    else
        [khigh, ~] = cdf(phigh, beta, coeffs, 200, 4);
        fprintf(LogPtr, '%.3e\t', [phigh, khigh(:).']); fprintf(LogPtr, '\n');
        if khigh(1) < k1lowlimit
            highflag = 2;
        elseif khigh(1) > khigh(2) && khigh(2) > khigh(3)  % phigh is too low
            highflag = 3;
        else
            highflag = 0;
        end
    end
    if 2==highflag
        while (highflag)
            phighprev = phigh;
            phigh = (phigh + plow) / 2;
            [khigh, ~] = cdf(phigh, beta, coeffs, 200, 4);
            fprintf(LogPtr, '%.3e\t', [phigh, khigh(:).']); fprintf(LogPtr, '\n');
            if khigh(1) > k1lowlimit
                % plow = phigh;
                phigh = phighprev;
                break;
            end
        end
    elseif 3==highflag
        while (highflag)
            phigh = phigh*2;
            [khigh, ~] = cdf(phigh, beta, coeffs, 200, 4);
            fprintf(LogPtr, '%.3e\t', [phigh, khigh(:).']); fprintf(LogPtr, '\n');
            if ~(khigh(1) > khigh(2) && khigh(2) > khigh(3))
                break;
            end
        end
    end
    %}
    plowlimit = plow;
    phighlimit = phigh;
    epsp = inf;
    while(epsp > tol)
        pcen = (plow+phigh)*0.5;
        [ks, ~] = cdf(pcen, beta, coeffs, 200, 4);
        if ks(1) > ks(2) && ks(2) > ks(3)  && ks(3) >ks(4)
            plow = pcen;
        else
            phigh = pcen;
            khigh = ks;
        end
        fprintf(LogPtr, '%.3e\t', [pcen, ks(:).']); fprintf(LogPtr, '\n');
        epsp = (phigh-plow)/pcen;
    end 
    rho = findrho(pcen, beta, coeffs);
    if khigh(1) < khigh(2) % transition type
        ttype = 2;      % repulsion induced (k(1)<k(2))
    else
        if khigh(2) < khigh(3)
            if khigh(3) < khigh(4)
                ttype = 4;  % repulsion induced (k(2)<k(3))
            else
                ttype = 3;  % cluster
            end
        else
            ttype = 1;  % condensation
        end
    end
    fprintf(FilePtr, '%.3f\t%.3e\t%.3e\t%d\n', T, pcen, rho, ttype);
end

if writefile
    fclose(FilePtr);
end
fclose(LogPtr);

end