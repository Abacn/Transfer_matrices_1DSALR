% print secondB and thirdB to a file
writeflag = true;
coeffs = [1, 2.5, 4, 1];
xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5];
Ts = 0.2:0.01:0.5;
if(writeflag)
% Pick up a filename
suffix='.dat';
nm = sprintf('B3_%.1f_%.1f_%.1f_%.1f',coeffs(1), coeffs(2), coeffs(3), ...
    coeffs(4));
% Scan if file exist
flag = true;
count = 0;
savedir = '../data_Bs';
if ~exist(savedir, 'dir')
    mkdir(savedir);
end
while flag
    if count==0
        fn = strcat(savedir,'/',nm,suffix);
    else
        fn = sprintf('%s/%s_%d%s',savedir,nm,count,suffix);
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

fprintf(FilePtr, 'Xi\tT\tB_3\n');
for xi=xis
    results = zeros(size(Ts,2), 3);
    coeff = [coeffs, xi];
    % sb = secondB(Ts, coeff);
    tb = thirdB(Ts, coeff);
    A=[];
    A(tb>0)=1;
    A(tb<0)=0;
    index = find(diff(A)~=0);
    records = [Ts.', tb.'];
    if(index)  % sign change exists
        Tlow = Ts(index);
        Thigh = Ts(index+1);
        Tol = Thigh-Tlow;
        while(Tol > 1e-4)
            Tcen = (Thigh+Tlow)/2;
            B3tmp = thirdB(Tcen, coeff);
            if(B3tmp>0)
              Thigh = Tcen;
            else
              Tlow=Tcen;
            end
            Tol = Thigh-Tlow;
            records(end+1,:) = [Tcen, B3tmp];
        end
        [~, si] = sort(records(:,1));
        records = records(si, :);
    end
    
    for rp=1:size(records, 1)
        fprintf(FilePtr, '%.3f\t%.3e\t%.3e\n', xi,records(rp,1), records(rp,2));
    end
end
if(writeflag)
fclose(FilePtr);
end