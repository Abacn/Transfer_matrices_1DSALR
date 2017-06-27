% print secondB and thirdB to a file

coeffs = [1, 2.2, 4, 1];
xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5];
Ts = linspace(0.06, 0.4, 18);

% Pick up a filename
suffix='.dat';
nm = sprintf('%.1f_%.1f_%.1f_%.1f',coeffs(1), coeffs(2), coeffs(3), ...
    coeffs(4));
% Scan if file exist
flag = true;
count = 0;
savedir = '../data_Bs'
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
records = [];

fprintf(FilePtr, 'T\tB_2\tB_3\n');
for xi=xis
    results = zeros(size(Ts,2), 3);
    coeff = [coeffs, xi];
    sb = secondB(Ts, coeff);
    tb = thirdB(Ts, coeff);
    fprintf(FilePtr, 'Xi=%.3f\n', xi);
    for rp=1:length(Ts)
        fprintf(FilePtr, '%.3f\t%.3e\t%.3e\n', Ts(rp), sb(rp), tb(rp));
    end
end
fclose(FilePtr);