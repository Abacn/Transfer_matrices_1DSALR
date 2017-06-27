% search critical T according to the sign of B3
%xis=[0];
xis=[0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5].';
num_of_points = length(xis);
coeffs = [1,2.2,4,1,1];
coeffss=ones(num_of_points,1) * coeffs;
coeffss(:, 5) = xis;
Trange = [0.01, 0.5];
% Pick up a filename
suffix='.dat';
nm = sprintf('Tc_%.1f_%.1f_%.1f_%.1f',coeffs(1), coeffs(2), coeffs(3), ...
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

records = criticalT_thirdB(Trange, coeffss);

FilePtr = fopen(fn, 'w');
fprintf(FilePtr, 'Xi\tTc\n');
for rp=1:length(xis)
	fprintf(FilePtr, '%.3f\t%.4f\n', xis(rp), records(rp));
end

fclose(FilePtr);
