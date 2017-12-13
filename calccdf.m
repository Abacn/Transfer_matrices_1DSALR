% calccdf.m
% Calculate cluster distribution function
writefile = false;
coeffs=[1,2.5,4,1,1];
%Ts=linspace(0.1,1,10);
Ts=[0.25];
N = 20;

option = 'rho'; % 'p' or 'rho'
if 'rho'==option
    rhos = [1e-2,5e-2,1e-1,5e-1,linspace(0.2,0.8,4)];
else
    ps = [2e-2,3e-2];
end
foldername = 'data_cdf';
nm = sprintf('%.1f_%.1f_%.1f_%.1f_%.2f',coeffs(1), coeffs(2), coeffs(3), ...
             coeffs(4), coeffs(5));
         
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
rq=1;
fprintf(FilePtr, '%s\t%s\t%s\tN=%d\n', 'T', 'rho', 'p', N);
for T=Ts
  rp=1;
  beta = 1/T;
  if 'rho'==option
      for rho=rhos
        p = findp(rho, T, coeffs, 1e-5);
        %p = rho*T;
        fprintf(FilePtr, '%.3f\t%.3e\t%.3e\n', T, rho, p);
        [ks, ~] = cdf(p, beta, coeffs, 200, N);
        fprintf(FilePtr, '%.3e\t', ks);
        fprintf(FilePtr, '\n');
        rp=rp+1;
      end
  else
      for p=ps
        rho = findrho(p, 1/T, coeffs);
        %p = rho*T;
        fprintf(FilePtr, '%.3f\t%.3e\t%.3e\n', T, rho, p);
        [ks, ~] = cdf(p, beta, coeffs, 200, N);
        fprintf(FilePtr, '%.3e\t', ks);
        fprintf(FilePtr, '\n');
        rp=rp+1;
      end
  end
  rq=rq+1;
end
