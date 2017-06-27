% draw p vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,2.5,4,1,2];
betas=[1 2 3 4 5];
Foldername = 'data_low';
suffix = '.mat';
% use pressure, not beta*pressure here
pr=[linspace(0.01,1,10), 1.5, linspace(2, 10, 9)];
%pr=linspace(0.1, 8, 20).';
%pr=10.^linspace(-4,-1,20).';
rq=1;
Nbeta = size(betas,2);
rho = zeros(size(pr,1),Nbeta); % NNN
for beta=betas
  rp=1;
  for p=pr
    rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 100))),p);
    rp=rp+1;
  end
  rq=rq+1;
end

y=zeros(size(rho));
for i=1:size(betas,2)
    y(:,i) = (pr(:)*betas(i)./rho(:,i)-1)./rho(:,i);
end

if(~exist(Foldername, 'dir'))
   mkdir(Foldername);
end
FilePtr = fopen(strcat(Foldername, '/outputs'),'a');
fprintf(FilePtr, '%s\n', datestr(now));
fprintf(FilePtr, '%s\n', num2str(coeffs,'%.1f '));

nm = num2str(coeffs,'%.1f_');
% Scan if file exist
flag = true;
count = 0;
while flag
    if count==0
        fn = strcat(Foldername,'/',nm(1:end-1),suffix);
    else
        fn = sprintf('%s/%s%d%s',Foldername, nm(1:end),count,suffix);
    end
    if exist(fn, 'file')
        count = count + 1;
    else
        flag = false;
    end
end
fprintf(FilePtr, 'Saveto: %s\n', fn);
fprintf(FilePtr, '=====\n');
fclose(FilePtr);
save(fn, 'betas', 'pr', 'coeffs', 'rho', 'y');
