% draw (\beta p - rho)/rho^2 vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
close all;
coeffs=[1,2.5,4,1,2];
betas=[1];
suffix = '.mat';
%targetp = 2.0;ps = [-0.02 -0.01 0 0.01 0.02]+targetp;
%ps=linspace(0.1, 8, 20).';
ps=10.^linspace(-4,-1,20).';
rq=1;
Nbeta = size(betas,2);
rho = zeros(size(ps,1),Nbeta); % NNN
for beta=betas
  rp=1;
  for p=ps.'/beta
    rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p);
    rp=rp+1;
  end
  rq=rq+1;
end

y=zeros(size(rho));
for i=1:size(betas,2)
    y(:,i) = (ps./rho(:,i)-1)./rho(:,i);
end
FilePtr = fopen('data/outputs','a');
fprintf(FilePtr, '%s\n', datestr(now));
fprintf(FilePtr, '%s\n', num2str(coeffs,'%.1f '));

nm = num2str(coeffs,'%.1f_');
% Scan if file exist
flag = true;
count = 0;
while flag
    if count==0
        fn = strcat('data/',nm(1:end-1),suffix);
    else
        fn = sprintf('data/%s%d%s',nm(1:end),count,suffix);
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
save(fn, 'betas', 'ps', 'coeffs', 'rho', 'y');
