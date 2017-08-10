% draw p vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
function calcrho()
close all;
writefile = false;
coeffs=[1,2.5,4,1,1];
Ts=[0.3 0.5 1];
betas=1 ./ Ts;
%ps=[linspace(0.01,0.1,10), linspace(0.1,1,10), 1.5, linspace(2, 10, 9)];
ps=[1.25 2.25];

if(writefile)
suffix='.dat';

% Pick up a filename
suffix='.dat';
foldername = 'data_rho';
nm = sprintf('%.1f_%.1f_%.1f_%.1f_%.2f',coeffs(1), coeffs(2), coeffs(3), ...
             coeffs(4), coeffs(5));
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
% use pressure, not beta*pressure here
rq=1;
Nbeta = size(betas,2);
rho = zeros(length(ps),Nbeta); % NNN
for beta=betas
  rp=1;
  for p=ps/beta
    rho(rp,rq) = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p);
    rp=rp+1;
  end
  rq=rq+1;
end

ys = zeros(length(ps),Nbeta);
for i=1:Nbeta
    ys(:,i) = (ps(:)./rho(:,i)-1)./rho(:,i);
end


fprintf(FilePtr, '%s\t%s\t%s\t%s\n', 'T', 'beta*p', 'rho', 'y');
for i=1:Nbeta
    raw_record(FilePtr, Ts(i), ps, rho(:,i), ys(:,i));
end

%save(fn, 'betas', 'pr', 'coeffs', 'rho', 'y');
end

%% Functions
% record the raw data
function raw_record(FilePtr, T, ps, rho, y)
rp = 1;
for pbeta=ps
    fprintf(FilePtr, '%.3f\t%.3e\t%.3e\t%.3e\n', T, ps(rp), rho(rp), y(rp));
    rp = rp + 1;
end
end
