% Calculate entropy and heat capacity
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
% Option: 'cp' (isobaric, default) ot 'cv' (isovolume)
function calcentropy(option)
if(nargin<1)
    option = 'cv';  % Default, calculate Cv
end
writefile = true;

%% Initialize A
coeffs=[1,2.2,4,1,1];
Ts=linspace(0.1,1,19);
switch lower(option)
    case 'cp'
        ps=[0.1];
    case 'cv'
        rhos = [0.001 0.01 0.1 0.2 0.5 0.9];
end

foldername = 'data_S';
nm = sprintf('%.1f_%.1f_%.1f_%.1f_%.2f',coeffs(1), coeffs(2), coeffs(3), ...
             coeffs(4), coeffs(5));
         

switch lower(option)
%% Cp specific
case 'cp'
%Ts=linspace(0.1,1,19);
if(writefile)
% Pick up a filename
suffix='_cp.dat';

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
Np = size(ps,2);
S = zeros(length(Ts),Np); 
Cp = zeros(length(Ts),Np); 
rhos = zeros(length(Ts),Np);
fprintf(FilePtr, '%s\t%s\t%s\t%s\t%s\n', 'p', 'T', 'rho', 'S', 'Cp');
for p=ps
  rp=1;
  for T=Ts
    [Cp(rp,rq), S(rp,rq)] = deriv2(@(bT)(bT*log(Pfunc_isobaric_3NN(p, 1/bT, coeffs,200))),T);
    rhos(rp, rq) = 1/deriv(@(bp)(-T*log(Pfunc_isobaric_3NN(bp, 1/T, coeffs, 200))),p);
    Cp(rp,rq) =Cp(rp,rq)*T;
    fprintf(FilePtr, '%.3e\t%.3g\t%.4e\t%.4e\t%.4e\n', p,Ts(rp), rhos(rp,rq),S(rp,rq), Cp(rp,rq));
    rp=rp+1;
  end
  rq=rq+1;
end

case 'cv' 
%% Cv specific
if(writefile)
% Pick up a filename
suffix='_cv.dat';

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
Nr = length(rhos);
ps = zeros(length(Ts),Nr);
S = zeros(length(Ts),Nr);
Cp = zeros(length(Ts),Nr);
Cv = zeros(length(Ts),Nr); 
fprintf(FilePtr, '%s\t%s\t%s\t%s\t%s\t%s\n', 'rho','T','p', 'S', 'Cv', 'Cp');
for rho=rhos
  rp=1;
  for T=Ts
    p=findp(rho, T, coeffs);
    [Cp(rp,rq), deno, nomi, S(rp,rq)] = ...
        deriv2s(@(bT, bp)(bT*log(Pfunc_isobaric_3NN(bp, 1/bT, coeffs, 200))),T, p);
    Cp(rp,rq) = Cp(rp,rq)*T;
    nomi = nomi^2;
    ps(rp,rq) = p;
    Cv(rp, rq) = Cp(rp,rq) - T*nomi/deno;
    fprintf(FilePtr, '%.3e\t%.3g\t%.4e\t%.4e\t%.4e\t%.4e\n', ...
            rho,Ts(rp), ps(rp,rq),S(rp,rq), Cv(rp,rq), Cp(rp,rq));
    rp=rp+1;
  end
  rq=rq+1;
end

end
end

