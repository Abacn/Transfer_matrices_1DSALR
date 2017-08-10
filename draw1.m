% draw p vs. rho lines in different temperatures
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]

close all;
writefile = false;
coeffs=[1,1.5,1.5,0,1];
ps=[ 0.001];
%Ts=linspace(0.1,1,10);
Ts=linspace(0.01,1,50);
if(writefile)
suffix='.dat';

% Pick up a filename
suffix='.dat';
foldername = 'data_S';
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
Np = size(ps,2);
S = zeros(length(Ts),Np); 
Cp = zeros(length(Ts),Np); 
Pfunc = zeros(length(Ts),Np);
for p=ps
  rp=1;
  for T=Ts
    %[Cp(rp,rq), S(rp,rq)] = deriv2(@(bT)(bT*log(Pfunc_isobaric_3NN(p, 1/bT, coeffs, 200))),T);
    %Cp(rp,rq) =Cp(rp,rq)*T;
    [~, Pfunc(rp,rq)] = deriv2(@(bT)(bT*log(Pfunc_isobaric_HS(p, 1/bT, coeffs, 200))),T);
    Cp(rp,rq) = log(T^1.5/p)+1.5;
    rp=rp+1;
  end
  rq=rq+1;
end


%fprintf(FilePtr, '%s\t%s\t%s\t%s\n', 'p', 'T', 'S', 'Cp');
for i=1:Np
    %raw_record(FilePtr, ps(i), Ts, S(:,i), Cp(:,i));
    hold on;
    plot(Ts, Pfunc(:,i));
    plot(Ts, Cp(:,i), 'x');
end
xlabel('T');ylabel('S');

%% Functions
% record the raw data
function raw_record(FilePtr, p, Ts, S, Cp)
rp = 1;
for pbeta=1:length(S)
    fprintf(FilePtr, '%.3e\t%.3g\t%.4e\t%.4e\n', p,Ts(rp), S(rp), Cp(rp));
    rp = rp + 1;
end
end
