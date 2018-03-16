% calccorlen.m
% Calculate correlation length = 1/ln(\lambda1/\lambda2)
function result = calccorlen(coeffs, Ts)
    writefile = false;
    if nargin<1
        coeffs=[1,2.5,4,1,0.6];
    end
    if nargin<2
        %Ts=0.15:0.01:0.5;
        Ts=[0.135, 0.145, 0.155, 0.165];
    end
    % Ts = 0.18;

    option = 'rho'; % 'p' or 'rho'
    if strcmp('rho',option)
        %rhos = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1];
        rhos=[1e-4];
        nx = length(rhos);
    else
        ps = [1e-6 1e-5 1e-4 1e-3];
        nx = length(ps);
    end
    foldername = 'data_corlen';
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
    logfile = strcat(foldername,'/calccorlen.log');
    LogPtr = fopen(logfile, 'a');
    else
       FilePtr = 1;
       LogPtr = 1;
    end
    rq=1;

    fprintf(LogPtr, '-- %s --\n', datestr(clock));
    fprintf(LogPtr, 'parameters: %.2f,%.2f,%.2f,%.2f,%.2f\n', coeffs);
    if writefile
        fprintf(LogPtr, 'T\trho\tp\tlambda1,2,3\n');
    end
    fprintf(FilePtr, 'T\trho\tp\txiL\tlambda2\n');
    for xx=1:nx
        rp=1;
      for T=Ts
        beta=1/T;
        l2status = 0; % 0: real, 1: not clear, 2: complex
          if strcmp('rho',option)
              rho = rhos(xx);
              p = findp(rho, T, coeffs, 1e-5);
          else
              p = ps(xx);
              rho = findrho(p, beta, coeffs);
          end
          %p = rho*T;
          condition = true;
          loopcount = 1;
          while condition
            [qs, ds] = corlen(p, beta, coeffs);
            if l2status==2 % already complex
                break;
            elseif abs(imag(ds(2,2)))>1e-10 % meet a complex
                if loopcount > 5
                    l2status = 2;
                    break;
                else
                    loopcount = loopcount + 1;
                end
            else
                break;  % meet a real
            end
          end
          ks = 1/log(ds(1,1)/abs(ds(2,2)));
          if writefile
            fprintf(LogPtr, '%.3f\t%.3e\t%.3e\t', T, rho, p);
          end
          if abs(imag(ds(2,2)))<1e-10
            fprintf(FilePtr, '%.3f\t%.3e\t%.3e\t%.3e\t%.3e\n', T, rho, p, ks, ds(2,2));
            if abs(imag(ds(3,3)))<1e-10
              fprintf(LogPtr, 'Lambda: %.3e\t%.3e\t%.3e\n', ds(1,1), ds(2,2), ds(3,3));
            else
              fprintf(LogPtr, 'Lambda: %.3e\t%.3e\t%.3g%+.3gi\n', ds(1,1), ds(2,2), real(ds(3,3)), imag(ds(3,3)));
            end
          else
            fprintf(FilePtr, '%.3f\t%.3e\t%.3e\t%.3e\t%.3g%+.3gi\n', T, rho, p, ks, real(ds(2,2)), imag(ds(2,2)));
            fprintf(LogPtr, 'Lambda: %.3e\t%.3g%+.3gi\t%.3g%+.3gi\n', ds(1,1), real(ds(2,2)), imag(ds(2,2)), real(ds(3,3)), imag(ds(3,3)));
          end
          rp=rp+1;
      end
      rq=rq+1;
    end
end
