function result=calcsq(rho, grs)
    if(nargin<2)
        rhostr = '0.01';
        Tstr = '0.2';
        rho = str2num(rhostr);
        grs=dlmread(strcat('../simulationdata/vmmc/outputs/1,2.5,4,1,1/T',Tstr,'/gofr/',rhostr,'.dat'),'',1,0);
        %xxx = linspace(1,100,100).';
        %grs=[xxx, (sin(xxx)).^2];
    end
    %plot2(grs)
    L = size(grs,1);
    X = grs(:,2);
    Y = X;
    dr = grs(2,1)-grs(1,1);
    SQ = real(fft(Y))*dr; % 
    SQ = SQ(1:L/2+1);
    SQ(2:end-1) = 2*SQ(2:end-1);
    SQ = 1 + rho*SQ;
    ws = 2*pi/dr/L*(0:(L/2)).';
    if(nargin<2)
        %close all;plot2(grs);xlim([0 12]);figure;
        plot(ws, SQ);
        %plot2(grs(2:end,:));
        xlim([0 20]);
    else
        result = [ws, SQ];
    end
end