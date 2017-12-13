% draw_rho_h_v2.m

function draw_rho_h_v2(option)
% draw rho-h(rho)
% coeffs = [1, 2.5, 4, 1, 1];
close all;
addpath('..');
if nargin<1
    option = 1;
end

 fnames = {'../data_min/1.0_2.5_4.0_1.0_1.0.dat';
 '../data_min/1.0_2.2_4.0_1.0_1.0.dat';
 '../data_min/1.0_2.0_4.0_1.0_1.00.dat';
}; %  '../data_min/1.0_2.5_4.0_1.0_0.05.dat'
labelA = {'(a)';'(b)';'(c)'};
labelB = {'(d)';'(e)';'(f)'};
rhodatas = {};
for rp=1:length(fnames)
	dn = chopnum(fnames{rp});
    rhodatas{rp}=dlmread(fnames{rp},'',[dn(1) 0 dn(2) 3]);
end
if 1==option
pickedT = [0.1, 0.12, 0.14, 0.16, 0.2, 0.24, 0.3];
lt = length(pickedT);

widthpix = 150;
heightpix = 120;
leftmargin = 42;
bottommargin = 38;
topmargin = 10;
rightmargin = 10;

widthtotal = widthpix*3+leftmargin+rightmargin;
heighttotal = heightpix*2+bottommargin+topmargin;
left = leftmargin/widthtotal;
bottom = bottommargin/heighttotal;
width = widthpix/widthtotal;
height = heightpix/heighttotal;

fig = figure('rend','painters','Position',[100 100 widthtotal heighttotal]);
set(gcf, 'DefaultLineLineWidth', 1);

axbottom = axes('Position',[left bottom 3*width 2*height], 'FontName','Times New Roman','fontsize',14);
plot(nan, nan);
set(gca,'XTick',[0,1]);
set(gca,'XTickLabel',[' ';' ']);
set(gca,'YTick',[]);
% axbottom = axes('Position',[0 0 1 1],'Visible','off');
xlabel('\rho', 'FontSize', 14);

psize=3;%figure; hold on;
for rp=1:length(rhodatas)
    rhodata = rhodatas{rp};
    % draw h(rho)
    axnow=axes('Position', [left+(rp-1)*width, bottom, width, height], 'fontname','Times New Roman','fontsize',14);
    minpoint = zeros(lt,3);
    for ind=1:lt
        tdata = rhodata(rhodata(:,1) == pickedT(ind), 2:4); % data of specific temperature
        [~, inds] = sort(tdata(:,2));
        stdata = tdata(inds, :);     % sorted
        minindex = plots(stdata(:,2), stdata(:,3));
        minpoint(ind,:) = stdata(minindex,:);
        if 1==ind
            hold on;
        end
    end
    if(3==rp)
        set(axnow,'XTick',[1e-8 1e-4 1e0],'FontName','Times New Roman');
    else
        set(axnow,'XTick',[1e-8 1e-4],'FontName','Times New Roman');
    end
    set(axnow,'YTick',[-1e6 -1e3],'FontName','Times New Roman');
    if(1==rp)
        ;
    else
        set(axnow,'YTickLabel',[]);
    end
    % xlabel('\rho'); 
    % mark the minimum
    set(0,'defaultfigurecolor',[0 0 0]);
    for ind=1:lt
        if(rp==1)
            tthres = 0.28;
        elseif(rp==2)
            tthres = 0.18;
        else
            tthres = 0;
        end
        if(pickedT(ind) <=tthres)
            plot(minpoint(ind,2), minpoint(ind,3), '*');
        end
    end
    xlim([1e-8 1e0]);
    ylim([-1e6 -1]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    % box(axnow,'off');
    if(1==rp)
        ylabel('$h(\rho)$', 'Interpreter', 'LaTeX', 'FontSize', 14);
    end
    text(2e-8, -3, labelB{rp}, 'Parent', axnow, 'fontname','Times New Roman','FontSize',14);
    
    % draw p(rho)
    axnow=axes('Position', [left+(rp-1)*width, bottom+height, width, height], 'FontName','Times New Roman','fontsize',14);
    for ind=1:lt
        tdata = rhodata(rhodata(:,1) == pickedT(ind), 2:4); % data of specific temperature
        [~, inds] = sort(tdata(:,2));
        stdata = tdata(inds, :);     % sorted
        plots(stdata(:,2), stdata(:,1));
        if 1==ind
            hold on;
        end
    end
    set(0,'defaultfigurecolor',[0 0 0]);
    if(rp==1)
        tthres = 0.22 ;
        for ind=1:lt
            if(pickedT(ind) <=tthres)
                plot(minpoint(ind,2), minpoint(ind,1), '*');
            end
        end
    end
    xlim([1e-8 1e0]);
    ylim([1e-8 1e-0]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    % box(axnow,'off');
    set(axnow,'YTick',[1e-8 1e-4 1],'FontName','Times New Roman');
    if(1==rp)
        ;
    else
        set(axnow,'YTickLabel',[]);
    end
    set(axnow,'XTick',[1e-8 1e-4]);
    set(axnow,'XTickLabel',[]);
    if(1==rp)
        ylabel('$\beta p$', 'Interpreter', 'LaTeX', 'FontSize', 14);
    end
    text(2e-8, 2e-1, labelA{rp}, 'Parent', axnow, 'fontname','Times New Roman','FontSize',14);
end
%{
leg = legend('T=0.10','T=0.12','T=0.14','T=0.16','T=0.20','T=0.24', 'T=0.30', 'Location','None');
legend boxoff;
set(leg, 'Position', [0.88,0.21,0.1,0.1])
set(leg,'FontName','Times New Roman');
set(leg,'FontSize',12);
%}
% print(fig,'render/rho_h.eps','-depsc');  % has problem 
elseif 2==option
% TODO
end
end

function minindex = plots(x, y)  % log scale (-inf, -10) ~ linear (-10, 10) ~ log(10, inf)
    f1 = y >= -10 & y <=10;
    %y(f1) = y(f1)/10;
    %y(y < -10) = -log10(-y(y < -10));
    %y(y > 10) = log10(y(y > 10));
    [~, ind] = min(y);
    plot(x, y);
    minindex = ind;
end