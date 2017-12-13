%function unifyfigure(size)
    if 1== psize % do not divide
        fig = figure('rend','painters','pos',[100 100 560 420]);
        set(gca,'fontsize',16);
    elseif 2==psize  % 2 subfigure in a row
        fig = figure('rend','painters','pos',[100 100 400 300]);
        set(gca,'fontsize',16);
    elseif 3==psize  % three subfigure in a row
        fig = figure('rend','painters','pos',[100 100 300 200]);
        set(gca,'fontsize',16);
    elseif 12==psize % 2 subfigure but double width
        fig = figure('rend','painters','pos',[100 100 500 300]);
        set(gca,'fontsize',16);
    elseif 21==psize % 2 subfigure but double height
        fig = figure('rend','painters','pos',[100 100 400 450]);
        set(gca,'fontsize',16);
    end
    set(gcf, 'DefaultLineLineWidth', 1.5);
%end