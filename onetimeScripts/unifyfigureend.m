%function unifyfigureend
    if 1== psize % do not divide
        set(gca,'fontsize',16);
    elseif 2==psize  % 2 subfigure in a row
        set(gca,'fontsize',16);
    elseif 3==psize  % three subfigure in a row
        set(gca,'fontsize',16);
    end
    set(gca,'fontname','Times New Roman');
    box off;
%end