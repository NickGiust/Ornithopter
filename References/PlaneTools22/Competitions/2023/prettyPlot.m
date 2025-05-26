function prettyPlot(xDecimal, yDecimal, zDecimal)
%%PRETTYPLOT makes plot latex formatted with good sizes and whatnot
    xtickString = strcat('%.', num2str(xDecimal), 'f');
    xtickformat(xtickString);
    if nargin > 1
        ytickString = strcat('%.', num2str(yDecimal), 'f');
        ytickformat(ytickString);
    end
    if nargin > 2
        ztickString = strcat('%.', num2str(zDecimal), 'f');
        ztickformat(ztickString);
    end
    set(gca, 'box', 'off');
    set(gca,'TickLabelInterpreter','Latex', 'FontSize', 16);
end




