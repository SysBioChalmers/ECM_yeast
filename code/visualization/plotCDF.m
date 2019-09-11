%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCDF(distribution1,distribution2,name,xAxis,filename)
    figure
    [y, stats]=cdfplot(distribution1);
    set(y, 'LineWidth', 4)
    hold on
    [y, stats]=cdfplot(distribution2);
    set(y, 'LineWidth', 4)
    hold on
    title(name,'FontSize',26)
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel(xAxis,'FontSize',30,'FontWeight','bold');
    legend('Substrates','Products')
    set(gca, 'XScale', 'log')
    savefig(filename)
end