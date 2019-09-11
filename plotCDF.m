%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCDF(distribution1,distribution2,name,xAxis)
    figure
    [y, stats]=cdfplot(distribution1);
    hold on
    [y, stats]=cdfplot(distribution2);
    hold on
    title(name)
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel(xAxis,'FontSize',30,'FontWeight','bold');
    legend('Substrates','Products')
end