function MakeDeformPlot( deformData, deformName, TL, Xticks, Yticks, varargin )
if nargin == 6
    deformLims = varargin{1};
else
    deformLims = [min(deformData(:)), max(deformData(:))];
end
if size(deformData,2) == 1
    plot( deformData );
    ylabel(deformName);
    ylim(deformLims)
else
    imagesc( deformData' );
    axPos = get(gca,'Position');
    caxis(deformLims);
    if ~any(isinf(deformLims) | isnan(isinf(deformLims))), colormap(gca, bluewhitered); end
    CB = colorbar('EastOutside'); CB.Label.String = deformName; CB.Label.FontWeight = 'bold';
    set(gca,'Ytick',Yticks, 'Position',axPos);
    ylabel('Plane');
end
set(gca,'TickLength',TL, 'Xtick',Xticks, 'XTickLabel',[], 'TickDir','out', 'FontSize',8)
end