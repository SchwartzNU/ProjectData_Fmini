function plotSpatialData(mapX, mapY, d)
    
    dif = abs(mapX(1,1) - mapX(2,1))/2;
    surface(mapX - dif, mapY - dif, zeros(size(mapY)), d, 'LineStyle','none');
    grid off
    set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')