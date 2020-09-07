%% 
overridePositions = 0
rfTypes = {'dim60_off0.375', 'dim60_off0'};
rfPlotTypes = {'offset','nonOffset'};
for cti = 1:2
    cellType = rfTypes{cti};
    ooNames = {'onOnly','offOnly'};


    ooi = 1;
    expe = mcm_Experiment();
    expe.retina.cellType = cellType;

    expe.retina.density = 250;
    expe.stimulus.angle = 0;
    expe.stimulus.location = [0,30];
    expe.setup()

    if cti == 1 && ~overridePositions
        cellLocations = expe.retina.locationsByCell;
    else
        for i = 1:size(cellLocations,1)
            expe.retina.cells{i}.location = cellLocations(i,:);
        end
    end

    expe.run();
    decoded = expe.decodedLocation;
    responses = [];

    % responses has columns [standard, On, Off]
    responses(:,1) = expe.retina.responsesByCellWithNoise;


    for ooi = 1:2
        expe.retina.cellType = [cellType '_' ooNames{ooi}];
        expe.setup();
        for i = 1:size(cellLocations,1)
            expe.retina.cells{i}.location = cellLocations(i,:);
        end
        expe.run();
        responses(:,ooi + 1) = expe.retina.responsesByCellWithNoise;
    end

    expe.decodedLocation = decoded;


    %%
    
    expe.retina.responsesByCellWithNoise = responses(:,1);
    expe.show(10+cti*5+1,'responses', 1)
    set(gca,'xcolor','none')
    set(gca,'ycolor','none')    

    
    % expe.show(21,'nonOffset')
    expe.retina.responsesByCellWithNoise = responses(:,2:3);
    expe.show(10+cti*5+2,rfPlotTypes{cti}, 0)
    % 
    x = 120;
    xlim(x*[-1,1])
    ylim(x*[-1,1]+expe.stimulus.location)
%     set(gca,'xtick',[])
    set(gca,'xcolor','none')
    set(gca,'ycolor','none')
end