%% dendrites and somas model


% load FmON network soma positions
load('cellLocationsTable.mat')
sel = contains(cellstr(cellLocationsTable.injectedSomaType), 'on');
cellLocationsTable = cellLocationsTable(sel,:);

% pre-process soma maps
figure(31);clf;
for ci = 1:size(cellLocationsTable,1)

    injectedSoma = cellLocationsTable{ci, 'injectedSomaLocation'};
    somaLocations = cellLocationsTable{ci, 'otherSomaLocations'}{1};
    somaOrders = cellLocationsTable{ci, 'otherSomaOrders'}{1};
    
    selOrder = (somaOrders ~= 2)';
    somaLocations = somaLocations - injectedSoma;
    somaLocations(:,2) = somaLocations(:,2) * -1; % fiji Y flip
    somaLocations = somaLocations(selOrder(1:size(somaLocations,1)),:);
    cellLocationsTable{ci,'otherSomaLocations'} = {somaLocations};
    
    subplot(4,4,ci)
    hold on
    l = 20;
    line([-l, l],[0,0])
    line([0,0],[-l,l])    
    scatter(somaLocations(:,1), somaLocations(:,2));
    axis equal
    names = cellstr(cellLocationsTable.cellName);
    title(names{ci})
end




%% load dendrite polygons
load('dendritePolygonDatabase.mat')
dendritePolygonDatabase_local = dendritePolygonDatabase;

% pre-process dendrite polygons
figure(32);clf;
ha = tight_subplot(10,9);
for ci = 1:size(dendritePolygonDatabase_local,1)

    scalingFactor = dendritePolygonDatabase_local{ci, 'scalingFactor'};
    soma = dendritePolygonDatabase_local{ci, 'soma'} * scalingFactor;
    dendriticPolygon = dendritePolygonDatabase_local{ci, 'polygon'}{1} * scalingFactor;
    dendriticPolygon = dendriticPolygon - soma;
    dendritePolygonDatabase_local{ci,'polygon'} = {dendriticPolygon};
    dendritePolygonDatabase_local{ci,'centroid'} = centroid(dendriticPolygon);
    
%     dendriticPolygonResampled = resamplePolygon(dendriticPolygon,200);
    
%     subplot(5,5,ci)
    axes(ha(ci))
    hold on
    l = 20;
    line([-l, l],[0,0])
    line([0,0],[-l,l])    
    fill(dendriticPolygon(:,1), dendriticPolygon(:,2),'k','FaceAlpha',.2);
    axis equal
    names = dendritePolygonDatabase_local.Row;
    title(names{ci})

end

linkaxes(ha)


%% get FmON centroid offsets
sel = contains(cellstr(dendritePolygonDatabase_local.cellType), 'F-mini ON');
offsets_FmON = dendritePolygonDatabase_local.centroid(sel,:);


%% display plot dendrite polygons
figure(39);clf;
ha = tight_subplot(7,8, -.06, 0,0);
set(gcf,'color','w');
% sel = ones(size(dendritePolygonDatabase_local
% sel = contains(cellstr(dendritePolygonDatabase_local.cellType), 'F-mini ON');
plotOffsets = [0, 40, 60];
cellTypes = {'F-mini ON','F-mini OFF'};

colors = [1, 0, 1;
          0, 1, 1;
          1, 1, 0];

for ai = 1:length(ha)
    set(ha, 'visible', 'off')
end
      
for cti = 1:length(cellTypes)

    sel = contains(cellstr(dendritePolygonDatabase_local.cellType), cellTypes{cti});
    fprintf('num displayed %g\n', sum(sel));
    subset = dendritePolygonDatabase_local(sel,:);
    for ci = 1:size(subset,1)
        axes(ha(ci + plotOffsets(cti)))

        dendriticPolygon = subset{ci,'polygon'}{1};

        hold on

        fill(dendriticPolygon(:,1), dendriticPolygon(:,2),colors(cti,:),'FaceAlpha',.7,'LineStyle','none');
        l = 25;
        line([-l, l],[0,0],'color','k','linewidth',2)
        line([0,0],[-l,l],'color','k','linewidth',2)   
        
        axis equal
        names = dendritePolygonDatabase_local.Row;
%         title(names{ci})
        box off
        set(gca, 'visible', 'off')
%         set(findall(gca, 'type', 'text'), 'visible', 'on')

    end
end

linkaxes(ha)
% xlim(

%% big loop:

% % filter to just FmOFF in table
sel = contains(cellstr(dendritePolygonDatabase_local.cellType), 'F-mini OFF');
dendritePolygonDatabase_local = dendritePolygonDatabase_local(sel,:);

spatialDimN = 200;
spatialDimRange = 250;
X = linspace(-spatialDimRange,spatialDimRange,spatialDimN);
Y = X;
Xl = [];
Yl = [];
% [X,Y] = meshgrid(,linspace(-250,250,N));
for xi = 1:spatialDimN
    for yi = 1:spatialDimN

        Xl(end+1,1) = X(xi);
        Yl(end+1,1) = Y(yi);
        
    end
end



figure(2);clf;
ha = tight_subplot(3,3);

figure(3);clf;
hb = tight_subplot(3,3);

drawPlots = 0;

numLoops = 5000;
allMaps = zeros(spatialDimN*spatialDimN, numLoops);
allCenters = zeros(numLoops, 2);

for li = 1:numLoops
% choose random position set

% if li

    r1 = randi(size(cellLocationsTable, 1));
    r2 = randi(size(offsets_FmON,1));
    somaPositions = cellLocationsTable{r1, 'otherSomaLocations'}{1} + -1*offsets_FmON(r2,:); % move somas by FmON ON dendrites offset
    numSomas = size(somaPositions, 1);
    
% choose random polygon set
    polygons = {};
    for si = 1:numSomas
        r2 = randi(size(dendritePolygonDatabase_local,1));
        polygons(si,1) = dendritePolygonDatabase_local{r2,'polygon'};
    end
    
    if drawPlots
        axes(ha(li))
        colors = distinguishable_colors(numSomas);
        scatter(somaPositions(:,1), somaPositions(:,2), 70, colors)
        hold on
        %draw soma
        l = 50;
        line([-l, l],[0,0])
        line([0,0],[-l,l])
        for pi = 1:numSomas

            %draw polygon
            dendriticPolygon = polygons{pi,1} + somaPositions(pi,:);
            fill(dendriticPolygon(:,1), dendriticPolygon(:,2), colors(pi,:),'FaceAlpha',.3);

        end
        axis equal
    end
    
    % convert polygons into a composed flat map
    outs = [];
    for pi = 1:numSomas
        dendriticPolygon = polygons{pi,1} + somaPositions(pi,:);

        outs(:,pi) = inpolygon(Xl, Yl, dendriticPolygon(:,1),dendriticPolygon(:,2));
        
    end
    outMap = sum(outs, 2);
    outMap = outMap ./ max(outMap(:));
    allMaps(:,li) = outMap;
    
    outMap = reshape(outMap,[spatialDimN,spatialDimN]);
    
    cent = (fliplr(centerOfMass(outMap)) / spatialDimN) - 0.5;
    cent = cent * spatialDimRange * 2;
    allCenters(li,:) = cent;
    
    if drawPlots
        axes(hb(li))
        surf(X,Y,outMap);
        view([0,0,1])
        axis equal
        shading interp
        l = 50;
        line([-l, l],[0,0],'Color',[1,1,1])
        line([0,0],[-l,l],'Color',[1,1,1])    
    end
end
%%
rfMapTogether = mean(allMaps,2);
rfMapTogether = reshape(rfMapTogether,[spatialDimN,spatialDimN]);
figure(4);clf;
surf(X,Y,-1*ones(spatialDimN,spatialDimN),rfMapTogether);
view([0,0,1])
axis equal
shading interp
l = 50;
line([-l, l],[0,0],'Color',[1,1,1])
line([0,0],[-l,l],'Color',[1,1,1])   

% centroid
% A = rfMapTogether;
% s = regionprops(true(size(A)), A, 'WeightedCentroid');
% cent = cat(1, s.WeightedCentroid)
cent = (fliplr(centerOfMass(rfMapTogether)) / spatialDimN) - 0.5;
cent = cent * spatialDimRange * 2;
fprintf('center of mass: x %g, y %g, length %g, angl %g\n', cent(1), cent(2), sqrt(sum(cent.^2)), atand(cent(2)/cent(1)))

line([-l, l] + cent(1),[0,0] + cent(2),'Color',[1,0,1])
line([0,0] + cent(1),[-l,l] + cent(2),'Color',[1,0,1])



%% igor export map
% outStruct = struct();
% outStruct.X = X;
% outStruct.model_rfMapTogether = rfMapTogether;
% 
% fname = 'igorExport/model_rf_map.h5';
% dataLabel = sprintf('model_rfMapTogether');
% 
% try
%     exportStructToHDF5(outStruct, fname, dataLabel)  
%     disp('file written')
% catch
%     warning('no write')
% end