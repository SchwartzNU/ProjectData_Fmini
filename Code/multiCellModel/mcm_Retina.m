classdef mcm_Retina < handle
    
    properties
        cells = {}
        extent = 0;
        responsesByCell = [];
        responsesByCellWithNoise = [];
        locationsByCell = [];
        
        density
        cellType
        
        randomAngleMode = 0;
        
        responseDisplayScale = 0.3;
        rfScaleFactor = 80;
        gaussianGridJitter = 10;
        noiseModelSD = 0;
    end
    
    methods
        function obj = mcm_Retina(extent)
            obj.extent = extent;
            
            %             obj.noiseModelSDFit = load('noiseModels.sfit');
        end
        
        function setup(obj)
            fprintf('set up retina\n')
            obj.generateCells('jitteredHexGrid', obj.gaussianGridJitter);
        end
        
        function generateCells(obj, arrangement, args)
            area = (obj.extent * 2 * .001) ^ 2;
            numCells = round(obj.density * area);
                   
            
            obj.locationsByCell = zeros(numCells, 2);
            
            % make locations
            switch arrangement
                case 'random'
                    for ci = 1:numCells
                        loc = 2* obj.extent * rand(1,2) - obj.extent;
                        obj.locationsByCell(ci,:) = loc;
                        
                    end
                case 'jitteredGrid'
                    ci = 0;
                    jitterSigma = args(1);
                    numRowCols = round(1 + sqrt(numCells));
                    for xi = 1:numRowCols
                        for yi = 1:numRowCols
                            loc = ([xi, yi] - numRowCols/2 - 0.5)/numRowCols * obj.extent * 2;
                            loc = loc + randn(1,2) * jitterSigma;
                            
                            ci = ci + 1;
                            obj.locationsByCell(ci,:) = loc;
                        end
                    end
                    
                    
                case 'jitteredHexGrid'
                    %                     function [X,Y] = hexgrid(planeSize,spacing,randPhase)
                    planeSize = obj.extent * 2;
                    
                    spacing = sqrt(planeSize*planeSize*2 / (obj.density * sqrt(3)));
                    
                    randPhase = 1;
                    
                    gridSize = ceil((planeSize./spacing) * 1.5);
                    if rem(gridSize,2) == 1, gridSize = gridSize+1; end
                    Rad3Over2 = sqrt(3) / 2;
                    [X, Y] = meshgrid(1:1:gridSize);
                    n = size(X,1);
                    X = Rad3Over2 * X;
                    Y = Y + repmat([0 0.5],[n,n/2]);
                    %set spacing
                    X = X * spacing;
                    Y = Y * spacing;
                    Ind = (X-spacing*2 <= planeSize & Y-spacing*2 <= planeSize);
                    X = X(Ind);
                    Y = Y(Ind);
                    if randPhase
                        %adding random noise
                        %rand(?seed?,1);
                        L = length(X);
                        X = X+spacing./2.*(rand(L,1)-0.5);
                        Y = Y+spacing./2.*(rand(L,1)-0.5);
                    end
                    %X,Y are now coordinates of bipolar centers
                    obj.locationsByCell = [X, Y] - obj.extent - spacing;
                    % nBipolars = numel(X);
                    
            end
            if size(obj.locationsByCell, 1) > numCells
                obj.locationsByCell = obj.locationsByCell(1:numCells,:);
            end
            numCells = size(obj.locationsByCell, 1);
            
            global SAM_RESEARCH
            s = load([SAM_RESEARCH sprintf('datasets/shapeModelOutputMatrix/shapeModelOutputMatrix_%s.mat', obj.cellType)]);
            responseMatrix = struct();
            responseMatrix.paramGridPoints1 = s.paramGridPoints1;
            responseMatrix.paramGridPoints2 = s.paramGridPoints2;
            responseMatrix.paramGridPoints3 = s.paramGridPoints3;
            s.valueMatrix = s.valueMatrix * obj.rfScaleFactor;%./ prctile(s.valueMatrix(:), 95) * 20; % 90% peak response 20 spikes
            %             s.valueMatrix = s.valueMatrix - min(s.valueMatrix(:));
            responseMatrix.valueMatrix = s.valueMatrix;
            
            % put a cell at each location
            obj.cells = cell(numCells,1);
            for ci = 1:numCells
                newCell = mcm_Cell();
                newCell.type = obj.cellType;
                newCell.location = obj.locationsByCell(ci,:);
                newCell.noiseModelSD = obj.noiseModelSD;
                newCell.responseMatrix = responseMatrix;
                if obj.randomAngleMode
                    newCell.rfRotationAngle = rand(1) * 360;
                end
                %                 newCell.setup();
                obj.cells{ci} = newCell;
            end
            
            fprintf('retina has %g cells in %g**2 area = %g\n', length(obj.cells), obj.extent*2, area);
            
            
        end
        
        function show(obj, showNoise, rfDisplayMode)
            
            locations = obj.locationsByCell;
            if isempty(obj.responsesByCell)
                scatter(locations(:,1), locations(:,2), 20, 'green', 'filled');
            else
                %                 disp(obj.responsesByCell);
                %                 scale = 200/max(obj.responsesByCell)
                if ~showNoise
                    responses = obj.responsesByCell(:,1);
                else
                    responses = obj.responsesByCellWithNoise(:,1);
                end
                
                responseDisplayMin = 0.7;
                %                 obj.responseDisplayScale = .0008;
                
                if strcmp(rfDisplayMode, 'responses')
                    scatter(locations(:,1), locations(:,2), 2, 'black', 'filled');
                    hold on
                    scatter(locations(:,1), locations(:,2), obj.responseDisplayScale*100 * responses + .1, 'green');
                    
                elseif strcmp(rfDisplayMode, 'offset')
                    
                    D = 60;
                    F = 0.4;
                    for ci = 1:size(locations,1)
                        loc = locations(ci,:);
                        
                        %                         w = 1;
                        line([loc(1),loc(1)], loc(2)+[(D*F)/2,-(D*F)/2],'Color',[1,.5,0],'LineWidth',3)
                        eOff = ellipse(D/2, D*(1-F)/2, 0, loc(1), loc(2)-(D*F)/2, 'cyan');
                        eOn = ellipse(D/2, D*(1-F)/2, 0, loc(1), loc(2)+(D*F)/2, 'magenta');
                        
                        if size(obj.responsesByCellWithNoise, 2) == 1
                            w = responses(ci)*obj.responseDisplayScale + 1;
                            eOn.LineWidth = w;
                            eOff.LineWidth = w;
                        else
                            % plot ON and OFF RFs with different activity levels
                            % use generate_activation_figure.m
                            responses = obj.responsesByCellWithNoise;
                            wOff = responses(ci,2)*obj.responseDisplayScale + responseDisplayMin;
                            eOff.LineWidth = wOff;
                            wOn = responses(ci,1)*obj.responseDisplayScale + responseDisplayMin;
                            eOn.LineWidth = wOn;
                        end
                    end
                    
                elseif strcmp(rfDisplayMode, 'nonOffset')
                    
                    for ci = 1:size(locations,1)
                        loc = locations(ci,:);
                        scatter(locations(:,1), locations(:,2), 20, [1,.5,0], 'filled');
                        
                        %                         w = 1;
                        D = 60;
                        eps = 2;
                        eOn = ellipse(D/2+eps, D/2+eps, 0, loc(1), loc(2), 'magenta');
                        eOff = ellipse(D/2-eps, D/2-eps, 0, loc(1), loc(2), 'cyan');
                        
                        if size(obj.responsesByCellWithNoise, 2) == 1 % one response per cell (ON+OFF)
                            w = responses(ci)*obj.responseDisplayScale + responseDisplayMin;
                            eOn.LineWidth = w;
                            eOff.LineWidth = w;
                        else
                            % plot ON and OFF RFs with different activity levels
                            % use generate_activation_figure.m
                            responses = obj.responsesByCellWithNoise;
                            wOn = responses(ci,1)*obj.responseDisplayScale + responseDisplayMin;
                            eOn.LineWidth = wOn;
                            wOff = responses(ci,2)*obj.responseDisplayScale + responseDisplayMin;
                            eOff.LineWidth = wOff;
                        end
                    end
                    
                end
            end
            
        end
        
        function responses = applyStimulus(obj, stimulus)
            stimulus.angle = mod(stimulus.angle, 360);
            responses = zeros(length(obj.cells), 1);
            for ci = 1:length(obj.cells)
                responses(ci,1) = obj.cells{ci}.applyStimulus(stimulus);
                obj.locationsByCell(ci,1:2) = obj.cells{ci}.location;
            end
            
            obj.responsesByCell = responses;
            
        end
        
        function responsesWithNoise = applyNoise(obj)
            responsesWithNoise = zeros(size(obj.cells,1), 1);
            for ci = 1:size(obj.cells,1)
                responsesWithNoise(ci,1) = obj.cells{ci}.applyNoise();
            end
            
            obj.responsesByCellWithNoise = responsesWithNoise;
            
        end
    end
    
end

