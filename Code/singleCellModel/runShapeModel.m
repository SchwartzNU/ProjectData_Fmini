%% SHAPE MODEL
% Sam Cooler 2016

% cellName = '060716Ac2';
% acName = '348';

cellName = '060216Ac2'; % original good WFDS
acName = '1032';


% cellName = '033116Ac2'; % nice RF with edges and bars, but missing bars spikes and inhibitory temporal align
% acName = '263';


% on alpha
% cellName = '051216Ac9';
% acName = '933';



useRealRf = 0;
useRealFilters = 0;
useStimulusFilter = 0;
useSubunits = 0;
useOffFilters = 1;
rectifyExOn = 1;
rectifyExOff = 1;
useInhibition = 1;
makeLightInputMonotonic = 0;

plotSpatialRF = 0;
plotStimulus = 0;
plotStimulusFramesOverParameterSets = 0;
plotFilters = 0;
plotSubunitCurrents = 0;
plotOutputCurrents = 1;
plotCellResponses = 0;
plotOutputNonlinearity = 0;
plotResultsByOptions = 1;
paramResponsePlotMode2d = 'lines';
plotSummaryVarOverParamsets = 1;

runInParallelPool = 0;

saveOutputSignalsToHDF5 = 0;
outputHDF5Name = sprintf('shapeModelOutput_%s.h5', cellName);
outputStruct = struct();

% imgDisplay = @(X,Y,d) imagesc(X,Y,flipud(d'));
% imgDisplay2 = @(mapX, mapY, d) (surface(mapX, mapY, zeros(size(mapY)), d), grid off);
normg = @(a) ((a+eps) / max(abs(a(:))+eps));
plotGrid = @(row, col, numcols) ((row - 1) * numcols + col);
calcDsi = @(angles, values) abs(sum(exp(sqrt(-1) * angles) .* values) / sum(values));


if ~exist('loadedCellName', 'var') || ~strcmp(loadedCellName, cellName) || ~strcmp(loadedAcName, acName)
    disp('Loading cell filters and responses')
%     extractMovingBarResponses
    loadedCellName = cellName;
    loadedAcName = acName;

    if useRealFilters
        extractFilters
    end
end


stim_mode = 'driftingGrating';
% stim_mode = 'movingBar';
% stim_mode = 'driftingGrating'; 
% stim_mode = 'flashedEdge';
% stim_mode = 'edgeGradient'; % mcm

%% Big parameter loop around overall parameters



paramColumnNames = {'rfOffset', 'barSpeeds'};
% col_barSpeed = 2;
col_rfSize = 2;
col_rfOffset = 1;

% paramColumnsNames = {'edge polarity','angle'};
% col_edgeFlip = 1;
% col_edgeAngle = 2;
%'spatial offset', 'ex delay', };

% col_filterDelay = 1;

% paramValues  = {0;
%                 10;
%                 50};

% paramValues  = {195;
%                 585;
%                 975;
%                 1365};


% paramValues = [.1, 20];

% paramValues = [0;
%                2;
% %                8;
%                16;
% %                24;
%                35;
%                52;
%                75];
            

% paramValues = [.04; 
%                .05;
%                .06;
%                .07;
%                .08];

paramValues = linspace(0.4, 1, 8)';
paramColumnsNames = {'rf aspect'};
col_rfAspect = 1;
           
% p1 = [.04,.07,.1];
% p2 = [2,12,32];
% 
% paramValues = [0, 0]; 

% paramValues = [0, 0
%                1, 0
%                0, 90
%                1, 90];

% 
% p1 = [250, 500, 1000, 2000];
% p1 = [250];
% % p2 = [6 12 18 24 28 34];
% % p2 = round(linspace(6,36,5));
% p2 = [0, 10, 20, 30, 40];
% p2 = [30];
% paramValues = [];
% for i1 = 1:length(p1)
%     for i2 = 1:length(p2)
%         paramValues(end+1,1:2) = [p1(i1),p2(i2)];
%     end
% end

% paramValues = [0, 10, 20, 30, 40]'*3;

% parameters for 3D map: stimulus X, Y, angle
% paramColumnsNames = {'locX','locY','edgeAngle'};
% inner = [linspace(-100, -10, 6), 0, linspace(10, 100, 6)];
% p1 = [linspace(-220, -120, 3), inner, linspace(120, 220, 3)];
% p2 = p1;
% p3 = linspace(0, 360, 16+1);
% 
% % %test:
% % p1 = linspace(-220,220,8);
% % p2 = p1;
% % p3 = linspace(0,360,8+1);
% 
% p1 = 0;
% p2 = 0;
% p3 = 0;
% 
% % p3(end) = [];
% paramValues = [];
% paramIndices = [];
% for i1 = 1:length(p1)
%     for i2 = 1:length(p2)
%         for i3 = 1:length(p3)
%             paramValues(end+1,1:3) = [p1(i1),p2(i2),p3(i3)];
%             paramIndices(end+1,1:3) = [i1, i2, i3];
%         end
%     end
% end


% multiCellModel setup:
% exportRfOffsets = [0,15,30,45,60];
% rfOffsets = linspace(0.0, 0.5, 5);
% rfOffsets = [0, .375];
% rfOffsets = exportRfOffsets;
rfOffsets = 0;
% rfSizes = linspace(30,90,5);
rfSizes = 60;

[numParamSets,numParams] = size(paramValues);

summaryVarByParamSet = [];

tic
fprintf('Running full simulation with %g paramsets\n', numParamSets);


for sizeI = 1:length(rfSizes)
    for offsetI = 1:length(rfOffsets)
        fprintf('dim %.2g, offset %.3g\n', rfSizes(sizeI), rfOffsets(offsetI))
        
        for paramSetIndex = 1:numParamSets
            if numParamSets < 100 || (mod(paramSetIndex, 500) == 0 || paramSetIndex == 0)
                fprintf('Param Set %d: %s\n', paramSetIndex, mat2str(paramValues(paramSetIndex,:)));
            end

            %% Setup
            
            shapeModelParameterSetup

            if ~useRealFilters
               shapeModelParameterizedFilters 
            end

            shapeModelStimSetup

            if useStimulusFilter
                shapeModelStimulusFilter
            end

            if plotStimulusFramesOverParameterSets
                shapeModelPlotStimFrames
            end

            %% Run simulation
            shapeModelSim

            %% Rescale currents and combine, then extract parameters
            shapeModelAnalyzeOutput

        end

        disp('Model run complete')

        shapeModelAnalyzeFinal
        
        toc
    end
end
toc



%% process output nonlinearity

% speeds = cell2mat(paramValues)';
% anglesRads = deg2rad(stim_barDirections);
% 
% inputs = [0,20.5,30];
% outputs = [0,7,45];
% 
% 
% figure(97)
% plot(inputs, outputs)
% 
% for ps = 1:numParamSets
%     lin = valuesByParamSet(ps,:);
%     
%     nonlin = interp1(inputs, outputs, lin);
%     nonlinValuesByParamSet(ps,:) = nonlin;
%     nonlinsummaryVarByParamSet(ps,1) = calcDsi(anglesRads, nonlin);
% end
% 
% 
% figure(110);
% clf;
% ha = tight_subplot(1,2,.1,.1);
% for ps = 1:numParamSets
%     axes(ha(1))
%     polar(anglesRads,valuesByParamSet(ps,:))
%     hold on
%     % ylim([0,max(valuesByParamSet(:))+2])
%     
%     axes(ha(2))
%     vals = nonlinValuesByParamSet(ps,:);
%     polar(anglesRads, vals./max(vals))
%     hold on
%     % ylim([0,max(nonlinValuesByParamSet(:))+2])
% end
% axes(ha(1))
% legend('195','585','975','1365','Location','best');
% % ylim([0,max(valuesByParamSet(:))+2])
% title('before nonlin')    
% 
% 
% axes(ha(2))
% title('after nonlin')
% legend('195','585','975','1365','Location','best');  
% 
% 
% sanes = [0.2
% 0.35
% 0.12
% 0.05];
% 
% figure(111);clf;
% plot(cell2mat(paramValues)', sanes);
% hold on
% plot(cell2mat(paramValues)', summaryVarByParamSet);
% plot(cell2mat(paramValues)', nonlinsummaryVarByParamSet);
% 
% legend('sanes','model','modelNonlin')
% hold off

%% Save HDF5 file
if saveOutputSignalsToHDF5
    delete(outputHDF5Name);
    exportStructToHDF5(outputStruct, outputHDF5Name, '/');
end