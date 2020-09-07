%% multiCellModel

% sam cooler
% in the lab of 
% greg w schwartz

% end of 2019, for paper 1

%% Generate model parameters
% Models will be created a set having N for each combination of parameter values
% Each of N models with those parameters will have a different randomly
% generated set of RGCs
% Model code incorporates parameter values into processing
% See: mcm_ExperimentSet.setup()
%
%

% cellTypeNames = {'FmON_largeOffset','ONOFF'};
% cellTypeNames = {'dim6030_off30_norm','dim6060_off0_norm'};
% cellTypeNames = {'dim6060_off0_inhib','dim6060_off0_norm'};
% cellTypeNames = {'dim6030_off30_inhib','dim6030_off30_norm'};
% cellTypeNames = {'dim60_off0.375', 'dim60_off0'}; % standard
% cellTypeNames = {'dim60_off0','dim60_off0.375','dim60_off0.375_noShrink'};
% cellTypeNames = {'surr_dim60_off0.375','surr_dim60_off0'};
cellTypeNames = {'surr_dim60_off0','surr_noShrink_dim60_off0.375','surr_dim60_off0.375'};
cellTypeNames_igor = {'not','offset','shrink'};

% cellTypeNames = {'surr_dim60_off0'};

% % comparison over RF properties:
% cellTypeNames = {};
% rfOffsets = linspace(0.0, 0.5, 5);
% rfSizes = linspace(30, 90, 5);
% 
% for sizeI = 1:length(rfSizes)
%     for offsetI = 1:length(rfOffsets)
%         cellTypeNames{end+1,1} = sprintf('surr_dim%.2g_off%.3g', rfSizes(sizeI), rfOffsets(offsetI));
%     end
% end

% fname = [SAM_RESEARCH sprintf('datasets/shapeModelOutputMatrix_dim%.2g_off%.3g.mat', exportRfSizes(sizeI), exportRfOffsets(offsetI))];

% run a set of experiments, varying cell/RF type, in this configuration

% run experiments on both cell types
% params vary: stimX, stimY, stimAngle, cellDensity

fprintf('Model parameter config\n')
parameterNames = {'cellType','stimX','stimY','stimAngle','cellDensity','noiseSD','randomRotation'};

cellTypeIndices = 1:length(cellTypeNames);

numAngles = 7;
stimulusAngles = linspace(0,180, numAngles);
% stimulusAngles(end) = [];
% stimulusAngles = [0];

% cellDensities = linspace(150, 450, 7);
cellDensities = 250;


% responseNoises = linspace(0,25,5);
% responseNoises = [0 1.5 3 4.5 6 12 24];
responseNoises = [2];

randomRotations = [0];

numPositions = 100;
randomPositionRange = 150;
positions = randi(2 * randomPositionRange, [numPositions, 2]) - randomPositionRange;
% positions = randn([numPositions, 2]) * randomPositionRange / 2;

parameterValues = [];
for cti = cellTypeIndices
    for pi = 1:numPositions
        for ai = 1:length(stimulusAngles)
            for di = 1:length(cellDensities)
                for ni = 1:length(responseNoises)
                    for rri = 1:length(randomRotations)
                        parameterValues(end+1,:) = [cti, positions(pi,:), stimulusAngles(ai), cellDensities(di), responseNoises(ni), randomRotations(rri)];
                    end                    
                end
            end
        end
    end
end
numParameters = size(parameterValues,2);
numParameterSets = size(parameterValues,1);

fprintf('Made %g parameter sets\n', size(parameterValues,1));


%%

fprintf('Running the paramset of models\n');
fprintf('may our hearts be sturdy and our math be true\n\n');

fprintf('SETUP\n');
exps = mcm_ExperimentSet(parameterValues, parameterNames);
exps.cellTypeNames = cellTypeNames;
tic
exps.setup()
toc
%
fprintf('RUN\n');
tic
exps.run()
toc
exps.experiments{1}.show(21,'responses',1)



%% Post-process analysis and display accuracy results
% looks through results along a single parameter axis
% figure(1);clf;

fprintf('\n\nRunning post-analysis\n');

compareByParam = [];
errorByParamType = [];
errorStdByParamType = [];
loopParameterValues = [];

if length(cellTypeNames) == 25
    loopParameterColumn = 1;
    loopParameterValues = 1:25;
    loopParameterName = 'RF type';
end

if length(cellDensities) > 1
    loopParameterColumn = 5;
    loopParameterValues = cellDensities;
    loopParameterName = 'Density';
end

if length(stimulusAngles) > 1
    loopParameterColumn = 4;
    loopParameterValues = stimulusAngles;
    loopParameterName = 'Angle';
end

if length(responseNoises) > 1
    loopParameterColumn = 6;
    loopParameterValues = responseNoises;
    loopParameterName = 'NoiseSD';
end

% SET ERROR MODE:
errorMode = 'X'; % X, Y, mag

if ~isempty(loopParameterValues)

figure(2);clf;

d1 = ceil(sqrt(length(loopParameterValues)));
d2 = ceil(length(loopParameterValues) / d1);
ha = tight_subplot(d1, d2, .1, .1, .1);


for pi = 1:length(loopParameterValues)
    fprintf('\nLoop param %s = %g\n', loopParameterName, loopParameterValues(pi))
%     
    
    compareByCellType = [];%zeros(numParameterSets, length(cellTypes));
    for cti = cellTypeIndices
        for rri = randomRotations

            exps.resetFilter();
            exps.filterResults(1, cti); % filter cell (RF) type
%             exps.filterResults(4, stimulusAngles(pi));
            exps.filterResults(loopParameterColumn, loopParameterValues(pi));
%             exps.numFilterResults()
            exps.filterResults(7, rri);
%             exps.numFilterResults()
            exps.analyzeResults();
            fprintf('%s\n', exps.resultsFilterString)

    %         figure(1);
    %         exps.showErrorVectors()
    %         hold on

            ai = cti + rri * length(cellTypeNames)
    
            axes(ha(pi))
            compareByCellType(:,ai) = exps.showErrorMagnitude(errorMode, 0); % Var extracted for below analysis across params
            hold on
    %         fprintf('\n');

            errorByParamType(pi,ai) = mean(compareByCellType(:,ai));
            errorStdByParamType(pi,ai) = std(compareByCellType(:,ai)) ./ sqrt(length(compareByCellType(:,ai)));
        end
    end

    axes(ha(pi))
    legend(cellTypeNames,'Interpreter','none')
    title(sprintf('%s: %g', loopParameterName, loopParameterValues(pi)));

    
    if length(cellTypeIndices) > 1
        compareByParam(pi,1) = mean(compareByCellType(:,2)) - mean(compareByCellType(:,1));

        [h,p] = ttest2(compareByCellType(:,1), compareByCellType(:,2),'Vartype','unequal');
        if h == 0
            fprintf('No stat mean difference, p = %.2g, diff = %g\n', p, compareByParam(pi,1));
        else
            fprintf('**Stat mean difference! p = %.2g, diff = %g\n', p, compareByParam(pi,1));
        end
    end
    
    
end
end
%%
figure(40);clf;
for ai = 1:size(errorByParamType,2)
    errorbar(loopParameterValues, errorByParamType(:,ai), errorStdByParamType(:,ai))
    hold on
end
xlabel(loopParameterName)
ylabel('mean error magnitude')
title(sprintf('Comparison of %s error across %s', errorMode, loopParameterName));
legend(horzcat(cellTypeNames,cellTypeNames),'Interpreter','none')
% ylim([0, max(errorByParamType(:))+3])
% ylim([15,25])
%

%%
% os = struct();
randNames = {'fixed','random'};
for ai = 1:size(errorByParamType,2)
    cti = mod(ai-1,3)+1;
    fprintf('%s over %s, "%s" %s (%s)\n',errorMode, loopParameterName, cellTypeNames_igor{cti}, randNames{1+(ai>=4)}, cellTypeNames{cti})
    d = (errorByParamType(:,1) - errorByParamType(:,ai)) ./ errorByParamType(:,1);
    disp(d)
    
%     errorMode = 'Abs';
    nam = sprintf('error%s_angle_%s', errorMode, cellTypeNames_igor{cti});
    if ai >= 4
        nam = [nam '_rr'];
    end
    os.(nam) = d;%errorByParamType(:,ai);
end

%% plot across 2 dims

return
errorByCellType = [];
    
for cti = cellTypeIndices

    exps.resetFilter();
    exps.filterResults(1, cti); % filter cell (RF) type
    exps.analyzeResults();

    errorByCellType(cti) = mean(exps.showErrorMagnitude(errorMode, 0))


end    

errorMap = reshape(errorByCellType,[5,5])'; % transposed from igor version
fig = figure(29);clf;
ax = axes;
imagesc(rfOffsets, rfSizes, errorMap);
ax.YDir = 'normal';
ylabel('RF size (µm)');
xlabel('RF offset (ratio)');
h = colorbar;
ylabel(h, 'Y Error')
colormap gray

mycmap = get(fig,'Colormap');
set(fig,'Colormap',flipud(mycmap))
% save([SAM_RESEARCH 'datasets/multiCellExperimentSet.mat'], 'errorMap')

% s = struct();
% s.errorY_dimOff = errorByParamType_noRotation;
% s.errorY_dimOff_rot = errorByParamType_rotation;
% s.errorY_dimOff_rotDiff = errorByParamType;
% % s.density = cellDensities;
% exportStructToHDF5(s, '~/analysis/igorExport/multiCellModel_2d.h5', 'errorMapRotation');
%% save to Igor
% s = struct();
% s.errorAcrossDensity = errorByParamType;
% s.density = cellDensities;
% exportStructToHDF5(s, '~/analysis/igorExport/multiCellModel.h5', 'multiCellModel');


%% save things
% global SAM_RESEARCH
% save([SAM_RESEARCH 'datasets/multiCellExperimentSet.mat'], 'exps')