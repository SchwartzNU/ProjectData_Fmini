classdef mcm_ExperimentSet < handle
    properties
        parameterValues
        numParamSets
        
        parameterNames
        cellTypeNames
        
        experiments
        decodedLocationByParam
        analysisResults
        resultsFilter
        resultsFilterString
        
    end
    
    methods
        function obj = mcm_ExperimentSet(parameterValues,parameterNames)
            obj.parameterNames = parameterNames;
            obj.numParamSets = size(parameterValues,1);
            obj.parameterValues = parameterValues;
            
            % init one experiment for each param set
            obj.experiments = cell(obj.numParamSets, 1);
            for paramIndex = 1:obj.numParamSets
                
                experiment = mcm_Experiment();          
                obj.experiments{paramIndex,1} = experiment;
            end
        end
        
        
        function setup(obj)
            % make one experiment for each parameter set
            fprintf('Experiment setup\n');
            obj.resetFilter();
            
            experimentsOut = cell(obj.numParamSets,1);
            
            for paramIndex = 1:obj.numParamSets
                fprintf('Setup experiment from param set %g / %g\n',paramIndex, obj.numParamSets);
                fprintf('%g,',obj.parameterValues(paramIndex,:));
                fprintf('\n');
                experiment = obj.experiments{paramIndex,1};
                experiment.retina.cellType = obj.cellTypeNames{obj.parameterValues(paramIndex,1)};
                experiment.retina.density = obj.parameterValues(paramIndex,5);
                experiment.retina.noiseModelSD = obj.parameterValues(paramIndex,6);
                experiment.retina.randomAngleMode = obj.parameterValues(paramIndex,7);
                stimulusAngle = mod(obj.parameterValues(paramIndex,4), 360);
                experiment.stimulus.angle = stimulusAngle;
                experiment.stimulus.location = obj.parameterValues(paramIndex, 2:3);
                experiment.setup();
                
                experimentsOut{paramIndex,1} = experiment;
            end
            obj.experiments = experimentsOut;
        end
        
        function run(obj)
            fprintf('Experiment run\n');
            obj.decodedLocationByParam = zeros(obj.numParamSets,2);
%             experimentsOut = obj.experiments;
            
            for paramIndex = 1:obj.numParamSets
                fprintf('Run experiment from param set %g / %g\n',paramIndex, obj.numParamSets);
                
                experiment = obj.experiments{paramIndex,1};
                experiment.run()
                
                obj.decodedLocationByParam(paramIndex,:) = experiment.decodedLocation;
%                 experimentsOut{paramIndex,1} = experiment;
            end
            
            fprintf('run result in %g outs\n', length(obj.decodedLocationByParam))
%             obj.experiments = experimentsOut;
        end
        
        function filterResults(obj, column, value)
            obj.resultsFilter = obj.resultsFilter & obj.parameterValues(:,column) == value;
            if column == 1
                obj.resultsFilterString = horzcat(obj.resultsFilterString, sprintf('%s: %s, ', obj.parameterNames{column}, obj.cellTypeNames{value}));
            else
                obj.resultsFilterString = horzcat(obj.resultsFilterString, sprintf('%s: %g, ', obj.parameterNames{column}, value));
            end
        end
        
        function c = numFilterResults(obj)
            c = sum(obj.resultsFilter);
        end
        
        function resetFilter(obj)
            obj.resultsFilter = true(obj.numParamSets,1);
            obj.resultsFilterString = '';
        end
        
        
        function analyzeResults(obj)
            obj.analysisResults = struct;
            errorVectors = obj.decodedLocationByParam - obj.parameterValues(:,2:3);
            errorVectorMean = mean(errorVectors,1);
%             errorVectors = errorVectors - errorVectorMean; % remove systematic offset
            
            errorMagnitude = horzcat(abs(errorVectors), sqrt(sum(power(errorVectors, 2),2)));
            
            obj.analysisResults.errorMagnitude = errorMagnitude;
            obj.analysisResults.errorVectors = errorVectors;
%             obj.analysisResults.errorMagnitudeMean = mean(errorMagnitude);
        end
        
        
        function showErrorVectors(obj)
%             rf = obj.resultsFilter;
%             x = obj.parameterValues(obj.resultsFilter,2)
%             u = obj.analysisResults.errorVectors(obj.resultsFilter,1)
            quiver(obj.parameterValues(obj.resultsFilter,2),obj.parameterValues(obj.resultsFilter,3), ...
                obj.analysisResults.errorVectors(obj.resultsFilter,1), obj.analysisResults.errorVectors(obj.resultsFilter,2))
            xlabel('offset (µm)');
        end
        
        function data = showErrorMagnitude(obj,direction,drawPlot)
            switch direction
                case 'X'
                    col = 1;
                case 'Y'
                    col = 2;
                case 'mag'
                    col = 3;
            end
            
            data = obj.analysisResults.errorMagnitude(obj.resultsFilter,col);
            if drawPlot
                histogram(data,10, 'Normalization', 'pdf');
    %             plot(cumsum(data))
                xlabel(sprintf('error %s (µm)', direction));
            end

            fprintf('mean %g, std %g, median %g\n', mean(data), std(data), median(data))
        end
    end
end