classdef mcm_Cell < handle
    %mcm_Cell A single cell in the retina
    %   
    
    properties
        location % (x,y) in µm
        type % mostly used to choose an RF model

        response % response after applying stimulus
        responseWithNoise % the above with noise added
        
        responseMode = 'matrix' % from NL RF model, using interpolation
        responseMatrix % response by stimulus pos x, y, angle
        
        rfRotationAngle = 0; % random offset angles
        
        noiseModel = 'gaussian'% 'gaussian' usually
        noiseModelSD = 0;%[.0027 1.2866] % linear fit of response standard deviation over spike count
    end
    
    methods
        function obj = mcm_Cell()
        end
        
        function setup(obj)

            
        end
        
        function response = applyStimulus(obj, stimulus)
            
            stimulusY = stimulus.getOffset(obj.location);
            
            switch obj.responseMode
                case 'distance'
                    response = max(200 - sqrt(sum(power(stimulus.location - obj.location,2))), 0);
                
                case 'distanceAbove'
                    response = max(200 - sqrt(sum(power(stimulusY - obj.location(2),2))), 0);
                    
                    if obj.location(2) < stimulusY
                        response = 0;
                    end
                    
                    
                case 'matrix'
                    
                    % generate angle and position to stimulus
                    offset = stimulus.location - obj.location;
                    if obj.rfRotationAngle == 0
                        stimAngle = mod(stimulus.angle, 360);
                    else
                        R = [cosd(obj.rfRotationAngle), -sind(obj.rfRotationAngle); sind(obj.rfRotationAngle), cosd(obj.rfRotationAngle)];
                        offset = (R * offset')';
                        
                        stimAngle = mod(stimulus.angle + obj.rfRotationAngle, 360);
                    end
                    
                    response = interp3(obj.responseMatrix.paramGridPoints1, obj.responseMatrix.paramGridPoints2, ...
                        obj.responseMatrix.paramGridPoints3, obj.responseMatrix.valueMatrix, ...
                        offset(1), offset(2), stimAngle, 'linear', 0);
                    
%                     if response > 0
%                         fprintf('offset %g %g resp %g angle %g\n', offset(1), offset(2), response, angle);
%                     end

            end
%             response = response;
            obj.response = response;
            
        end
        
        function responseWithNoise = applyNoise(obj)
%             disp(obj.noiseModelSD)

%             if obj.response > 0.001 % only add noise if cell is spiking already
                noiseSD = polyval(obj.noiseModelSD, obj.response);

                switch obj.noiseModel
                    case 'gaussian'
                        noise = randn(1,1) * noiseSD;
    %                     obj.response = obj.response + noise;

    %                 case 'poisson'
    %                     noise = 
                end

                % crop responses to zero
                % NB: a nonlinearity

                responseWithNoise = obj.response + noise;
                responseWithNoise = max([responseWithNoise, 0], [], 2);
%             else
%                 responseWithNoise = obj.response;
%             end
            obj.responseWithNoise = responseWithNoise;
            
        end        
        
        function showResponseMatrix(obj)
            clf
            slice(obj.responseMatrix.paramGridPoints1, obj.responseMatrix.paramGridPoints2, obj.responseMatrix.paramGridPoints3, obj.responseMatrix.valueMatrix,0,0,[0,180])
            xlabel('X offset')
            ylabel('Y offset')
            zlabel('angle')
            colormap parula
            colorbar

        end
    end
end

