


% a word about meanLevel:
% this is a contrast-based encoding of the light stimulus. So the meanLevel
% will be zero.


switch stim_mode
    case 'movingBar'
        numAngles = 12;
        stim_directions = linspace(0,360,numAngles+1);
        stim_directions(end) = [];
%         stim_directions = [0, 270, 90];
        
        stim_numOptions = length(stim_directions);
        
        stim_barSpeed = paramValues(paramSetIndex,col_barSpeed);
        stim_barLength = stim_barSpeed * 0.8*1.1;
        
        stim_barWidth = 200;
        stim_intensity = 0.5;
        stim_meanLevel = 0;
        
    case 'flashedSpot'
        % SMS
        numSizes = 12;
        stim_spotDiams = logspace(log10(30), log10(1200), numSizes);
        stim_numOptions = length(stim_spotDiams);
        stim_spotDuration = 0.7;
        stim_spotStart = 0.1;
        stim_intensity = 0.5;
        stim_spotPosition = [0,0];
        stim_meanLevel = 0; % doesn't work well
        
        
    case 'driftingGrating'
        numAngles = 12;
        stim_directions = linspace(0,360,numAngles+1);
        stim_directions(end) = [];
%         stim_directions = [0, 30];
        
        stim_numOptions = length(stim_directions);
        
        stim_texSpeed = 500;
        stim_moveTime = .5;
        stim_meanLevel = 0;
        stim_contrast = 1;
        stim_textureScale = 250;
        stim_movementDelay = 0;
        stim_timeOffset = 0.3;
        
    case 'flashedEdge'
        % stim_edgeSpacing = 20;
        % stim_positions = [-120, -90, -60, -30, 0, 30, 60, 90, 120];
        stim_positions = linspace(-130, 130, 7);
        % stim_positions = [140];
        
        stim_numOptions = length(stim_positions);
        stim_edgeAngle = paramValues(paramSetIndex, col_edgeAngle);
        stim_edgeFlip = false;
        stim_contrastSide1 = 1;
        stim_contrastSide2 = 0;
        % stim_meanLevel = 1; % assume the mean is constant and equal to 1
        stim_startTime = 0.01;
        stim_stimTime = 0.3;
        stim_fullField = 1;
        
    case 'edgeGradient'
        
        % numOffsets = 20;
        stim_numOptions = 1;
%         stim_offset = [0, 0];
        
        stim_contrast = 1;
        stim_meanLevel = 0;
        stim_startTime = 0.02;
        stim_stimTime = 0.3;
%         stim_angle = 0;
        stim_width = 150;
        stim_singleSideLength = 150;
        
        stim_offset = paramValues(paramSetIndex, 1:2);
        stim_angle = paramValues(paramSetIndex, 3);
        
end


% stim_mode = 'flashedSpot';
% stim_numOptions = 1;
% stim_spotDiam = 200;


% stim_mode = 'driftingTexture';
% numAngles = 9;
% stim_directions = linspace(0,360,numAngles+1);
% stim_directions(end) = [];
% stim_numOptions = length(stim_directions);
%
% stim_texSpeed = 500;
% stim_moveTime = sim_endTime + 1.0;
% stim_meanLevel = 0.5;
% stim_uniformDistribution = 1;
% stim_resScaleFactor = 2;
% stim_randomSeed = 1;
% stim_textureScale = 30;
% stim_movementDelay = 0.5;

% stim_mode = 'driftingGrating';
% numAngles = 1;
% stim_directions = linspace(0,360,numAngles+1);
% stim_directions(end) = [];
% stim_numOptions = length(stim_directions);
% %
% stim_texSpeed = 500;
% stim_moveTime = sim_endTime + 1.0;
% stim_meanLevel = 0.5;
% stim_textureScale = 30;
% stim_movementDelay = 0.5;

stim_lightMatrix_byOption = {};

for optionIndex = 1:stim_numOptions
    
    %% Setup stimulus
    center = [0,0];
    
    stim_lightMatrix = zeros(sim_dims); %#ok<*PFTUS>
    % stim_lightMatrix_byOptions = cell(num
    
    switch stim_mode
        case 'flashedSpot'
            % flashed spot
            stim_spotDiam = stim_spotDiams(optionIndex);
            
            
            pos = stim_spotPosition + center;
            for ti = 1:sim_dims(1)
                t = T(ti);
                if t >= stim_spotStart && t < stim_spotStart + stim_spotDuration
                    
                    for xi = 1:sim_dims(2)
                        x = X(xi);
                        for yi = 1:sim_dims(3)
                            y = Y(yi);
                            
                            val = stim_intensity;
                            % circle shape
                            rad = sqrt((x - pos(1))^2 + (y - pos(2))^2);
                            if rad < stim_spotDiam / 2
                                stim_lightMatrix(ti, xi, yi) = val;
                            end
                        end
                    end
                end
            end
            
        case 'flashedEdge'
            if ~(stim_edgeFlip)
                c1 = stim_contrastSide1;
                c2 = stim_contrastSide2;
            else
                c1 = stim_contrastSide2;
                c2 = stim_contrastSide1;
            end
            if stim_fullField
                for ti = 1:sim_dims(1)
                    t = T(ti);
                    if t > stim_startTime && t < stim_startTime + stim_stimTime
                        
                        if stim_edgeAngle == 0
                            for yi = 1:sim_dims(3)
                                y = Y(yi);
                                if y > stim_positions(optionIndex)
                                    stim_lightMatrix(ti, :, yi) = c1;
                                else
                                    stim_lightMatrix(ti, :, yi) = c2;
                                end
                            end
                        elseif stim_edgeAngle == 90
                            for xi = 1:sim_dims(2)
                                x = X(xi);
                                if x > stim_positions(optionIndex)
                                    stim_lightMatrix(ti, xi, :) = c1;
                                else
                                    stim_lightMatrix(ti, xi, :) = c2;
                                end
                            end
                        end
                    end
                end
            else
                error('no non full field code yet')
            end
            
        case 'movingBar'
            
            stim_barDirection = stim_directions(optionIndex); % degrees
            
            
            % make four corner points
            l = stim_barLength / 2;
            w = stim_barWidth / 2;
            corners = [-l,w;l,w;l,-w;-l,-w];
            
            % rotate corners
            theta = stim_barDirection; % not sure if this should be positive or negative... test to confirm
            R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            for p = 1:size(corners, 1)
                corners(p,:) = (R * corners(p,:)')';
            end
            
            
            % translate corners
            movementVector = stim_barSpeed * [cosd(stim_barDirection), sind(stim_barDirection)];
            for ti = 1:sim_dims(1)
                
                barCenter = center + movementVector * (T(ti) - sim_endTime / 2 + .2);
                cornersTranslated = bsxfun(@plus, corners, barCenter);
                
                % what a nice thing to find just the right MATLAB function #blessed
                stim_lightMatrix(ti,:,:) = inpolygon(mapX, mapY, cornersTranslated(:,1), cornersTranslated(:,2));
                
            end
            
            %     elseif strcmp(stim_mode, 'driftingTexture')
            %         direction = stim_directions(optionIndex);
            %         canvasSize = sim_dims(2:3);
            %         sigmaPix = 0.5 * (stim_textureScale / sim_spaceResolution);
            %         distPix = (stim_texSpeed/sim_spaceResolution) * stim_moveTime; % um / sec
            %         moveDistance = distPix;
            %         res = [max(canvasSize) * 1.42,...
            %             max(canvasSize) * 1.42 + distPix]; % pixels
            %         res = round(res);
            %
            % %         M = cos();
            
            
        case {'driftingTexture','driftingGrating'}
            direction = stim_directions(optionIndex);
            canvasSize = sim_dims(2:3);
            sigmaPix = 0.5 * (stim_textureScale / sim_spaceResolution);
            distPix = (stim_texSpeed/sim_spaceResolution) * stim_moveTime; % um / sec
            moveDistance = distPix;
            res = [max(canvasSize) * 1.42, max(canvasSize) * 1 + distPix]; % pixels
            res = round(res);
            
            switch stim_mode
                case 'driftingTexture'
                    stream = RandStream('mt19937ar','Seed',stim_randomSeed);
                    M = randn(stream, res);
                    defaultSize = 2*ceil(2*sigmaPix)+1;
                    M = imgaussfilt(M, sigmaPix, 'FilterDomain','frequency','FilterSize',defaultSize*2+1);
                    
                    
                    if stim_uniformDistribution
                        bins = [-Inf prctile(M(:),1:1:100)];
                        M_orig = M;
                        for i=1:length(bins)-1
                            M(M_orig>bins(i) & M_orig<=bins(i+1)) = i*(1/(length(bins)-1));
                        end
                        M = M - min(M(:)); %set mins to 0
                        M = M./max(M(:)); %set max to 1;
                        M = M - mean(M(:)) + 0.5; %set mean to 0.5;
                    else % normal distribution
                        M = zscore(M(:)) * 0.3 + 0.5;
                        M = reshape(M, res);
                        M(M < 0) = 0;
                        M(M > 1) = 1;
                    end
                case 'driftingGrating'
                    
                    xg = 0:res(1); %in µm
                    yg = ((cos(xg * 2 * 3.141 / sigmaPix)) * stim_contrast + 1) - 1;
                    M = repmat(yg, [res(2),1]);
                    
                    % square wave
                    %                 M(M>stim_meanLevel) = stim_meanLevel * (1+stim_contrast);
                    %                 M(M<=stim_meanLevel) = stim_meanLevel * (1-stim_contrast);
                    
            end
            
            for ti = 1:sim_dims(1)
                if T(ti) < stim_movementDelay
                    tMove = 0;
                else
                    tMove = T(ti) - stim_movementDelay;
                end
                translation = (stim_texSpeed/sim_spaceResolution) * (tMove - stim_timeOffset - sim_endTime / 2);
                M_translate = imtranslate(M, [translation, 0],'FillValues',stim_meanLevel);
                M_rot = imrotate(M_translate, direction - 90);
                
                extend = floor(sim_dims(2:3)/2);
                x = (-extend(1):extend(1)) + round(size(M_rot,1)/2);
                y = (-extend(2):extend(2)) + round(size(M_rot,2)/2);
                y = y(1:sim_dims(3));
                x = x(1:sim_dims(2));
                stim_lightMatrix(ti,:,:) = M_rot(x,y);
            end
            
            
        case 'edgeGradient'
            % interpolate along points (x, y)
%             xUpper = [0 + stim_length, 0];
%             yUpper = [0, stim_contrast];
%             xLower = [0, 0 - stim_length];
%             yLower = [-stim_contrast, 0];
            
            pos = [-stim_singleSideLength, 0, 1e-10, stim_singleSideLength];
            cont = [0, -stim_contrast, stim_contrast, 0];
            
            xSel = X > -stim_width/2 & X < stim_width/2;
            
            
            % make pattern, centered
            s = zeros(sim_dims(2:3));
            
            s(xSel, :) = repmat(interp1(pos, cont, Y, 'linear', 0), sum(xSel), 1);
            
%             for yi = 1:sim_dims(3)
%                 y = Y(yi);
%                 if y > 0
%                     s(xSel, yi) = interp1(xUpper, yUpper, y, 'linear', 0);
%                 else
%                     s(xSel, yi) = interp1(xLower, yLower, y, 'linear', 0);
%                 end
%             end
            
%             s = s + 100;
            sr = imrotate(s, stim_angle, 'bilinear', 'crop');
%             sr = sr - 100;
%             sr(sr < 0) = 0;
            
            translation = stim_offset ./ sim_spaceResolution;
            sr = imtranslate(sr, translation, 'FillValues', 0);
            
            for ti = 1:sim_dims(1)
                t = T(ti);
                if t > stim_startTime && t < stim_startTime + stim_stimTime
                    
                    stim_lightMatrix(ti, :, :) = sr;
                    
                end
            end
    end
    
    stim_lightMatrix_byOption{optionIndex} = stim_lightMatrix;
end

%% plot movie of stimulus
% plotStimulus = 1
if plotStimulus
    optionIndex = 2;
    figure(101);
    set(gcf, 'Name','Stimulus Movie Display','NumberTitle','off');
    clf;
    stim_lightMatrix = stim_lightMatrix_byOption{optionIndex};
    for ti = 1:length(T)
        sim_light = squeeze(stim_lightMatrix(ti, :, :));
        clf
        plotSpatialData(mapX,mapY,sim_light);
        colormap gray
                caxis([-1,1])
%         caxis([0,1])
        colorbar
        title(sprintf('stimulus at %.3f sec, optionIndex %g', T(ti), optionIndex));
        axis tight
        drawnow
        pause(sim_timeStep)
        
    end
end