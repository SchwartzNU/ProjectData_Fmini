% shapeModelSetup


% sim_endTime = 2.4;
% sim_timeStep = 0.020;
% sim_spaceResolution = 2; % um per point
% s_sidelength = 300;%max(cell_rfPositions);
c_extent = 0; % start and make in loop:

switch stim_mode
    case 'flashedSpot'
        sim_endTime = 1.4;
        sim_timeStep = 0.02;
        sim_spaceResolution = 6; % um per point
        s_sidelength = 1200;%max(cell_rfPositions);

    case 'movingBar'
        sim_endTime = 2.8;
        sim_timeStep = 0.001;
        sim_spaceResolution = 1; % um per point
        s_sidelength = 400;%max(cell_rfPositions);
        
    case 'driftingGrating'
        sim_endTime = 2.0;
        sim_timeStep = 0.004;
        sim_spaceResolution = 2; % um per point
        s_sidelength = 600;%max(cell_rfPositions);        
        
    case 'edgeGradient'
        sim_endTime = .3;
        sim_timeStep = 0.01;
        sim_spaceResolution = 1; % um per point
        s_sidelength = 500;%max(cell_rfPositions);

end
        
% subunit locations, to generate positions
% c_subunitSpacing = [20 20; 35 35];
% c_subunit2SigmaWidth = [40 40; 60 60;];
% c_subunit2SigmaWidth_surround = [80 80; 120 120];
% c_subunitSurroundRatio = [0.15 0.15; 0.0 0.0];

% paramValues(paramSetIndex,col_rfOffset)
% positionOffsetByVoltageDim = [6, 8; 0, 0];
% positionOffsetByVoltageDim = [paramValues(paramSetIndex,col_rfOffset), 0; 0, 0];
% positionOffsetByVoltageOnoffDim = zeros(2,2,2);
% positionOffsetByVoltageOnoffDim(1, :, 2) = 7;

% gaussianSigmaByVoltageOnoffDim = 40*ones(2,2,2);

% good F mini On offsets for edge:
positionOffsetByVoltageOnoffDim = zeros(2,2,2);
% positionOffsetByVoltageOnoffDim(:,1,2) = 0; % ON dorsal
% positionOffsetByVoltageOnoffDim(:,2,2) = -0; % OFF ventral
% positionOffsetByVoltageOnoffDim(:,1,2) = paramValues(paramSetIndex, col_rfOffset) / 2; % ON dorsal
% positionOffsetByVoltageOnoffDim(:,2,2) = -paramValues(paramSetIndex, col_rfOffset) / 2; % OFF ventral

% center
gaussianSigmaByVoltageOnoffDim = zeros(2,2,2);
% gaussianSigmaByVoltageOnoffDim(1,:,:) = 50; % 50 is the one that fits SMS

% surround
gaussianSigmaByVoltageOnoffDim(2,:,:) = 120;

% gaussianSigmaByVoltageOnoffDim(1,1,:) = [50, 30]; % on ex
% gaussianSigmaByVoltageOnoffDim(1,2,:) = [50, 30]; % off ex


% good F mini On offsets for edge:
offsetRF = 1;
D = rfSizes(sizeI);
if offsetRF
    F = rfOffsets(offsetI);
    % vertical spatial offset
%     positionOffsetByVoltageOnoffDim(:,1,2) = (F * D)/2;
%     positionOffsetByVoltageOnoffDim(:,2,2) = -(F * D)/2;%-28; % OFF Offset
    
    % vertical size shortening
%     rfY = D * (1 - F);
    rfY = paramValues(paramSetIndex, col_rfAspect) * D
    gaussianSigmaByVoltageOnoffDim(1,:,:) = D; % horizontal (both) size
    gaussianSigmaByVoltageOnoffDim(1,1,2) = rfY; % on ex  (x,y)
    gaussianSigmaByVoltageOnoffDim(1,2,2) = rfY; % off ex (x,y)
end


% generate RF map for EX and IN
% import completed maps 

if useRealRf
    load(sprintf('rfmaps_%s_%s.mat', cellName, acName));
    ephys_data_raw = data;
else
    voltages = [-60;20];
    intensities = [1];
end

%% Setup cell data from ephys

e_positions = {};
e_voltages = sort(voltages);
e_numVoltages = length(e_voltages);
e_intensities = intensities;
e_numIntensities = length(intensities);
clear voltages
clear intensities
s_voltageLegend = {};
for vi = 1:e_numVoltages
    s_voltageLegend{vi} = num2str(e_voltages(vi));
end
s_voltageLegend = {'ex','in'};
s_voltageLegend{end+1} = 'Combined';

T = 0:sim_timeStep:sim_endTime;

% dims for: time, X, Y
sim_dims = round([length(T), s_sidelength / sim_spaceResolution, s_sidelength / sim_spaceResolution]);
e_map = nan * zeros(sim_dims(2), sim_dims(3), e_numVoltages, 2);

ii = 1; % just use first intensity for now
if useRealRf
    for vi = 1:e_numVoltages
        e_vals(vi,:) = ephys_data_raw(vi, ii, 2);
        pos = ephys_data_raw{vi, ii, 1:2};
        e_positions{vi, ii} = pos; %#ok<*SAGROW>
    end
end



X = linspace(-0.5 * sim_dims(2) * sim_spaceResolution, 0.5 * sim_dims(2) * sim_spaceResolution, sim_dims(2));
Y = linspace(-0.5 * sim_dims(3) * sim_spaceResolution, 0.5 * sim_dims(3) * sim_spaceResolution, sim_dims(3));
[mapY, mapX] = meshgrid(Y,X);
distanceFromCenter = sqrt(mapY.^2 + mapX.^2);

                    %  ex  in
% shiftsByDimVoltage = [-30,-30;  % x
%                       -30,-30]; % y
% shiftsByDim = analysisData.positionOffset;
% positionOffset = paramValues{paramSetIndex,col_positionOffset};


% Set up spatial RFs for voltages and OO
for vi = 1:e_numVoltages
    for ooi = 1:2
    
    %     c = griddata(e_positions{vi, ii}(:,1), e_positions{vi, ii}(:,2), e_vals{vi,ii,:}, mapX, mapY);
    %     e_map(:,:,vi) = c;

        % add null corners to ground the spatial map at edges
        if useRealRf
            positions = e_positions{vi, ii};
            positions = bsxfun(@plus, positions, [positionOffsetByVoltageOnoffDim(vi,ooi,1),positionOffsetByVoltageOnoffDim(vi,ooi,2)]);
            vals = e_vals{vi,ii,:};
        %     positions = vertcat(positions, [X(1),Y(1);X(end),Y(1);X(end),Y(end);X(1),Y(end)]);
        %     vals = vertcat(vals, [0,0,0,0]');
            F = scatteredInterpolant(positions(:,1), positions(:,2), vals,...
                'linear','none');
            m_rf = F(mapX, mapY) * sign(e_voltages(vi));
        else
            % Simple gaussian RF approximation
            d = sqrt((mapY-positionOffsetByVoltageOnoffDim(vi,ooi,2)).^2 / gaussianSigmaByVoltageOnoffDim(vi,ooi,2)^2 + ...
                (mapX-positionOffsetByVoltageOnoffDim(vi,ooi,1)).^2 / gaussianSigmaByVoltageOnoffDim(vi,ooi,1)^2);
            m_center = exp(-d.^2);
            
            s_enableSurroundGaussian = 1;
            if s_enableSurroundGaussian & vi == 1
                d = sqrt((mapY-positionOffsetByVoltageOnoffDim(vi,ooi,2)).^2 / gaussianSigmaByVoltageOnoffDim(vi,ooi,2)^2 + ...
                    (mapX-positionOffsetByVoltageOnoffDim(vi,ooi,1)).^2 / gaussianSigmaByVoltageOnoffDim(vi,ooi,1)^2);
                m_surround = .2 * exp(-d.^2 / 4);
                m_rf = m_center - m_surround;
            else
                m_rf = m_center;
            end
            
        end

        m_rf(isnan(m_rf)) = 0;
%         m(m < 0) = 0;
%         m_rf = m_rf ./ max(m_rf(:)); % normalize to peak 1
        m_rf = m_rf ./ sum(m_rf(:)) * 1000; % normalize sum to 1000
        e_map(:,:,vi,ooi) = m_rf;
    %     e_map(:,:,vi) = e_map(:,:,vi) - min(min(e_map(:,:,vi)));

        c_extent = max(c_extent, max(distanceFromCenter(m_rf > 0)));

%         if plotSpatialGraphs
%             axes(axesSpatialData((vi-1)*3+1))
%     %         imgDisplay(X,Y,e_map(:,:,vi))
%             plotSpatialData(mapX,mapY,e_map(:,:,vi))
%             title(s_voltageLegend{vi});
%             colormap parula
%             colorbar
%             axis equal
%             axis tight
% %             title('cell RF')
%     %         xlabel('µm')
%     %         ylabel('µm')
%         %     surface(mapX, mapY, zeros(size(mapX)), c)
%         end
    end
end

%% Subunit code

if useSubunits
    
    c_subunitSigma = c_subunit2SigmaWidth / 2;
    c_subunitSigma_surround = c_subunit2SigmaWidth_surround / 2;
    c_subunitCenters = {};
    for vi = 1:2
        for ooi = 1:2
            c_subunitCenters{vi,ooi} = generatePositions('triangular', [c_extent, c_subunitSpacing(vi,ooi), 0]);
            c_numSubunits(vi,ooi) = size(c_subunitCenters{vi,ooi},1);
        end
    end

    % subunit RF profile, using gaussian w/ set radius (function)
    c_subunitRf = {};
    for vi = 1:2
        for ooi = 1:2
            c_subunitRf{vi,ooi} = zeros(sim_dims(2), sim_dims(3), c_numSubunits(vi,ooi));
            for si = 1:c_numSubunits(vi)
                center = c_subunitCenters{vi,ooi}(si,:);
                dmap = (mapX - center(1)).^2 + (mapY - center(2)).^2; % no sqrt, so
                rf_c = exp(-(dmap / (2 * c_subunitSigma(vi,ooi) .^ 2))); % no square
                rf_s = exp(-(dmap / (2 * c_subunitSigma_surround(vi,ooi) .^ 2))); % no square

                rf = rf_c - c_subunitSurroundRatio(vi,ooi) * rf_s;
                rf = rf ./ max(rf(:));
                c_subunitRf{vi,ooi}(:,:,si) = rf;
            end
        end
    end

    % calculate connection strength for each subunit, for each voltage
    s_subunitStrength = {};
    for vi = 1:e_numVoltages
        for ooi=1:2
            s_subunitStrength{vi,ooi} = zeros(c_numSubunits(vi,ooi),1);
            for si = 1:c_numSubunits(vi,ooi)

                rfmap = e_map(:,:,vi,ooi);
                sumap = c_subunitRf{vi,ooi}(:,:,si);
                [~,I] = max(sumap(:));
                [x,y] = ind2sub([sim_dims(2), sim_dims(3)], I);

                s_subunitStrength{vi,ooi}(si) = rfmap(x,y);

        %         todo: change it to a regression between map and each subunit as a predictor
        %         s_subunitStrength{vi}(si) = sum(rfmap(:) ./ sumap(:));
            end
        end
    end

    % remove unconnected subunits
    for vi = 1:e_numVoltages
        for ooi=1:2
            nullSubunits = s_subunitStrength{vi,ooi} < eps+.1;
            c_subunitRf{vi,ooi}(:,:,nullSubunits) = [];
            s_subunitStrength{vi,ooi}(nullSubunits) = [];
            c_subunitCenters{vi,ooi}(nullSubunits',:) = [];
            c_numSubunits(vi,ooi) = size(s_subunitStrength{vi,ooi},1);
        end
    end
    
else % single complete subunit
    c_subunitRf = {};
    s_subunitStrength = {};
    c_numSubunits = [];
    for vi = 1:2
        for ooi = 1:2
            c_subunitRf{vi,ooi}(:,:,1) = e_map(:,:,vi,ooi);
            s_subunitStrength{vi,ooi}(1) = 1;
            c_numSubunits(vi,ooi) = 1;
        end
    end
       
end

%% plot the spatial RF graphs

os = struct();
os.X = X;
os.Y = Y;
oo = {'ON','OFF'};
vv = {'EX','IN'};

if plotSpatialRF
    figure(90);clf;
    set(gcf, 'Name','Spatial Graphs','NumberTitle','off');
    [~, axesSpatialData] = tight_subplot(e_numVoltages, 3, .05, .15);
    for vi = 1:2%:e_numVoltages
        for ooi = 1:2%:2
            axes(axesSpatialData(vi, ooi))
            d = zeros(sim_dims(2), sim_dims(3));
            for si = 1:c_numSubunits(vi,ooi)
                d = d + c_subunitRf{vi, ooi}(:,:,si) * s_subunitStrength{vi, ooi}(si);
            end
            plotSpatialData(mapX,mapY,d)
            axis equal
            axis tight
            a = {'on','off'};
            title(a{ooi})
    %         title('all subunits scaled by maps')
            hold on
            
%             e = ellipse(60, 30, 0, 0, 0);
%             set(e, 'Color', 'm')
%             e.LineWidth = 4;
%             
%             plot(0,0,'.','MarkerSize',50,'Color','m');
            
            hold off
            
    %         xlabel('µm')
    %         ylabel('µm')
            % plot points at the centers of subunits
            if c_numSubunits(vi,ooi) > 1
                plot(c_subunitCenters{vi}(:,1), c_subunitCenters{vi}(:,2),'r.')
            end
            
            nam = sprintf('rf_%s_%s', oo{ooi}, vv{vi});
            os.(nam) = d';
        end
    end
%     xlim([-150,150])
%     ylim([-150,150])
    drawnow
end
exportStructToHDF5(os, 'igorExport/model_maps.h5', 'model_rfmaps');

