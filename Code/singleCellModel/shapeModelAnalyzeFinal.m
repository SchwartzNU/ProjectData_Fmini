

%% display DSI over parameters

if plotSummaryVarOverParamsets
    if numParamSets >= 3
        figure(901);clf;
        %     clf;

        if numParams == 1
            %         d1 = linspace(min(paramValues(:,1)), max(paramValues(:,1)), 10);
            %         m = interp1(paramValues(:,1), summaryVarByParamSet, d1);
            %         plot(d1, m)
            switch stim_mode
                case 'flashedSpot'

                case 'movingBar'
                    plot(paramValues(:, col_barSpeed), summaryVarByParamSet)
                    hold on
                    plot(paramValues(:, col_barSpeed), [0.25, 0.12, 0.07, 0.07])
                    legend('Model','FmON','Control')
                    xlabel(paramColumnNames{1})
                    ylabel('DSI/OSI')

                case 'driftingGrating'
                    plot(paramValues(:, col_rfOffset), summaryVarByParamSet(:,1))
                    hold on
                    plot(paramValues(:, col_rfOffset), summaryVarByParamSet(:,2))
                    ylim([0,1])
                    
                    line(0.7+[0,0], [0,1])

                    legend('DSI','OSI')
                    xlabel(paramColumnNames{1})
                    ylabel('DSI/OSI')
            end

            %% Plot FmON cells
            %         hold on
            %         plot(abs(dtab{selectWfdsOn,'spatial_exin_offset'}), dtab{selectWfdsOn,'best_DSI_sp'}, '.', 'MarkerSize', 20)
            %         plot(abs(dtab{selectWfdsOff,'spatial_exin_offset'}), dtab{selectWfdsOff,'best_DSI_sp'}, '.', 'MarkerSize', 20)
            %         plot(abs(dtab{selectControl,'spatial_exin_offset'}), dtab{selectControl,'best_DSI_sp'}, '.', 'MarkerSize', 20)
            %         hold off




        elseif numParams == 2
            switch paramResponsePlotMode2d
                case 'lines'
                    speeds = unique(paramValues(:,1));
                    offsets = unique(paramValues(:,2));
                    x = speeds;


                    for offi = 1:length(offsets)
                        curOffset = offsets(offi);

                        % make a speed, DSI plot
                        y = [];
                        for spi = 1:length(speeds)
                            curSpeed = speeds(spi);
                            y(spi) = summaryVarByParamSet(paramValues(:,1) == curSpeed & paramValues(:,2) == curOffset, 2);
                        end
                        plot(x,y)
                        hold on
                        outputStruct.(sprintf('dsi_offset%g',curOffset)) = y;
                    end

                    plot(speeds, [0.25, 0.12, 0.07, 0.07], '--', 'Color', 'k', 'LineWidth',3)
                    plot(speeds, [0.2 0.35 0.12 0.05], '-.', 'Color', 'g', 'LineWidth',2)

                    legend(vertcat(cellstr(num2str(offsets, 'offset=%-d')), {'real data'; 'sanes data'}))

                case 'spatial'
                    d1 = linspace(min(paramValues(:,1)), max(paramValues(:,1)), 10);
                    d2 = linspace(min(paramValues(:,2)), max(paramValues(:,2)), 10);

                    [d1q,d2q] = meshgrid(d1, d2);
                    c = griddata(paramValues(:,1), paramValues(:,2), summaryVarByParamSet(:,1), d1q, d2q);
                    surface(d1q, d2q, zeros(size(d1q)), c)
                    hold on
                    %         plot(paramValues(:,1), paramValues(:,2), '.', 'MarkerSize', 40)
                    hold off
                    xlabel(paramColumnNames{1})
                    ylabel(paramColumnNames{2})
                    colorbar
            end

        elseif numParams == 3
            valueMatrix = [];
            for i = 1:size(paramIndices,1)
                pi = paramIndices(i,:);
                valueMatrix(pi(1), pi(2), pi(3)) = summaryVarByParamSet(i, 1);
            end

            [paramGridPoints1,paramGridPoints2,paramGridPoints3] = meshgrid(p1, p2, p3);

            slice(paramGridPoints1, paramGridPoints2, paramGridPoints3, valueMatrix,0,0,[0,180])
            xlabel('X offset')
            ylabel('Y offset')
            zlabel('angle')
            
            colormap parula
            colorbar
        end
    end
end

%% play with 3d output

return

valueMatrix = [];
for i = 1:size(paramIndices,1)
    pi = paramIndices(i,:);
    valueMatrix(pi(1), pi(2), pi(3)) = summaryVarByParamSet(i, 1);
end

[paramGridPoints1,paramGridPoints2,paramGridPoints3] = meshgrid(p1, p2, p3);
% 
% slice(paramGridPoints1, paramGridPoints2, paramGridPoints3, valueMatrix,0,0,[0,180])
% xlabel('X offset')
% ylabel('Y offset')
% zlabel('angle')
% out = interp3(paramGridPoints1, paramGridPoints2, paramGridPoints3, valueMatrix, 1,1,1);

fname = [SAM_RESEARCH sprintf('datasets/shapeModelOutputMatrix_surr_dim%.2g_off%.3g.mat', rfSizes(sizeI), rfOffsets(offsetI))];
save(fname,'paramGridPoints1','paramGridPoints2','paramGridPoints3','valueMatrix');
fprintf('Saved 3D output data to %s\n',fname)
% out = griddedInterpolant(paramValues, valueMatrix);
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
%     nonlinDsiByParamSet(ps,1) = calcDsi(anglesRads, nonlin);
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
% plot(cell2mat(paramValues)', dsiByParamSet);
% plot(cell2mat(paramValues)', nonlinDsiByParamSet);
%
% legend('sanes','model','modelNonlin')
% hold off