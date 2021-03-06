% if plotSubunitCurrents
%     figure(103);clf;
%     set(gcf, 'Name','Subunit signals','NumberTitle','off');
%     if useSubunits
%         axesSignalsBySubunit = tight_subplot(c_numSubunits, 2);
%     end
% end


sim_responseSubunitsCombinedByOption = {};

% Main stimulus change loop

if runInParallelPool
    parforArg = Inf;
else
    parforArg = 0;
end

% parfor (optionIndex = 1:stim_numOptions, parforArg)
for (optionIndex = 1:stim_numOptions)
    
%     fprintf('Running option %d of %d\n', optionIndex, stim_numOptions);
    
    stim_lightMatrix = stim_lightMatrix_byOption{optionIndex};
    
    %% Simulation starts here
    
    %% Run through time to calculate light input signals to each subunit
    %     s_lightSubunit = zeros(c_numSubunits, sim_dims(1));
    s_lightSubunit = {};
    area = (sim_dims(2) * sim_dims(3));
    
    prevLight = zeros(2);
    offsetBy = zeros(2);
    for ti = 1:length(T)
        sim_light = squeeze(stim_lightMatrix(ti, :, :));
        
        %% Calculate illumination for each subunit
        for vi = 1:e_numVoltages
            for ooi = 1:2
                
                for si = 1:c_numSubunits(vi)
                    lightInput = sum(sum(sim_light .* c_subunitRf{vi,ooi}(:,:,si))) / area;
                    
                    if makeLightInputMonotonic
                        % if it's going in the wrong direction, cancel it out
                        % by holding an offset value over time
                        if (ooi == 1 & lightInput < prevLight(vi,ooi)) | (ooi == 2 & lightInput > prevLight(vi,ooi))
                            % ON mode and decrease or OFF mode and increase
                            offsetBy(vi,ooi) = prevLight(vi,ooi) - lightInput;
                        end
                    end
                    lightOrigOffset = lightInput + offsetBy(vi,ooi);
                    
%                     if ooi == 1 & vi == 1
%                         fprintf('prev %g,  input %g, holding offset: %g\n', prevLight(vi,ooi), lightInput, offsetBy(vi,ooi))
%                     end
%                     
                    prevLight(vi,ooi) = lightInput;
                    
                    s_lightSubunit{vi,ooi}(si,ti) = lightOrigOffset;
                    
                end
            end
        end
    end
    
    %% temporal filter each subunit individually
    
    
    
    %         a = sim_lightSubunit(si,:);
    %         d = diff(a);
    %         d(end+1) = 0;
    %         d(d < 0) = 0;
    %         lightOnNess = cumsum(d);
    
    %%      for each subunit & voltage, do a bunch of signal processing
    
    s_responseSubunitScaledByRf = {};
    sim_responseSubunitsCombined = [];
    
    if plotSubunitCurrents
        figure(143);clf;
        [~,hg] = tight_subplot(2,2,0.05, 0.05);
    end
    
    for vi = 1:e_numVoltages
        for ooi = 1:2

            % loop through subunits
            s_responseSubunit = [];
            s_responseSubunitScaledByRf{vi,ooi} = [];

            for si = 1:c_numSubunits(vi,ooi)

                % linear convolution
                temporalFiltered = conv(s_lightSubunit{vi,ooi}(si,:), c_filtersResampled{vi,ooi});
                temporalFiltered = temporalFiltered(1:length(s_lightSubunit{vi,ooi}(si,:)));

                % nonlinear effects
                %             sel = [];
                if e_voltages(vi) < 0 % excitation
                    sel = temporalFiltered < 0;
%                     s = -1;
                else
                    sel = temporalFiltered > 0;
%                     s = 1;
                end

                % rectify and saturate (replace with a hill function)
                rectified = temporalFiltered;
                rectifiedSignalMultiplier = 0.2;
                if vi == 1
                    if ooi == 1 && rectifyExOn
                        rectified(sel) = rectified(sel) * rectifiedSignalMultiplier;
                    end
                    if ooi == 2 && rectifyExOff
%                         rectified(sel) = 0;
                        rectified(sel) = rectified(sel) * rectifiedSignalMultiplier;
                    end
                elseif vi == 2
                    rectified(sel) = 0;
                end

                rectThresh = .005;
                rectMult = 1;
                %             saturThresh = .04;
                sel = abs(rectified) > rectThresh;
                nonlin = rectified;
                nonlin(sel) = nonlin(sel) * rectMult;
                %             nonlin(abs(nonlin) > saturThresh) = sign(nonlin(abs(nonlin) > saturThresh)) * saturThresh;

                %             filterTimeConstant = 0.17;
                %             filterLength = 0.4; %sec
                %             tFilt = 0:sim_timeStep:filterLength;
                %             filter_subunitTemporalDecay = exp(-tFilt/tau);
                %             filter_subunitTemporalDecay = filter_subunitTemporalDecay ./ sum(filter_subunitTemporalDecay);
                %             decayed = conv(nonlin, filter_subunitTemporalDecay, 'same');
                %

                s_responseSubunit(si,:) = nonlin(1:sim_dims(1));

                warning('off','MATLAB:legend:IgnoringExtraEntries')

                if plotSubunitCurrents
                    % plot individual subunit inputs and outputs

                    axes(hg(vi, ooi));
                    %                 axes(axesSignalsBySubunit(((si - 1) * 2 + vi)))
                    plot(normg(s_lightSubunit{vi,ooi}(si,:)));
                    hold on
                    %             plot(normg(lightOnNess))
                    plot(normg(c_filtersResampled{vi,ooi}) - 0.5)
                    %                 plot(normg(filter_subunitTemporalDecay));
                    plot(normg(temporalFiltered));
                    plot(normg(rectified));
                    plot(normg(nonlin));
                    %                 plot(normg(decayed));
                    a = {'on','off'};
                    title(sprintf('subunit %d light convolved with filter v = %d %s', si, e_voltages(vi), a{ooi}))
                    hold off
                    legend('light','temporalFilter', 'filtered', 'rectified', 'nonlinear')
                    drawnow
%                     if ooi + vi == 4
%                         pause
%                     end
                end


                % Multiply subunit response by RF strength (connection subunit to RGC)

                strength = s_subunitStrength{vi,ooi}(si);
                s_responseSubunitScaledByRf{vi,ooi}(si,:) = strength * s_responseSubunit(si,:);

            end
        % combine current across subunits for this voltage and polarity
        sim_responseSubunitsCombined(vi,ooi,:) = sum(s_responseSubunitScaledByRf{vi,ooi}, 1);    
        end
    end
    
    sim_responseSubunitsCombinedByOption{optionIndex} = sim_responseSubunitsCombined;
    
end % end of options and response generation
