classdef mcm_Experiment < handle
    properties
        retina
        stimulus
        brain
        extent = 500

        cellType
       
        decodedLocation = [nan,nan]
    end
    
    
    methods
        
        function obj = mcm_Experiment()
            obj.retina = mcm_Retina(obj.extent);
            
%             fprintf('Init experiment\n');
            obj.stimulus = mcm_Stimulus();
            obj.brain = mcm_Brain(obj.retina);
        end
        
        
        function setup(obj)
            fprintf('set up experiment\n')
            obj.retina.setup()
            
        end
        
        function run(obj)
%             locs = obj.retina.locationsByCell;
            obj.retina.applyStimulus(obj.stimulus);
            obj.retina.applyNoise();

            decode = obj.brain.decode();
            obj.decodedLocation = decode;
            
            fprintf('Decoded position (%.3g, %.3g) from (%.3g, %.3g)\n', decode(1), decode(2), obj.stimulus.location(1), obj.stimulus.location(2));
         
        end 
        
        function show(obj, fignum, rfmodel, drawDecodedPosition)
            fprintf('showing\n');
            figure(fignum); clf;
            obj.stimulus.show()
            hold on
            obj.retina.show(1,rfmodel)
            colormap gray
%             colorbar
            caxis([0,1])
            
            xlim([-obj.extent, obj.extent])
            ylim([-obj.extent, obj.extent])
            axis square
            xlabel('Dorsal-Ventral (µm)');
            ylabel('Rostral-Nasal (µm)');
            
            if drawDecodedPosition
                c = obj.decodedLocation;
                l = 60;
                line(c(1) + [-l, l]/2, c(2) * [1,1], 'LineWidth', 2, 'Color', 'r');
                line(c(1) * [1,1], c(2) + [-l, l]/2, 'LineWidth', 2, 'Color', 'r');     
            end
            drawnow
            set(gcf,'color','w');
%             c = colorbar();
%             ylabel(c, 'Y Error')
        end
    end
    
end
