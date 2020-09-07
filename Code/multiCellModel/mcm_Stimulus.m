classdef mcm_Stimulus < handle
    
    properties
        meanLevel = 0;
        contrast = 1;
        location = [0,0]; % X, Y
        width = 150;
        singleSideLength = 150;
        angle = 0;
    end
    
    methods
        function show(obj)
            x = [-obj.width/2*cosd(obj.angle), obj.width/2*cosd(obj.angle)];
            y = sind(obj.angle) * -obj.width/2 + obj.location(2);
            x = x + obj.location(1);
%             plot(x,y, 'k', 'LineWidth', 2)
            ax = gca();
            ax.Color = [.5, .5, .5];
            
            x0 = x(1);
            y0 = y(1);

            xk = -sind(obj.angle)*obj.singleSideLength;
            xj = cosd(obj.angle)*obj.width;
            pX = [x0, x0+xj, x0+xj+xk, x0+xk];
            
            yk = cosd(obj.angle)*obj.singleSideLength;
            yj = sind(obj.angle)*obj.width;
            pY = [y0, y0+yj, y0+yk+yj, y0+yk];
            pC = [1, 1, .5, .5]';
            pUpper = patch('XData', pX, 'YData', pY, 'FaceColor', 'interp', 'FaceVertexCData', pC, 'EdgeColor', 'none');
            
            
            pX = [x0, x0+xj, x0+xj-xk, x0-xk];
            pY = [y0, y0+yj, y0-yk+yj, y0-yk];
            pC = [0, 0, .5, .5]';
            pLower = patch('XData', pX, 'YData', pY, 'FaceColor', 'interp', 'FaceVertexCData', pC, 'EdgeColor', 'none');            
            
        end
        
        function yOffset = getOffset(obj, cellLocation)
            yOffset = sind(obj.angle) * cellLocation(1) + obj.location;
        end
    end
end