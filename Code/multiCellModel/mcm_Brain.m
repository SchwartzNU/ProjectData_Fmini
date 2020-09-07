classdef mcm_Brain < handle

    properties
        retina
        decodeMethod
        decodedCells
    end
    
    methods
        
        function obj = mcm_Brain(retina)
            obj.retina = retina;
            obj.decodeMethod = 'com';
        end
        
        function location = decode(obj)
            locations = obj.retina.locationsByCell;
            responses = obj.retina.responsesByCellWithNoise;
            
            switch obj.decodeMethod
                case 'com'
                    exclude = responses <= max(responses) * 0.3;
                    
                    fprintf('using %g of %g = %g responses\n',sum(~exclude), length(responses), sum(~exclude)/length(responses));
                    
                    
                    responses(exclude) = [];
                    locations(exclude,:) = [];
%                     responses = responses .^ 2;
                    responses = responses ./ sum(responses);

                    x = sum(locations(:,1) .* responses);
                    y = sum(locations(:,2) .* responses);
                    
                    location = [x,y];
                    obj.decodedCells = ~exclude;
                    
                case 'peak'
                    [~,i] = max(responses);
                    location = locations(i,:);
            end
        end
    end
    
    
end