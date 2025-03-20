function [pos_integral_features] = get_pos_integrals(integrals, FramesPerHour, endFrame)

integrals_pos = zeros(size(integrals));
for i = 1:size(integrals_pos,1)
    if integrals(i,1) > 0
        integrals_pos(i,1) = integrals(i,1);
    end
    for j = 2:size(integrals_pos,2)
        if integrals(i,j)>= integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1) + (integrals(i,j) - integrals(i,j-1));
        %elseif integrals(i,j) = integrals(i,j-1)
        elseif integrals(i,j)< integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1);
        end
    end
end

pos_integral_features.integrals_pos =integrals_pos; 

pos_integral_features.max_pos_integral = nanmax(integrals_pos,[],2);

pos_integral_features.time2HalfMaxPosIntegral = get_time_to_half_max_integral_new(integrals_pos(:,1:endFrame), FramesPerHour);
end
