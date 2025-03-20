function  time_to_half_max_integral=get_time_to_half_max_integral_new(integrals, FramesPerHour)

halfMaxIntegral = nanmax(integrals,[],2)/2;
 
 distances = abs(integrals- halfMaxIntegral);
 [~, idx] = nanmin(distances,[],2);
 idx(idx==1) = NaN;
 time_to_half_max_integral = (idx-1)/FramesPerHour;

end

