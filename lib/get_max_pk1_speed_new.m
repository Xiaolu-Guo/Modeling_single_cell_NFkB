function [max_pk1_speed,max_pk1_speed_frame] = get_max_pk1_speed_new(pk1_time, derivatives, FramesPerHour)
smoothed = zeros(size(derivatives));
for j = 1:size(smoothed, 1)
    smoothed(j, :) = smooth(derivatives(j, :), 'lowess');
end
pk1_frame = pk1_time *FramesPerHour + 1;
max_pk1_speed = nan(size(pk1_frame)); max_pk1_speed_frame = nan(size(pk1_frame));
for i=1:length(pk1_frame)
    if ~isnan(pk1_frame(i))
        [max_pk1_speed(i), max_pk1_speed_frame(i)] = nanmax(smoothed(i,1:pk1_frame(i)),[],2);
    end
end