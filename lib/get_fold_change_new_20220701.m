function[output] = get_fold_change_new_20220701(time_series_no_base_ded, baseline,FramesPerHour)

%calculates max fold change in time_series matrix (non-baseline deducted)
output.fold_change = time_series_no_base_ded(:,1:end)./baseline;
output.max_fold_change = nanmax(output.fold_change,[],2);
[output.max_value, idx_max] = nanmax(time_series_no_base_ded - baseline, [], 2);

index_mat = ones(size(time_series_no_base_ded,1),1) * 1:size(time_series_no_base_ded,2) <= idx_max;

halfMaxIntegral = output.max_value/2;

distances = abs(time_series_no_base_ded - baseline - halfMaxIntegral);
distances = distances .* index_mat + 100000* ~index_mat; 

[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
output.time2HalfMaxValue = (idx-1)/FramesPerHour;

end