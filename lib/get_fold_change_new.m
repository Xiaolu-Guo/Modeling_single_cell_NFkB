function[output] = get_fold_change_new(time_series_no_base_ded, baseline)

%calculates max fold change in time_series matrix (non-baseline deducted)
output.fold_change = time_series_no_base_ded(:,1:end)./baseline;
output.max_fold_change = nanmax(output.fold_change,[],2);
output.max_value = nanmax(time_series_no_base_ded - baseline, [], 2);
output.time2
end
