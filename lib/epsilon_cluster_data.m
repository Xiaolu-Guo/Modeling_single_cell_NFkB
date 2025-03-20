
function clusterResults = epsilon_cluster_data(data,metrics,epsilon)
for i_data = 1:length(data.model_sim)
    
    data_info.info_ligand{i_data} = data.info_ligand{i_data};
    data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
    data_info.data_label{i_data} = 'SamplingComb';
    
end
metrics_cal = cell(1,length(metrics));
metric_names = fieldnames(metrics{1});

for i_metric_index = 1:length(metrics)
    for i_metric_name = 1:length(metric_names)
        metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
    end
end

%[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
[collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled

codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
%codon_list = {'TotalActivity'};

%% draw single cells
clear codon_single_cell
codon_single_cell=cell(2,1);
for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
    for i_ligand = 1:length(collect_feature_vects.info_ligand)
        for i_codon = 1:length(codon_list)
            codon_single_cell{i_cell}(i_ligand,i_codon) = collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
        end
    end
end

clusterResults = cell2mat(cellfun(@(matrix) applyDBSCAN(matrix,epsilon), codon_single_cell, 'UniformOutput', false));

end

function cluster_num = applyDBSCAN(matrix,epsilon)

% epsilon = 1; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end