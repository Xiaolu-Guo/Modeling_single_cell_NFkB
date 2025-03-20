function [] = Figure_S7A_single_cell_specifcity
% For Guo et al. Figure S7B: single-cell stimulus response specificity
% 
% tested 05/15/2024, Matlab 2020a

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

% Linearly interpolate between blue (at the start), white (in the middle), and red (at the end)
n =20;
blueToWhite = [linspace(0, 1, n/2)' linspace(0, 1, n/2)' ones(n/2, 1)];
whiteToRed = [ones(n/2, 1) linspace(1, 0, n/2)' linspace(1, 0, n/2)'];
customColormap = [blueToWhite; whiteToRed(2:end,:)];

load('../raw_data2023/Sim8_5_signle_ligand_codon_metric_r3.mat')

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

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
%[collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled

codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
% codon_list = {'TotalActivity'};

i_column = 1;
ligand_index = [1,5,3,2,4]; % reoder the stim

for i_ligand = 1:length(collect_feature_vects.info_ligand)
    for i_codon = 1:length(codon_list)
        codon_mat(:,i_column) = collect_feature_vects.(codon_list{i_codon}){ligand_index(i_ligand)}(:,:);
        i_column = i_column + 1;
    end
end

[rho,p] = corr(codon_mat, 'Type', 'Spearman');

figure(1)
h = heatmap(rho,'Colormap',customColormap);
caxis([-1,1])

YLabels = 1:size(codon_mat,2);
YCustomXLabels = string(YLabels);
YCustomXLabels(:) = " ";
h.YDisplayLabels = YCustomXLabels;
h.XDisplayLabels = YCustomXLabels;

saveas(gcf,strcat(fig_save_path,'SRS_codon_corr_wt'),'epsc')
close()

end


%%
function cluster_num = applyDBSCAN(matrix)
epsilon = 0.7; % Set the epsilon value
minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end