% goal: plot heatmap of supriya's data, ade's data, fitting data, sampling data
% in the order of total activitiy, freq, osc power, etc.
% need to set to 1 to draw save figs 
%% notes
% Feature_Name	specifier	Codeword
% duration	2	Duration
% time2HalfMaxPosIntegral	-1	EarlyVsLate
% oscpower	0	OscVsNonOsc
% max_value	0	PeakAmplitude
% pos_pk1_amp	0	PeakAmplitude
% max_pos_pk1_speed	0	Speed
% pos_pk1_time	-1	Speed
% derivatives	2	Speed
% max_pos_integral	0	TotalActivity


% TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
% Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
% CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
% LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
% PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15

% goal: supriya's data, ade's data, fitting data, and sampled data for
% benchamarking
% compare the prediction of smapled data with supriya's data

% TNF 10ng/mL Pam3CSK 100ng/mL supriya: 03 index: 3
% TNF 10ng/mL CpG 100nM supriya: 03 index: 5
% TNF 10ng/mL LPS 10ng/mL supriya: 04 index: 8
% TNF 10ng/mL PolyIC 100ug/mL supriya: 04 index: 9
% Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
% Pam3CSK 100ng/mL LPS 10ng/mL supriya: 23 index: 26
% Pam3CSK 100ng/mL PolyIC 100ug/mL supriya: 13 index: 14
% CpG 100nM LPS 10ng/mL supriya: 23 index: 25
% CpG 100nM PolyIC 100ug/mL supriya: 21 index: 18
% LPS 10ng/mL PolyIC 100ug/mL supriya: 21 index: 19


%% initalization 
% clear all
addpath('./lib/')
addpath('./src/')
addpath('./bin/')

fig_save_path = '../SubFigures2023/';

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';


%% plot heatmap of ade's data
% Figure S4A: tested 05/10/2024
if 1
    clear data_to_draw data metrics
    
    % Ade's data:
    
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    data_index = [3;18;11;6;16];
    
    vis_data_field = {'exp_ade'};
    
    order_name_vec = {'max_pos_integral'};
    data_to_draw.(vis_data_field{1}) = data.exp_mode_filter_nan(data_index);
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);
    
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)*2-1}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end


%% plot heatmap of sampling 
% Figure S3A: tested 05/10/2024
if 1
    clear data_to_draw data metrics
    
    % sampling  data
    data_filename = 'Sim5_codon_all5dose_metric.mat';
    load(strcat(data_save_file_path_1,data_filename))
    data_index = [13;12;11;14;15];
    vis_data_field = {'sampling'};
    order_name_vec = {'max_pos_integral'};

    
    metric_fields_name = fieldnames(metrics{1});
    for i_metric = 1:length(metrics)
        
        for i_metric_fields = 1:length(metric_fields_name)
            metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
        end
    end
    metrics = metrics_sample;
    
    clear data_to_draw
    data_to_draw.(vis_data_field{1}) = data.model_sim(data_index);
    for i_data = 1:length(data_to_draw.(vis_data_field{1}))
        data_to_draw.(vis_data_field{1}){i_data} = data_to_draw.(vis_data_field{1}){i_data}(1:9:end,:);
    end
    
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);

    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end

