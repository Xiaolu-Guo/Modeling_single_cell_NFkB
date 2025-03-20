function [] =  save_sampling_20doses_for_channel_capacity_202408(data_filename)

%This is for visualizing Monolix estimation

%% filepath


% data_save_file_path = '../SAEM_proj_2022_2/';
% data_proj_num_vec = 5;

% load('allsti_data_codon_0603.mat')

%%

sample_data = 1; % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
vis_data_field = {'pred_mode_amp'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};
vis_cv_field = {'vis_cv'};

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';


load(strcat(data_save_file_path_1,data_filename))


metrics_new = cell(1,2);
metric_names = fieldnames(metrics{1});
for i_metric_name = 1:length(metric_names)
    for i_metric = 1:length(metrics)
        
        metrics_new{i_metric}.(metric_names{i_metric_name}) = metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
        data_info.info_ligand{i_metric} = collect_feature_vects.info_ligand{i_metric};
        data_info.info_dose_str{i_metric} = collect_feature_vects.info_dose_str{i_metric};
        data_info.data_label{i_metric} = 'sampling';
    end
end

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_new); %,  parameter
collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);


%% save the signaling coodn data to mutual info format
clear nfkb
load('./example_data_format/mutual_info_cal_data_example.mat')

codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};


sti_all = cellfun(@strcat,collect_feature_vects.info_ligand,collect_feature_vects.info_dose_str,'UniformOutput', false) ;
ligand_inf = collect_feature_vects.info_ligand{1};

for i_sti = 1:length(sti_all)
    nfkb(i_sti).sc_metrics = struct();
    
    for i_codon =1:length(codon_list)
        nfkb(i_sti).id = sti_all{i_sti};
        nfkb(i_sti).ids = sti_all;
        nfkb(i_sti).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){i_sti};
        
    end
end


if 1
    switch ligand_inf
        case 'TNF'
            index_3_doses =[1,7,11,15];% d7,d11,d15
        case 'CpG' %: d9,d11,d13
            index_3_doses = [1,9,11,13];
        case 'LPS' %: d11,d13,d15
            index_3_doses = [1,11,13,15];
        case 'Pam3CSK' %: d7,d11,d15
            index_3_doses = [1,7,11,15];
        case 'PolyIC' %: d11,d13,d15
            index_3_doses = [1,11,13,15];
    end
    
    for i_index_3_doses = 1:length(index_3_doses)
        nfkb_4doses(i_index_3_doses) = nfkb(index_3_doses(i_index_3_doses));
        nfkb_4doses(i_index_3_doses).ids = nfkb(index_3_doses(i_index_3_doses)).ids(index_3_doses);
    end
    
    nfkb = nfkb_4doses;
    save(strcat('mutual_info_format_codon_',ligand_inf,'_3doses_zerodose_20230707.mat'),'nfkb')
    
end


