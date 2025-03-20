function [] =  draw_sampling_traj_heatmap_2023_05(fig_save_path,data_filename)
% For Guo et al. Figure 5A, extended doses of each ligand simulated heatmaps
% 
% tested 05/12/2024, Matlab 2020a

%%

sample_data = 1; % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
vis_data_field = {'pred_mode_amp'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};
vis_cv_field = {'vis_cv'};

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

load(strcat(data_save_file_path_1,data_filename))
  
% [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter

if sample_data % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
    data_fields_names = {'model_sim'};
    
    for i_data = 1:length(data.(data_fields_names{1}))
        for i_data_name = 1:length(data_fields_names)
            data.(vis_data_field{1}){i_data} = data.(data_fields_names{i_data_name}){i_data}(1:9:end,:);
        end
    end
else
    data_fields_names = {'pred_mode_filter_nan'};
    i_data_name =1;
    data.(vis_data_field{i_data_name}) = data.(data_fields_names{i_data_name});
    
end

% savepath='./draft_figures/';
if 1
    for i_sti = 1:length(data.(data_fields_names{i_data_name}))
        
        data.(vis_cv_field{1}){i_sti} = std(data.(vis_data_field{1}){i_sti},[],2)./mean(data.(vis_data_field{1}){i_sti},2);
        [~,data.cv_order{i_sti}] = sort(data.(vis_cv_field{1}){i_sti},'descend');
        % [~,data.pred_mode_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2},'descend');
        % [~,data.exp_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2-1},'descend');
        
    end
end
filter_TNF = 0;
data.info_ligand{1}
plot_sampling_data_seperate_ligand_2023_05(data,vis_data_field,fig_save_path,filter_TNF)

