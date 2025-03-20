function [] =  draw_traj_heatmap_2023_03(data_save_file_path,fig_save_path) 

%This is for visualizing Monolix estimation

%% filepath


% clear all

    % data_save_file_path = '../SAEM_proj_2022_2/';
data_proj_num_vec = [2:6];
% data_proj_num_vec = 5;

% load('allsti_data_codon_0603.mat')

%%
% data = read_monolix_data_to_matrix_2023(data_proj_num_vec,data_save_file_path);

vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'};%,'sample'};
data_label = {'experiments','prediction'};%,'sample'};
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));

% [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter

% savepath='./draft_figures/';
if 1
for i_sti = 1:length(data.pred_mode)   
    data.pred_mode_fano{i_sti} = std(data.pred_mode_filter_nan{i_sti},[],2).^2./mean(data.pred_mode_filter_nan{i_sti},2);
    data.pred_mode_cv{i_sti} = std(data.pred_mode_filter_nan{i_sti},[],2)./mean(data.pred_mode_filter_nan{i_sti},2);
    [~,data.pred_mode_fano_order{i_sti}] = sort(data.pred_mode_fano{i_sti},'descend');
    [~,data.pred_mode_cv_order{i_sti}] = sort(data.pred_mode_cv{i_sti},'descend');
    [~,data.pred_mode_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2},'descend');
    [~,data.exp_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2-1},'descend');
end
end

vis_data_field = {'pred_mode_filter_nan'};
vis_data_field_order = {'pred_mode_cv_order'};
vis_data_field_order = {'pred_mode_fano_order'};
vis_data_field_order = {'pred_mode_osc_order'};
vis_data_field = {'exp_mode_filter_nan'};
vis_data_field_order = {'exp_osc_order'};

%% trajectories
vis_data_field = {'exp','pred_mode'};
% 
% for i_sti = length(data.(vis_data_field{1}))
%  [~, data.order{i_sti}] =sort(collect_feature_vects.Speed{2*i_sti-1},'descend');%OscVsNonOsc PeakAmplitude
% 
% end
% % 
% plot_data(data,vis_data_field)

plot_data_seperate_ligand_2023_03(data,vis_data_field,fig_save_path)
% cell_num_vec = [1:5;6:10;11:15;16:20;21:25;26:30]
% plot_trajectory(data,vis_data_field,cell_num_vec)
