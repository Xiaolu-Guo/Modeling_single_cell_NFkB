function [] = ppt_ParameterEstimation_draw_para_distribution(data_save_file_path,fig_save_path) 

% data_proj_num_vec=[18,25:28];%TNF,LPS [27];%


%% para_meter_distribution
% 
% sti_vec = {'TNF','LPS','CpG','polyIC','Pam3CSK'};

sti_vec = {'Pam3CSK'};%,'Pam3CSK'
% for proj_num =[18,25:28]
proj_num = 7;
% [para_sample,para_name,estimates] = read_draw_individual_parameters(proj_num,data_save_file_path);
% plot_est_distribution(para_sample(:,1:6), para_name(1:6),estimates)


% [para_sample,para_name,estimates] = read_draw_individual_parameters_seperate_dose(proj_num,data_save_file_path);
data_save_file_path = '../SAEM_proj_2022_2/';

[para_sample,para_name,estimates] = read_draw_individual_parameters_2_seperate_dose(proj_num,data_save_file_path);
stimuli_info_tbl = get_stimuli_info_tbl();
i_sti = 1;
estimates.cell_num = stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == sti_vec{i_sti}));
% i_sti = 2;
% estimates.cell_num = [estimates.cell_num;stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == sti_vec{i_sti}))];
estimates.dose_color = {'b','g','r'};%,'b','g','r'
% estimates.dose_color = {'b','b','b','g','g','g'};

ppt_ParameterEstimation_plot_est_distribution(para_sample(:,1:end-1), para_name(1:end-1),estimates)

% plot_est_distribution_random_effects(para_sample(:,1:6), para_name(1:6),estimates)
% plot_est_distribution(para_sample, para_name)
% fig_num =1;
% para_sample_all{i_sti} = para_sample(:,1:4)
% plot_est_distribution_samefig(para_sample(:,1:4), para_name(1:4),fig_num)

figure(2)
% saveas(gcf, strcat(fig_save_path,'core_parameter_distrib_',sti_vec{i_sti}),'epsc')
saveas(gcf, strcat(fig_save_path,'core_parameter_distrib_',sti_vec{1},'_',sti_vec{2}),'epsc')

close
i_sti = i_sti+1;
% end
