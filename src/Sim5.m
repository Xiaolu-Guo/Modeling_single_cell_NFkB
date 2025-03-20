% co-sti simulation, using the original fitted parameters

clear all
vers = '202305';

run_sim_example = 0;
draw_fig_traj_heatmap = 0;
draw_fig_para_distr = 0;
run_save_codon = 0;
draw_codon = 0;
draw_spearman_corr = 0;
save_para_codon_corr_python = 0;
% visualize_para_codon_importance_python = 0;

run_sampling = 0;

%_supriya_data_sampling
run_co_sti = 1;
draw_co_sti_sampling = 0;

%% initializing
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

if run_co_sti
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    % run all conditions
    % Initializing_sampling_supriya_data_2023
    
    Initializing_sampling_supriya_data_2023_05
    % single: 1,2,4,6,7
    % dual: 3,5,9,10,14,15,20,26,27
    % tri 21,22
    
    dual_ligand_length = 10;
    single_ligand_length = 5;


   % input_paras.Num_sample = 1000;
    
    index_sim = [3;5;9;10;14;15;20;26;27;29;1;2;4;6;7];%1;2;6;7;11;
    input_paras.proj_num_vec = input_paras.proj_num_vec(index_sim);
    input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_sim);
    input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_sim);
    input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_sim);
    input_paras.Num_sample = [1*ones(dual_ligand_length,1);
        2*ones(single_ligand_length,1)]*300;
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    save_metric_name = 'Sim5_codon_all5dose_metric.mat';
    cal_codon =1;
    % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras)
    para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    
    vers_fig = '2023';
    single_or_dual = 'dual';
    %     para_sample_sim(data_save_file_path,Num_sample)
    %     draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path)
end
