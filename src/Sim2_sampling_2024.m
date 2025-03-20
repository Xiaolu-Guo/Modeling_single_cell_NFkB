
% Sim2: sampling the original dose and check whether the distribution are consistent. 
% (TNF simulation: first sample and simulate, then delete all the trajectories
% with CV smaller than a threshold.) 

clear all
vers = '202401';

run_sim_example = 0;
draw_fig_traj_heatmap = 0;
draw_fig_para_distr = 0;
run_save_codon = 0;
draw_codon = 0;
draw_spearman_corr = 0;
save_para_codon_corr_python = 0;
% visualize_para_codon_importance_python = 0;

run_sampling = 1;

%_supriya_data_sampling
run_co_sti = 0;
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

if run_sampling
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    %    run_sampling_simulation(monolix_data_save_file_path,fig_save_path)
    %    TNFo_dual_para_sample(vers,data_save_file_path,Num_sample)
    %     Num_sample = 500;
    
    %   para_sample_sim(data_save_file_path,Num_sample)
    %     var_fold = [2,4,8];
    %     ligand_str_ = 'TNF';
    % run all conditions
    %     Initializing_sampling_supriya_data_2023
    input_paras.Num_sample = 10;
    % input_paras.proj_num_vec ={[6,3,5]
    %     [6,3,7]};
    % input_paras.proj_ligand_vec = {{'polyIC','TNF','CpG'}
    %     {'polyIC','TNF','Pam3CSK'}};
    % input_paras.proj_dose_str_vec = { {'100ug','10ng','100nM'}
    %     {'100ug','10ng','100ng'}};
    % input_paras.proj_dose_val_vec = {{100000,10,100}
    %     {100000,10,100}};
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    input_paras.proj_num_vec = {[2,3,4,5,6]};
    input_paras.Num_sample = 200;% [30;10;10;10;10]/5;%50*
    
    input_paras.proj_ligand_vec = {{'TNF','LPS','CpG','PolyIC','Pam3CSK'}};
    input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
    input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
    
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    cal_codon =1;
    save_metric_name = 'Sim2_r3_srs_codon_metric_2024.mat';
    
    
    para_sample_srs_2024(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)
    
    
end

%% visualization
    %     for i_var_fold = 1:length(var_fold)
    %         para_sample_sim_sti_doses_small_var(data_save_file_path,Num_sample,ligand_str_,var_fold(i_var_fold))
    %     end
    %
    %     vers_fig = 'v1';
    %
    %     single_or_dual = 'dual';
    %     draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path)
