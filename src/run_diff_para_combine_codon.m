function [] = run_diff_para_combine_codon(monolix_data_save_file_path,data_save_file_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3
data_proj_nums ={[18,22,21,20,19,16,24,23];%TNF
    [25,30,34,35,36,29,31,33,32];%LPS
    [26,41,42,43,44,38,40,39];% CpG
    [28,46,47,48,49,45,51,50];
    [27,55,56,57,58,59,52,53,54]};% P3CSK

ligand_all = {'TNF';'LPS';'CpG';'Pam3CSK';'polyIC'};%polyIC

for i_ligand = 1:length(data_proj_nums)
    data_proj_num_vec = data_proj_nums{i_ligand};    
    data = read_monolix_data_to_matrix(data_proj_num_vec,monolix_data_save_file_path);
    vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'}; %,'sample'};
    data_label = {'experiments','prediction'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label); %,  parameter
    save(strcat(data_save_file_path,ligand_all{i_ligand},'_diff_para_combine_codon.mat'),'data','collect_feature_vects','metrics');
    
end



