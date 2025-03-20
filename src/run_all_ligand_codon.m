function [] = run_all_ligand_codon(monolix_data_save_file_path,data_save_file_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

    % data_proj_num_vec = [18,25:28]; % ../SAEM_proj_2022/
    data_proj_num_vec = [3:7]; 
    % data = read_monolix_data_to_matrix(data_proj_num_vec,monolix_data_save_file_path);
    data = read_monolix_data_to_matrix_2(data_proj_num_vec,monolix_data_save_file_path);

    vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'}; %,'sample'};
    data_label = {'experiments','prediction'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label); %,  parameter
    save(strcat(data_save_file_path,'All_ligand_codon_2.mat'),'data','collect_feature_vects','metrics');




