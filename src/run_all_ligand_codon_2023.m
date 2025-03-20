function [] = run_all_ligand_codon_2023(monolix_data_save_file_path,data_save_file_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

    % data_proj_num_vec = [18,25:28]; % ../SAEM_proj_2022/
    data_proj_num_vec = [2:6]; 
    % data = read_monolix_data_to_matrix(data_proj_num_vec,monolix_data_save_file_path);
   data = read_monolix_data_to_matrix_2023(data_proj_num_vec,monolix_data_save_file_path);

    % load(strcat(data_save_file_path,'All_ligand_codon_2.mat'))% data,collect_feature_vects,metrics));

    vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'}; %,'sample'};
    data_label = {'experiments','prediction'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon_2023(data,vis_data_field,data_label); %,  parameter
    collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
    collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
    collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
    save(strcat(data_save_file_path,'All_ligand_codon_2023.mat'),'data','collect_feature_vects','metrics');




