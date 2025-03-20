function [] = run_all_ligand_codon_2023_03(monolix_data_save_file_path,data_save_file_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

    % data_proj_num_vec = [18,25:28]; % ../SAEM_proj_2022/
    data_proj_num_vec = [2:6]; 
    % data = read_monolix_data_to_matrix(data_proj_num_vec,monolix_data_save_file_path);
data_save_file_path = '../raw_data2023/';%_fay_parameter/';

%     data = read_monolix_data_to_matrix_2023(data_proj_num_vec,monolix_data_save_file_path);
 load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));

for i_sti = 1:length(data.pred_mode_filter_nan)
    
    data.pred_mode_cv{i_sti} = std(data.pred_mode_filter_nan{i_sti},[],2)./mean(data.pred_mode_filter_nan{i_sti},2);
    [~,data.pred_mode_cv_order{i_sti}] = sort(data.pred_mode_cv{i_sti},'descend');
    [~,data.pred_mode_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2},'descend');
    [~,data.exp_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2-1},'descend');

end

% input
    thresh_TNF=0.33;
    thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
for i_data = 1:3
    index = data.(thresh_field{1}){i_data} (1: ceil(length(data.(thresh_field{1}){i_data}) * thresh_TNF));
    data.pred_mode_filter_nan{i_data} = data.pred_mode_filter_nan{i_data}(index ,:);
    data.exp_mode_filter_nan{i_data} = data.exp_mode_filter_nan{i_data}(index,:);
    data.parameters_mode_nan{i_data} = data.parameters_mode_nan  {i_data}(index,:);

    [~, data.order{i_data}]=sort(max(data.exp{i_data},[],2),'descend');
end

  vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'}; %,'sample'};
  data_label = {'experiments','prediction'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon_2023(data,vis_data_field,data_label); %,  parameter
    collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
    collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
    collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
%  load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'),'collect_feature_vects','metrics')% data,collect_feature_vects,metrics));

    save(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'),'data','collect_feature_vects','metrics');




