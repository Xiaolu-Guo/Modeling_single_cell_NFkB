% TNFo_dual_para_sample_main
% clear all
function [] = para_sample_srs_2024(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)
data_info.save_file_path = data_save_file_path;

if nargin < 4
    cal_codon = 0;
end
% the paramters that will be sampled

proj_num_vec = input_paras.proj_num_vec;
proj_ligand_vec = input_paras.proj_ligand_vec;
proj_dose_str_vec = input_paras.proj_dose_str_vec;
proj_dose_val_vec = input_paras.proj_dose_val_vec;

if isfield(input_paras,'Num_sample')
    Num_sample = input_paras.Num_sample;
else
    Num_sample = 500;
end

if length(Num_sample) == length(proj_num_vec)
    Num_sample = Num_sample;
else
    Num_sample = Num_sample * ones(size(proj_num_vec));
end

if isfield(input_paras,'var_fold_mat')
    if length(input_paras.var_fold_mat) == length(proj_num_vec)
        var_fold_mat = input_paras.var_fold_mat;
    else
        var_fold_mat = input_paras.var_fold_vec * ones(size(proj_num_vec));
    end
else
    var_fold_mat = ones(size(proj_num_vec));
end



gene_info.gene_type = {'wt'};
gene_info.gene_parameter_value_vec_genotype = cell(0);

for i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    var_fold = var_fold_mat(i_proj_num_vec);
    
    
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            [para_val,est]  = sample_est_parameters_2023(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),monolix_data_save_file_path,var_fold);
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';
        otherwise
            % sample bi-stimulation
            [para_sample_multi_ligand,estimates_2ligand] = sample_est_parameters_multi_receptor_2023(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),monolix_data_save_file_path,var_fold);
            gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.NFkB_index)';
    end
    
    % stimili info
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    save_filename = strcat('202304_para_sampled_NFkBinit_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    
    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end
    
    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB';'TNFR';'TLR4';'TLR2';'TLR3';'TLR9';'IKK';'TAK1';'IkBamRNA'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'TLR4'};{'TLR2'};{'TLR3'};{'TLR9'};{'IKK'};{'TAK1'};{'IkBat'}};
    data_info.save_file_name = save_filename; % only beginning, no .mat
    
    fieldnames_sim_info = fieldnames(sim_info);
    sim_data_tbl = SimDataTblInitialize();
    
    for i_ligand_stim = 1:length(sim_info.ligand)
        for i_fieldnames_sim_info = 1:length(fieldnames_sim_info)
            sim_info_single_ligand.(fieldnames_sim_info{i_fieldnames_sim_info}) = ...
                sim_info.(fieldnames_sim_info{i_fieldnames_sim_info})(i_ligand_stim);
        end
        if isfield(data_info,'save_file_name')
            data_info = rmfield(data_info,'save_file_name');
        end
        
        sim_data_tbl_tmpt = genotype_sim_save_2023(sim_info_single_ligand,data_info,gene_info);
          
        if cal_codon
            
            ligand_str= sim_info.ligand{i_ligand_stim};
            dose_str = sim_info.dose_str{i_ligand_stim};
            
            
            data.info_ligand{i_ligand_stim} = ligand_str;
            data.model_sim{i_ligand_stim} = sim_data_tbl_tmpt.trajectory(:,1:5:end);
            data.info_dose_index{i_ligand_stim} = 1;
            data.info_dose_str{i_ligand_stim} = dose_str;
            data.info_num_cells{i_ligand_stim} = size(sim_data_tbl_tmpt.trajectory,1);
            [~, data.order{i_ligand_stim}] = sort(max(sim_data_tbl_tmpt.trajectory,[],2),'descend');
            
        end
        
        sim_data_tbl = [sim_data_tbl;sim_data_tbl_tmpt];
    end
    % save(save_metric_name,'data');
end

if cal_codon
    
    data.exp = data.model_sim;
    
    % cd(codon_info.codon_path)
    vis_data_field = {'model_sim'};%,'sample'};
    data_label = {'simulation'};%,'sample'};
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
    
    % for i_metrics = length(metrics)
    % pk_time(i_metrics) = mean( metrics{i_metrics}.pk1_time );
    % osc_ratio(i_metrics) = sum(metrics{i_metrics}.oscpower>2e-5)/length(metrics{i_metrics}.oscpower);
    % end
    save(strcat(data_save_file_path,save_metric_name),'data','metrics','collect_feature_vects');
    
end


