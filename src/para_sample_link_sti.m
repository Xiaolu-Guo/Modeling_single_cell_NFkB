% TNFo_dual_para_sample_main
% clear all
function [] = para_sample_link_sti(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon)
data_info.save_file_path = data_save_file_path;

if nargin <4
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

if isfield(input_paras,'var_fold_vec')
    if length(input_paras.var_fold_vec) == length(proj_num_vec)
        var_fold_vec = input_paras.var_fold_vec;
    else
        var_fold_vec = input_paras.var_fold_vec * ones(size(proj_num_vec));
    end
else
    var_fold_vec = ones(size(proj_num_vec));
end


gene_info.gene_type = {'wt'};
gene_info.gene_parameter_value_vec_genotype = cell(0);
save_metric_name = strcat(data_save_file_path,'link_cell_metrics');


for i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    var_fold = var_fold_vec(i_proj_num_vec);
    
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            [para_val,est]  = sample_est_parameters_2(proj_num_vec{i_proj_num_vec},Num_sample,monolix_data_save_file_path,var_fold);
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';
        otherwise
            % sample bi-stimulation
            [para_sample_multi_ligand,estimates_2ligand] = sample_est_parameters_multi_receptor(proj_num_vec{i_proj_num_vec},Num_sample,monolix_data_save_file_path,var_fold);
            gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.NFkB_index)';
    end
    
    save_filename_ini = strcat('para_sampled_link_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    for i_ligand = 1:length(proj_ligand_vec{i_proj_num_vec})
        save_filename_ini = strcat(save_filename_ini,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand},...
            '_',replace(replace(proj_dose_str_vec{i_proj_num_vec}{i_ligand},'/',''),'.','p'));
    end
    
    for i_ligand = 1:length(proj_ligand_vec{i_proj_num_vec})
        
        % stimili info
        sim_info.ligand = proj_ligand_vec{i_proj_num_vec}(i_ligand);
        sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec}(i_ligand);
        sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec}(i_ligand);
        save_filename = save_filename_ini;
        for i_ligand_inner = 1:length(sim_info.ligand)
            save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand_inner},...
                '_',replace(replace(sim_info.dose_str{i_ligand_inner},'/',''),'.','p'));
        end
        
        % species that will be saved
        % must be r x 1, for each cell i must be ri x 1
        data_info.species_outputname = {'nucNFkB'};
        data_info.species_composition = {{'NFkBn';'IkBaNFkBn'}};
        data_info.save_file_name = save_filename; % only beginning, no .mat
        
        sim_data_tbl = genotype_sim_save(sim_info,data_info,gene_info);
        
        
        if cal_codon
            
            ligand_str= proj_ligand_vec{i_proj_num_vec}{i_ligand};
            dose_str = proj_dose_str_vec{i_proj_num_vec}{i_ligand};

            data.info_ligand{i_ligand} = ligand_str;
            data.model_sim{i_ligand} = sim_data_tbl.trajectory(:,1:5:end);
            data.info_dose_index{i_ligand} = 1;
            data.info_dose_str{i_ligand} = dose_str;
            data.info_num_cells{i_ligand} = size(sim_data_tbl.trajectory,1);
            [~, data.order{i_ligand}] = sort(max(sim_data_tbl.trajectory,[],2),'descend');
            save_metric_name = strcat(save_metric_name, '_',ligand_str);

        end
    end
end

if cal_codon
    
    data.exp = data.model_sim;
    
    % cd(codon_info.codon_path)
    vis_data_field = {'model_sim'};%,'sample'};
    data_label = {'simulation'};%,'sample'};
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
    
    save(save_metric_name,'data','metrics','collect_feature_vects');
end


