function [] = draw_figure6_cal_codon(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)

if nargin <4
    cal_codon = 0;
end

proj_num_vec = input_paras.proj_num_vec;
proj_ligand_vec = input_paras.proj_ligand_vec;
proj_dose_str_vec = input_paras.proj_dose_str_vec;
proj_dose_val_vec = input_paras.proj_dose_val_vec;

if isfield(input_paras,'Num_sample')
    Num_sample = input_paras.Num_sample;
else
    Num_sample = 500;
end
Num_sample = 2;

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


for i_proj_num_vec = 1:length(proj_num_vec)
    
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
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    save_filename = strcat('para_sampled_NFkBinit_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    
    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end
    
    for i_para_vec = 1:length(gene_info.parameter_name_vec)
        
        sim_info.parameter_name = gene_info.parameter_name_vec{i_para_vec};
        
        %para_val = gene_info.para_val_vec{i_para_vec};
        
        
        for i_para_name = 1:length(sim_info.parameter_name)
            save_filename = strcat(save_filename,'_',sim_info.parameter_name{i_para_name});
        end
        save_filename= strcat(save_filename,'.mat');
    end
    
    load(strcat(data_save_file_path,save_filename))
    
    
    
    if cal_codon
        
        i_ligand = 1;
        ligand_str= proj_ligand_vec{i_proj_num_vec}{i_ligand};
        dose_str = proj_dose_str_vec{i_proj_num_vec}{i_ligand};
        
        for i_ligand = 2:length(proj_ligand_vec{i_proj_num_vec})
            ligand_str = strcat(ligand_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
            dose_str = strcat(dose_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
        end
        
        data.info_ligand{i_proj_num_vec} = ligand_str;
        data.model_sim{i_proj_num_vec} = sim_data_tbl.trajectory(:,1:5:end);
        data.info_dose_index{i_proj_num_vec} = 1;
        data.info_dose_str{i_proj_num_vec} = dose_str;
        data.info_num_cells{i_proj_num_vec} = size(sim_data_tbl.trajectory,1);
        [~, data.order{i_proj_num_vec}] = sort(max(sim_data_tbl.trajectory,[],2),'descend');
        
    end
    
    
   
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
    save(save_metric_name,'data','metrics','collect_feature_vects');
end