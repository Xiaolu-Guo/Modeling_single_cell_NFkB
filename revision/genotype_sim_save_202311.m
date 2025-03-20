function sim_data_tbl = genotype_sim_save_202311(sim_info,data_info,gene_info) %,codon_info


for i_para_vec = 1:length(gene_info.parameter_name_vec)

    sim_data_tbl = SimDataTblInitialize();
    data_info.flag = '';
    sim_info.parameter_name = gene_info.parameter_name_vec{i_para_vec};
    if isfield(gene_info,'species_name_vec')
        sim_info.species_vec.name = gene_info.species_name_vec{i_para_vec};
    end
    %para_val = gene_info.para_val_vec{i_para_vec};

    for i_gene_type = 1:length(gene_info.gene_type)

        data_info.type = gene_info.gene_type{i_gene_type};
        sim_info.parameter_value_vec = gene_info.gene_parameter_value_vec_genotype{i_gene_type}{i_para_vec};
        sim_info.fold_change_vec = ones(size(sim_info.parameter_value_vec));
        if isfield(gene_info,'species_name_vec')

            sim_info.species_vec.val = gene_info.species_value_vec_genotype{i_gene_type}{i_para_vec};
        end
        %         if isfield(gene_info,'gene_fold_change_vec_genotype')
        %             sim_info.fold_change_vec = gene_info.gene_fold_change_vec_genotype{i_gene_type}{i_para_vec};
        %             sim_info.parameter_value_vec = diag(para_val) * sim_info.fold_change_vec;
        %
        %         elseif isfield(gene_info,'gene_parameter_value_vec_genotype')
        %             sim_info.parameter_value_vec = gene_info.gene_parameter_value_vec_genotype{i_gene_type}{i_para_vec};
        %             sim_info.fold_change_vec  = diag(1./para_val) * sim_info.parameter_value_vec;
        %
        %         else
        %             error('Paramter value/fold change vectors are not specified');
        %         end

        sim_data_tbl_tmpt = NFkB_signaling_para_conc_sim_2023(sim_info,data_info);
        sim_data_tbl = [sim_data_tbl;sim_data_tbl_tmpt];

    end

    % current_fold = pwd;

    %     for i_genotype = 1
    %         gtype = gene_info.gene_type{i_genotype};
    %         traj = ...
    %             sim_data_tbl.trajectory(sim_data_tbl.dose_val == 10 ...
    %             & sim_data_tbl.species =='nucNFkB'...
    %             & sim_data_tbl.type == gtype, :);
    %         data.info_ligand{i_genotype} = 'TNF';
    %         data.model_sim{i_genotype} = traj(:,1:5:end);
    %         data.info_dose_index{i_genotype} = 3;
    %         data.info_dose_str{i_genotype} = '10ng/mL';
    %         data.info_num_cells{i_genotype} = size(data.model_sim{i_genotype},1);
    %         data.order{i_genotype} = (1:data.info_num_cells{i_genotype})';
    %     end
    %
    %     data.exp = data.model_sim;

    % cd(codon_info.codon_path)
    %     vis_data_field = {'model_sim'};%,'sample'};
    %     data_label = {'simulation'};%,'sample'};
    %[collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter

    % cd(current_fold)

    if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
        save_file_name = data_info.save_file_name;
        for i_para_name = 1:length(sim_info.parameter_name)
            if length(save_file_name) < 114

                save_file_name = strcat(save_file_name,'_',sim_info.parameter_name{i_para_name});
            end
        end
        if isfield(sim_info,'species_vec')
            for i_spec_name = 1:length(sim_info.species_vec.name)
                save_file_name = strcat(save_file_name,'_',sim_info.species_vec.name{i_spec_name});
            end

        end
        save_file_name = strcat(save_file_name,'.mat');
        save(strcat(data_info.save_file_path,save_file_name),'sim_data_tbl');%,'collect_feature_vects','metrics','data');
    end

end



end