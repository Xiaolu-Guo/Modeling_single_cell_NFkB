% TNFo_dual_para_sample_main
% clear all
function [] = para_sample_fitting_match_sim_SRS_IkBao_202401(para6fold,data_save_file_path,input_paras,cal_codon,save_metric_name,fitting_dose)
data_info.save_file_path = data_save_file_path;

if nargin < 4
    cal_codon = 0;
end

if nargin < 6
    fitting_dose = 'alldose'; % middose, lowdose, alldose, highdose, originaldose
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

%% parameter from the data

if 1
    data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
    
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    
    data_fitting = data;
    clear data
    for i_data = 1:length(data_fitting.pred_mode)
        
        data_fitting.pred_mode_cv{i_data} = std(data_fitting.pred_mode_filter_nan{i_data},[],2)./mean(data_fitting.pred_mode_filter_nan{i_data},2);
        [~,data_fitting.pred_mode_cv_order{i_data}] = sort(data_fitting.pred_mode_cv{i_data},'descend');
        
    end
    
    % input
    thresh_TNF=0.33;
    thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
    data_fitting.parameters_mode_filter_TNF = cell(2,1);
    for i_data = 1:3
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(index ,:);
    end
    
    
    for i_data = 4:length(data_fitting.parameters_mode_nan)
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(: ,:);
    end
    % change the parameters for TNF only for high CV
    data_fitting.parameters_mode_nan = data_fitting.parameters_mode_filter_TNF;
    
    data_field_names = fieldnames(data_fitting);
    
    data_index = [1:6,10:12,14:19];
    
    for i_data_field = 1:length(data_field_names)
        data_new.(data_field_names{i_data_field}) = data_fitting.(data_field_names{i_data_field}) (data_index);
    end
    
    data_fitting_all = data_fitting;
    data_fitting = data_new;
    
end


%%

old_str = {'p','p','p','p','p','p','p','p','p','p'};
new_str = {'params','params','params','params','params','params','params','params','params','params'};


for i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    
    var_fold = var_fold_mat(i_proj_num_vec);
    if length(proj_num_vec{i_proj_num_vec})>1
        % sample multi-stimulation
        
        [para_sample_multi_ligand,estimates_2ligand] = sample_fitting_multi_receptor_2023_11(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),data_fitting,data_fitting_all,sim_info,'alldose');
        
        gene_info.parameter_name_vec = {{estimates_2ligand.name{estimates_2ligand.non_NFkB_index},'params6'}};
        gene_info.gene_parameter_value_vec_genotype{1}{1} = [para_sample_multi_ligand(:,estimates_2ligand.non_NFkB_index)';
            ones(1,size(para_sample_multi_ligand,1))*6e-05*para6fold];
        gene_info.species_name_vec= {{'NFkB'}};
        gene_info.species_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.NFkB_index)';
        
    else           
        error('at least dual ligand for SRS analysis!!')
    end
    
    % stimili info
    
    save_filename = strcat('202305_para_sampled_fitting_NFkBinit_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    
    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end
    
    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB';'TNFR';'TLR4';'TLR2';'TLR3';'TLR9';'IKK';'TAK1';'IkBamRNA'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'TLR4'};{'TLR2'};{'TLR3'};{'TLR9'};{'IKK'};{'TAK1'};{'IkBat'}};
    data_info.save_file_name = save_filename; % only beginning, no .mat
    
    % sim_info_single_ligand = sim_info;
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
        sim_data_tbl_tmpt = genotype_sim_save_2023(sim_info_single_ligand,data_info,gene_info);%sim_info_single_ligand
        
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
        %clear sim_info_single_ligand
        
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
    save(strcat(data_save_file_path,save_metric_name),'sim_data_tbl','data','metrics','collect_feature_vects');
    
end


