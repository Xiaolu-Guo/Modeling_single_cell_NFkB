% TNFo_dual_para_sample_main
% clear all
function [] = para_sample_all_scramble_fitting_sim_sti_doses_var_2024_02(data_save_file_path,input_paras,cal_codon,save_metric_name,fitting_dose)

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

ncpu=10;
pc=parcluster('local');
pc.NumThreads=2;%
parpool(pc,ncpu)

parfor i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    gene_info = struct;
    sim_info = struct;
    data_info = struct;
    est = struct;
    
    gene_info.gene_type = {'wt'};
    gene_info.gene_parameter_value_vec_genotype = cell(0);
    data_info.save_file_path = data_save_file_path;
    
    
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    
    var_fold = var_fold_mat(i_proj_num_vec);
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            %             [para_val,est]  = sample_est_parameters_2023(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),monolix_data_save_file_path,var_fold);
            %             NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            %             non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            %             gene_info.parameter_name_vec = {est.name(non_NFkB_index)};
            %             gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,non_NFkB_index)';
            %             gene_info.species_name_vec= {{'NFkB'}};
            %             gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';
            switch fitting_dose
                case 'middose'
                    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1,'last')-1;
                    para_val_ele = data_fitting.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
                case 'lowdose'
                    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1);
                    para_val_ele = data_fitting.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
                    
                case 'highdose'
                    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1,'last');
                    para_val_ele = data_fitting.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
                    
                case 'alldose'
                    i_data = find(strcmp(data_fitting_all.info_ligand,sim_info.ligand));
                    para_val_ele = [];
                    for i_i_data= 1:length(i_data)
                        para_val_ele = [para_val_ele;data_fitting_all.parameters_mode_nan{i_data(i_i_data)} ] ; % [cellnumber x parameter]
                    end
                    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1);
                    
                case 'originaldose'
                    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand) &...
                        strcmp(data_dose_str, sim_info.dose_str)) ;
                    if isempty(strcmp(data_dose_str, sim_info.dose_str))
                        error('dose_str is not valid for parameter estimation')
                    end
                    para_val_ele = data_fitting.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
                otherwise
                    error('wrong input for fitting_dose value!')
            end
            
            rpt_time = ceil(Num_sample(i_proj_num_vec)*5/size(para_val_ele,1));
            para_val_total = [];
            for i_rpt = 1:rpt_time
                para_val_total = [para_val_total;para_val_ele];
            end
            para_val = para_val_total(randperm(size(para_val_total,1),Num_sample(i_proj_num_vec)),:);
            est.name = cellfun(@replace,data_fitting.para_name{i_data},old_str(1:length(data_fitting.para_name{i_data})),new_str(1:length(data_fitting.para_name{i_data})),'UniformOutput',false);
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            shift_index = find(strcmp(est.name,'shift'));
            parameter_index = setdiff(setdiff(1:length(est.name),NFkB_index,'stable'),shift_index,'stable');
            
            % non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(parameter_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,parameter_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';
            
        otherwise
            % sample multi-stimulation
            
            [para_sample_multi_ligand,estimates_2ligand] = sample_all_scramble_fitting_multi_receptor_2024_02(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),data_fitting,data_fitting_all,sim_info,'alldose');
            
            gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.NFkB_index)';
    end
    
    % stimili info
    
    save_filename = strcat('202312_para_sampled_scramble_fitting_NFkBinit_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    
    switch length(sim_info.ligand)
        
        case 1
            save_filename = strcat(save_filename,'_',sim_info.ligand{1},...
                '_',replace(replace(sim_info.dose_str{1},'/',''),'.','p'));
            
        case 2
            save_filename = strcat(save_filename,'_',sim_info.ligand{1},...
                '_',replace(replace(sim_info.dose_str{1},'/',''),'.','p'),...
                '_',sim_info.ligand{2},...
                '_',replace(replace(sim_info.dose_str{2},'/',''),'.','p'));
            
        case 3
            save_filename = strcat(save_filename,'_',sim_info.ligand{1},...
                '_',replace(replace(sim_info.dose_str{1},'/',''),'.','p'),...
                '_',sim_info.ligand{2},...
                '_',replace(replace(sim_info.dose_str{2},'/',''),'.','p'),...
                '_',sim_info.ligand{3},...
                '_',replace(replace(sim_info.dose_str{3},'/',''),'.','p'));
        case 4
            save_filename = strcat(save_filename,'_',sim_info.ligand{1},...
                '_',replace(replace(sim_info.dose_str{1},'/',''),'.','p'),...
                '_',sim_info.ligand{2},...
                '_',replace(replace(sim_info.dose_str{2},'/',''),'.','p'),...
                '_',sim_info.ligand{3},...
                '_',replace(replace(sim_info.dose_str{3},'/',''),'.','p'),...
                '_',sim_info.ligand{4},...
                '_',replace(replace(sim_info.dose_str{4},'/',''),'.','p'));
        case 5
            save_filename = strcat(save_filename,'_',sim_info.ligand{1},...
                '_',replace(replace(sim_info.dose_str{1},'/',''),'.','p'),...
                '_',sim_info.ligand{2},...
                '_',replace(replace(sim_info.dose_str{2},'/',''),'.','p'),...
                '_',sim_info.ligand{3},...
                '_',replace(replace(sim_info.dose_str{3},'/',''),'.','p'),...
                '_',sim_info.ligand{4},...
                '_',replace(replace(sim_info.dose_str{4},'/',''),'.','p'),...
                '_',sim_info.ligand{5},...
                '_',replace(replace(sim_info.dose_str{5},'/',''),'.','p'));
            
        otherwise
            error('more than 5 ligand stimulation is not supported')
    end
    
    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB';'TNFR';'TLR4';'TLR2';'TLR3';'TLR9';'IKK';'TAK1';'IkBamRNA'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'TLR4'};{'TLR2'};{'TLR3'};{'TLR9'};{'IKK'};{'TAK1'};{'IkBat'}};
    data_info.save_file_name = save_filename; % only beginning, no .mat
    
    sim_data_tbl = genotype_sim_save_2023(sim_info,data_info,gene_info);
    
    
    if cal_codon
        
        
        switch length(proj_ligand_vec{i_proj_num_vec})
            
            case 1
                ligand_str= proj_ligand_vec{i_proj_num_vec}{1};
                dose_str = proj_dose_str_vec{i_proj_num_vec}{1};
                
            case 2
                ligand_str = strcat(proj_ligand_vec{i_proj_num_vec}{1},...
                    '_',proj_ligand_vec{i_proj_num_vec}{2});
                dose_str = strcat(proj_dose_str_vec{i_proj_num_vec}{1},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{2});
                
            case 3
                ligand_str = strcat(proj_ligand_vec{i_proj_num_vec}{1},...
                    '_',proj_ligand_vec{i_proj_num_vec}{2},...
                    '_',proj_ligand_vec{i_proj_num_vec}{3});
                dose_str = strcat(proj_dose_str_vec{i_proj_num_vec}{1},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{2},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{3});
                
            case 4
                ligand_str = strcat(proj_ligand_vec{i_proj_num_vec}{1},...
                    '_',proj_ligand_vec{i_proj_num_vec}{2},...
                    '_',proj_ligand_vec{i_proj_num_vec}{3},...
                    '_',proj_ligand_vec{i_proj_num_vec}{4});
                dose_str = strcat(proj_dose_str_vec{i_proj_num_vec}{1},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{2},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{3},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{4});
            case 5
                ligand_str = strcat(proj_ligand_vec{i_proj_num_vec}{1},...
                    '_',proj_ligand_vec{i_proj_num_vec}{2},...
                    '_',proj_ligand_vec{i_proj_num_vec}{3},...
                    '_',proj_ligand_vec{i_proj_num_vec}{4},...
                    '_',proj_ligand_vec{i_proj_num_vec}{5});
                dose_str = strcat(proj_dose_str_vec{i_proj_num_vec}{1},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{2},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{3},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{4},...
                    '_',proj_dose_str_vec{i_proj_num_vec}{5});
            otherwise
                error('more than 5 ligand stimulation is not supported')
        end
        
        
        
        info_ligand{i_proj_num_vec} = ligand_str;
        model_sim{i_proj_num_vec} = sim_data_tbl.trajectory(:,1:5:end);
        info_dose_index{i_proj_num_vec} = 1;
        info_dose_str{i_proj_num_vec} = dose_str;
        info_num_cells{i_proj_num_vec} = size(sim_data_tbl.trajectory,1);
        [~, order{i_proj_num_vec}] = sort(max(sim_data_tbl.trajectory,[],2),'descend');
        
    end
    % save(save_metric_name,'data');
end

delete(gcp)

data.info_ligand = info_ligand;
data.model_sim = model_sim;
data.info_dose_index = info_dose_index;
data.info_dose_str = info_dose_str;
data.info_num_cells = info_num_cells;
data.order = order;

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


