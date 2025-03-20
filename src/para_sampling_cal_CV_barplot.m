% TNFo_dual_para_sample_main
% clear all
function [rcp_cvs, adp_cvs,core_cvs] = para_sampling_cal_CV_barplot(var_input,data_save_file_path,input_paras,cal_codon,save_metric_name)


% input: paranames : parameter names that should be specified
% paravals : parameter value that should be specified: must be a p by 1

paranames = var_input.paranames;
paravals = var_input.paravals;
specienames = var_input.specienames;
specievals = var_input.specievals;


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

sim_data_tbl = SimDataTblInitialize();
rcp_cvs =[];
core_cvs = [];
adp_cvs = [];
for i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    
    var_fold = var_fold_mat(i_proj_num_vec);
    % sample multi-stimulation
    i_data = find(strcmp(data_fitting_all.info_ligand,sim_info.ligand));
    para_val_ele = [];
    for i_i_data= 1:length(i_data)
        para_val_ele = [para_val_ele;data_fitting_all.parameters_mode_nan{i_data(i_i_data)} ] ; % [cellnumber x parameter]
    end
    i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1);
    
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
    
    %est.name
    cv_vals = std(para_val)./mean(para_val);
    switch sim_info.ligand{1}
        case 'TNF'
            rcp_cvs = [rcp_cvs;cv_vals(4);cv_vals(5)];
            core_cvs = [core_cvs;cv_vals(2);cv_vals(3);cv_vals(6)];
            adp_cvs = [adp_cvs;cv_vals(1)];
        case 'Pam3CSK'
            rcp_cvs = [rcp_cvs;cv_vals(4);cv_vals(5)];
            core_cvs = [core_cvs;cv_vals(2);cv_vals(3);cv_vals(6)];
            adp_cvs = [adp_cvs;cv_vals(1)];
        otherwise
            rcp_cvs = [rcp_cvs;cv_vals(4);cv_vals(5);cv_vals(6)];
            core_cvs = [core_cvs;cv_vals(2);cv_vals(3);cv_vals(7)];
            adp_cvs = [adp_cvs;cv_vals(1)];
    end
    
end


