function [] = Module_parameters_sim(data_save_file_path,module_info)

% Module Senstive Analysis
% 
% addpath('../NFkB_common/')
% 
% addpath('../NFkB_sim/')
data_info.save_file_path = data_save_file_path;%  '../../NFkB_para_estm_project/NFkB_figures/ParameterScan/data/';


sim_data_tbl = SimDataTblInitialize();

% sim_info.para_setting_sheet = 'param_setting'; %default

sim_info.fold_change_vec = [10.^linspace(-1,1,21); 10.^linspace(-1,1,21)];


if nargin<2
    sim_info.ligand = {'TNF'};
    sim_info.dose_str = {'100pg/mL','1ng/mL','10ng/mL'};
    sim_info.dose_val = {0.1,1,10};
    module_info.parameter_sheet = 'param_setting';
    module_info.ligand = {'TNF'};
else
    sim_info.ligand = module_info.ligand;
    sim_info.dose_str = module_info.dose_str;
    sim_info.dose_val = module_info.dose_val;
end

data_info.species_outputname = {'nucNFkB';'TNFR';'IKK'};
data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'IKK'}}; % must
% be r x 1, for each cell i must be ri x 1

data_info.flag = '';
data_info.type = 'wt';

%% read parameters:

opts1 = detectImportOptions('parameter_setting.xlsx','Sheet',module_info.parameter_sheet);%sim_info.para_setting_sheet,'param_setting'
for jj = 1: length(opts1.VariableNames)
    opts1 = setvartype(opts1, opts1.VariableNames{jj}, 'char');
end
parameter_setting = readtable('parameter_setting.xlsx',opts1);

if isfield(module_info,'module')
    index_distribute_para = find(strcmp(module_info.module,parameter_setting.Module));
    data_info.save_file_name = strcat('Module_Sens_',module_info.module{1});

else
    index_distribute_para = find(strcmp(module_info.ligand,parameter_setting.Module));
    data_info.save_file_name = strcat('Module_Sens_',module_info.ligand{1});
end



for i_para=1:length(index_distribute_para)
    i_para
    %  parameter_setting.parameter
    % Specific Paramters main function
    sim_info.parameter_name = parameter_setting.parameter(index_distribute_para(i_para));
    
    para_val = str2num(parameter_setting.Value{index_distribute_para(i_para)});%[8.224e-06; 1889];
    sim_info.parameter_value_vec = diag(para_val) * sim_info.fold_change_vec;
    
    sim_data_tbl_tmpt = NFkB_signaling_para_value_sim(sim_info,data_info);
    
    % the parameter numbers has to be the same
    sim_data_tbl= [sim_data_tbl;sim_data_tbl_tmpt];
    
end
%
%     codon calculation
%     current_fold = pwd;
%
%     for i_genotype = 1:2
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
%
%     cd(codon_info.codon_path)
%     vis_data_field = {'model_sim'};%,'sample'};
%     data_label = {'simulation'};%,'sample'};
%     [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
%
%     cd(current_fold)

if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
    save(strcat(data_info.save_file_path,data_info.save_file_name),'sim_data_tbl');
end


% current_folder = pwd;
% cd(current_folder)

