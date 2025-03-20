function [] = ParameterScan_RMSD(data_save_file_path)

module_vec = {'TNF';'CpG';'Pam3CSK';'LPS';'polyIC';'coreup';'core'};%'TNF';


'_allsti_small_change'; % vers = '_allsti_large_change_brookmodel';
fold_change_vec = [0.1,1,10];
para_setting_sheet = 'param_setting_v0';% 'param_setting';

for i_mod = 1:length(module_vec)
    NFkB_para_analysis_Module(module_vec{i_mod},vers,fold_change_vec,para_setting_sheet,data_save_file_path)
end


vers = '_allsti_small_change';%'_allsti_small_change';
fold_change_vec = [0.9,1,1.1];
para_setting_sheet = 'param_setting';% 'param_setting';

for i_mod = 1:length(module_vec)
    NFkB_para_analysis_Module(module_vec{i_mod},vers,fold_change_vec,para_setting_sheet,data_save_file_path)
end

end