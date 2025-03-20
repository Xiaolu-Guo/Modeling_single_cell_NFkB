function [] = ParameterScan_Sens_traj(data_save_file_path)

% module_info.parameter_sheet = 'param_setting_onlysens';% 'param_setting','param_setting_v0'
% module_info.parameter_sheet = 'param_setting_50sens';% 'param_setting','param_setting_v0'
module_info.parameter_sheet = 'param_setting_FaysParameter';%

module_info.ligand = {'TNF'};
module_info.dose_str = {'100pg/mL','1ng/mL','10ng/mL'};
module_info.dose_val = {0.1,1,10};
Module_parameters_sim(data_save_file_path,module_info)
    
module_info.ligand = {'LPS'};
module_info.dose_str = {'1ng','3ng','10ng','33ng','100ng'};
module_info.dose_val = {1,3,10,33,100};
Module_parameters_sim(data_save_file_path,module_info)

module_info.ligand = {'CpG'};
module_info.dose_str = {'10nM','33nM','100nM','333nM','1uM'};
module_info.dose_val = {10,33,100,333,1000};
Module_parameters_sim(data_save_file_path,module_info)

module_info.ligand = {'polyIC'};
module_info.dose_str = {'10ug','33ug','100ug'};
module_info.dose_val = {1000*10,1000*33,1000*100};
Module_parameters_sim(data_save_file_path,module_info)

module_info.ligand = {'Pam3CSK'};
module_info.dose_str = {'10ng','100ng','1ug'};
module_info.dose_val = {10,100,1000};
Module_parameters_sim(data_save_file_path,module_info)

module_info.ligand = {'TNF'};
module_info.module = {'core'};
module_info.dose_str = {'100pg/mL','1ng/mL','10ng/mL'};
module_info.dose_val = {0.1,1,10};
Module_parameters_sim(data_save_file_path,module_info)