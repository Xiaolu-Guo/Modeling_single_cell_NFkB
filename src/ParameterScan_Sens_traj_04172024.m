function [] = ParameterScan_Sens_traj_04172024()
% For Guo et al. Figure S1, sensitivity analysis of NFkB signaling network
% 
% tested 05/09/2024, Matlab 2020a

data_save_file_path = '../raw_data2023/';
fig_save_path = '../SubFigures2023/';

% module_info.parameter_sheet = 'param_setting_onlysens';% 'param_setting','param_setting_v0'
% module_info.parameter_sheet = 'param_setting_50sens';% 'param_setting','param_setting_v0'
module_info.parameter_sheet = 'param_setting_FaysParameter';%
% ade's doses

%% TNF
module_info.ligand = {'TNF','TNF','TNF'};
module_info.dose_str = {'100pg/mL','1ng/mL','10ng/mL'};
module_info.dose_val = {0.1,1,10};
para_type_vec = {'NFkBtot','syn','deg','TAKac','timed1','timed2'}; % ,'NFkBtot'
for i_para_type = 1:length(para_type_vec)
    module_info.para_type = para_type_vec([i_para_type,i_para_type,i_para_type]);
    if 1 % run to get the data

    Module_parameters_sim_04172024(data_save_file_path,module_info)
    end
    
    if 1 % draw the data
        draw_traj_codon_04182024(data_save_file_path,fig_save_path, module_info)
    
    end
end

%% LPS
module_info.ligand = {'LPS','LPS','LPS'};
module_info.dose_str = {'1ng','3ng','10ng'};
module_info.dose_val = {1,3,10};
para_type_vec = {'NFkBtot','syn','deg','endo','TAKac','timed1','timed2'}; % ,'NFkBtot'

for i_para_type = 1:length(para_type_vec)
    module_info.para_type = para_type_vec([i_para_type,i_para_type,i_para_type]);
    if 1

    Module_parameters_sim_04172024(data_save_file_path,module_info)
    end
    
    if 1
        draw_traj_codon_04182024(data_save_file_path,fig_save_path, module_info)
    
    end
end


%% CpG
module_info.ligand = {'CpG','CpG','CpG'};
module_info.dose_str = {'33nM','100nM','333nM'};
module_info.dose_val = {33,100,333};
para_type_vec = {'NFkBtot','syn','deg','endo','TAKac','timed1','timed2'}; % ,'NFkBtot'
for i_para_type = 1:length(para_type_vec)
    module_info.para_type = para_type_vec([i_para_type,i_para_type,i_para_type]);
    if 1

    Module_parameters_sim_04172024(data_save_file_path,module_info)
    end
    
    if 1
        draw_traj_codon_04182024(data_save_file_path,fig_save_path, module_info)
    
    end
end

%% polyIC
module_info.ligand = {'polyIC','polyIC','polyIC'};
module_info.dose_str = {'10ug','33ug','100ug'};
module_info.dose_val = {1000*10,1000*33,1000*100};
para_type_vec = {'NFkBtot','syn','deg','endo','TAKac','timed1','timed2'}; % ,'NFkBtot'
for i_para_type = 1:length(para_type_vec)
    module_info.para_type = para_type_vec([i_para_type,i_para_type,i_para_type]);
    if 1

    Module_parameters_sim_04172024(data_save_file_path,module_info)
    end
    
    if 1
        draw_traj_codon_04182024(data_save_file_path,fig_save_path, module_info)
    
    end
end


%% Pam3CSK

module_info.ligand = {'Pam3CSK','Pam3CSK','Pam3CSK'};
module_info.dose_str = {'10ng','100ng','1ug'};
module_info.dose_val = {10,100,1000};
para_type_vec = {'NFkBtot','syn','deg','TAKac','timed1','timed2'}; % ,'NFkBtot'
for i_para_type = 1:length(para_type_vec)
    module_info.para_type = para_type_vec([i_para_type,i_para_type,i_para_type]);
    if 1

    Module_parameters_sim_04172024(data_save_file_path,module_info)
    end
    
    if 1
        draw_traj_codon_04182024(data_save_file_path,fig_save_path, module_info)
    
    end
end



% module_info.ligand = {'TNF'};
% module_info.module = {'core'};
% module_info.dose_str = {'100pg/mL','1ng/mL','10ng/mL'};
% module_info.dose_val = {0.1,1,10};
% Module_parameters_sim_04172024(data_save_file_path,module_info)