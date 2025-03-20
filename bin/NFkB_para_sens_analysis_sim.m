function [] = NFkB_para_sens_analysis_sim(sim_info, data_info)

    default_ratio = sim_info.default_ratio ;


% if isfield(sim_info,'Module')
%     switch sim_info.Module
%         case 'TNF'
%             sti_vec = {sim_info.Module};
%         case 'LPS'
%             sti_vec = {sim_info.Module};
%         case 'Pam3CSK'
%             sti_vec = {sim_info.Module};
%         case 'polyIC'
%             sti_vec = {sim_info.Module};
%         case 'CpG'
%             sti_vec = {sim_info.Module};
%         case 'core'
%             sti_vec = {'TNF','LPS','Pam3CSK','polyIC','CpG'};
%         case 'coreup'
%             sti_vec = {'TNF','LPS','Pam3CSK','polyIC','CpG'};
%     end
% end
% 
% if isfield(sim_info,'stimuli')
%     sti_vec = {sim_info.stimuli};
% end

if isfield(data_info,'species')
    % names = {'TNF','TNFR','TNFR_TNF','TTR','C1_off','C1','IkBa','IkBan','IKKIkBa','IKKIkBaNFkB','IkBaNFkB','IkBaNFkBn','IKK','NFkBn'};
    names = data_info.species;
else
    names = {'IKK','IkBaNFkBn','NFkBn'};
end




%% read parameters:

opts1 = detectImportOptions('parameter_setting.xlsx','Sheet',sim_info.para_setting_sheet);
for jj = 1: length(opts1.VariableNames)
    opts1 = setvartype(opts1, opts1.VariableNames{jj}, 'char');
end
parameter_setting = readtable('parameter_setting.xlsx',opts1);


%% set the index
% if the exact parameter specified, then calculate the specified paramters
% otherwishe if module is speicified, then calculate all the parameters in
% the module

if isfield(sim_info,'parameter')
    index_distribute_para = zeros(length(sim_info.parameter),1);
    for i_parameter = 1:length(sim_info.parameter)
        index_distribute_para(i_parameter) = find(strcmp(sim_info.parameter{i_parameter},parameter_setting.parameter));
    end
else
    index_distribute_para = find(strcmp(sim_info.Module,parameter_setting.Module));
    
end

%% initialize the sim_data
sim_data = struct();
sim_data.parameter_module = {};
sim_data.parameter_reac_num = {};
sim_data.parameter_para_num = {};
sim_data.parameter_name(:) = {};
sim_data.parameter_fold_change = {};
sim_data.parameter_value(:) = {};
sim_data.ligand = {};
sim_data.dose_str(:) = {};
sim_data.dose_val = {};
sim_data.species = {};
sim_data.flag = {};
sim_data.type = {};
sim_data.trajectory = {};
run('exp_info_initialization.m')

%% Calculate sim_data
for i_sti = 1:length(sim_info.stimuli)
    for i_para=1:length(index_distribute_para)
        i_sti
        i_para
        %% sim_info value setting
        sim_info_para.reac_num = str2double(parameter_setting.reaction_number(index_distribute_para(i_para)));
        sim_info_para.para_num = str2double(parameter_setting.parameter_number(index_distribute_para(i_para)));
        para_value = str2double(parameter_setting.Value(index_distribute_para(i_para)));
        sim_info_para.parameter_module =parameter_setting.Module{index_distribute_para(i_para)};
        
        sim_info_para.reac_num
        sim_info_para.para_num
        
        para_min_val = parameter_setting.min_val{index_distribute_para(i_para)};
        para_max_val = parameter_setting.max_val{index_distribute_para(i_para)};
        
        if isempty(para_min_val)
            para_vec_min = para_value * default_ratio(1:ceil(length(default_ratio)/2));
        else
            para_vec_min = linspace(str2double(para_min_val),para_value,ceil(length(default_ratio)/2));
        end
        
        if isempty(para_max_val)
            para_vec_max = para_value * default_ratio(ceil(length(default_ratio)/2):end);
        else
            para_vec_max = linspace(para_value,str2double(para_max_val),ceil(length(default_ratio)/2));
        end
            
        if ~isfield(sim_info,'para_fold')
            
            para_vec = [para_vec_min(1:end),para_vec_max(2:end)];
            fold_change_vec = para_vec/para_vec_min(end);
            
        elseif length(sim_info.para_fold)>1
            fold_change_vec = sim_info.para_fold{i_para};
            para_vec =  para_value*fold_change_vec;
        else
            fold_change_vec = sim_info.para_fold{1};
            para_vec =  para_value*fold_change_vec;
        end
        
        if (max(para_vec)>max(para_vec_max)) || (min(para_vec)<min(para_vec_min))
            para_vec = [para_vec_min(1:end),para_vec_max(2:end)];
            fold_change_vec = para_vec/para_vec_min(end);
        end
        para_vec
        sim_info_para.para_vec = para_vec;
        sim_info_para.fold_change_vec = fold_change_vec;
        sim_info_para.doses = cell2mat(dose_val_all{strcmp(ligand_all,sim_info.stimuli(i_sti))});
        sim_info_para.dose_scale = 1/dose_scale_all{strcmp(ligand_all,sim_info.stimuli(i_sti))}; % Convert to uM (from nM)
        sim_info_para.dose_field = dose_all{strcmp(ligand_all,sim_info.stimuli(i_sti))};
        sim_info_para.stimuli = sim_info.stimuli{i_sti};
        %%
        
        sim_data = NFkB_para_vec_sim(sim_info_para,data_info,sim_data);
    end
end

%,sim_data.parameter_i_para(:)
sim_table_data = table(sim_data.parameter_module(:),sim_data.parameter_reac_num(:),sim_data.parameter_para_num(:),sim_data.parameter_name(:),...
    sim_data.parameter_fold_change(:),sim_data.parameter_value(:),...
    sim_data.ligand(:),sim_data.dose_str(:),...
    sim_data.dose_val(:),sim_data.species(:),sim_data.flag(:),sim_data.type(:),sim_data.trajectory(:),'VariableNames',fieldnames(sim_data));

current_folder = pwd;

if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
    save(strcat(data_info.save_file_path,data_info.save_file_name),'sim_data','sim_table_data');
end

cd(current_folder)
