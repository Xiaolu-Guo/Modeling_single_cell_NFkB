% est_info_proj.model  
% 
% model.est_para_receptor = {'params53','params58'};  
% model.est_para_core = {'params52n2','params99','params100'}; 
% model.est_receptor = 'TNF';
% model.est_para_default  = {'NFkB_cyto_init','shift'};
% est_info_proj.model.sti_lag  = 0 ;
% 
% model.est_para = {'params52n2','params99','params100','params53','params58','NFkB_cyto_init','shift'};
% 
% model.filename = 'params52n2_params99_params101_params54_params58_NFkB_cyto_init_shift';



% est_info_proj.est.para_range_set = {[],[],[],[],[],[],[]};
% est_info_proj.est.distribution_type = {[],[],[],[],[],[],[]};
% est_info_proj.est.data_path_relative_to_proj  = '../../Experiments/202202_rescaled_byXiaolu/';

setpath_datafile

exp_info.sheetname = 'Ade_exp_data_2023'; % Talor_pub_model % Stefanie_tnfo_data

sti_ligand = 'TNF';

switch sti_ligand
    
    case 'LPS'
        data_info.ligand_vec = {sti_ligand,sti_ligand,sti_ligand,sti_ligand,sti_ligand};
        data_info.dose_index_vec = {1,2,3,4,5};
        data_info.SAEM_ending_time_mins_vec = {480,480,480,480,480}; %min
        est_info_proj.model.est_para_receptor = {'params35','params44','params36'};
        data_info.filename = 'SAEM_data_LPS480min.txt';
        proj_num = 3;
        
    case 'CpG'
        data_info.ligand_vec = {sti_ligand,sti_ligand,sti_ligand,sti_ligand,sti_ligand};
        data_info.dose_index_vec = {1,2,3,4,5};
        data_info.SAEM_ending_time_mins_vec = {480,480,480,480,480}; %min
        est_info_proj.model.est_para_receptor = {'params85','params93','params88'};
        data_info.filename = 'SAEM_data_CpG480min.txt';
        proj_num = 4;

    case 'PolyIC'     % to be fixed    
        data_info.ligand_vec = {sti_ligand,sti_ligand,sti_ligand};
        data_info.dose_index_vec = {1,2,3};
        % data_info.cell_num_vec = {100,100,100};
        data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min
        est_info_proj.model.est_para_receptor = {'params77','params83','params79'};
        data_info.filename = 'SAEM_data_PolyIC480min.txt';
        proj_num = 5;
        
    case 'TNF'
        data_info.ligand_vec = {sti_ligand,sti_ligand,sti_ligand};
        data_info.dose_index_vec = {1,2,3};
        % data_info.cell_num_vec = {100,100,100};
        data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min
        est_info_proj.model.est_para_receptor = {'params54','params58'};
        data_info.filename = 'SAEM_data_TNF480min.txt';
        proj_num = 2;

    case 'Pam3CSK'
        data_info.ligand_vec = {sti_ligand,sti_ligand,sti_ligand};
        data_info.dose_index_vec = {1,2,3};
        % data_info.cell_num_vec = {100,100,100};
        data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min
        est_info_proj.model.est_para_receptor = {'params68','params75'};
        data_info.filename = 'SAEM_data_Pam3CSK480min.txt';
                proj_num = 6;

end




% model info
est_info_proj.model.est_para_core = {'params52n2','params99','params101'};
est_info_proj.model.est_receptor = sti_ligand;

est_info_proj.model.est_para_default = {'NFkB_cyto_init','shift'};

%%%%% remember to fix the receptor sens_para setting!!!

est_info_proj.model.est_para = [est_info_proj.model.est_para_core(:)', est_info_proj.model.est_para_receptor(:)',est_info_proj.model.est_para_default(:)'];%


% estimation info
% both core and receptor! first core then receptor!
est_info_proj.est.para_range_set = cell(size(est_info_proj.model.est_para));
% est_info_proj.est.min_val = {'','0.1','','',''};
% est_info_proj.est.max_val = {'','0.1','','',''};

est_info_proj.est.distribution_type = cell(size(est_info_proj.model.est_para));
% est_info_proj.est.para_distribution_type = {'logitNormal','logitNormal','logitNormal','logitNormal','logitNormal'};
% est_info_proj.model.para_default_range_set = {[];[]};
% est_info_proj.model.para_default_distribution_type = {'logitNormal','logNormal'};

% 'logitNormal', normal? lognormal?


%% 0. name the project [done]

est_info_proj.data = subfunc_data_info_generate(data_info,exp_info.sheetname);


est_info_proj.proj_name = strcat('XGESN2023',sprintf( '%03d', proj_num  ));
est_info_proj.est.data_path_relative_to_proj = data_path_relative_to_proj;

%% 1. create folder [done]

est_info_proj.proj_path = '../../NFkB_para_estm_project/SAEM_proj_2023/';

est_info_proj.proj_path = strcat(est_info_proj.proj_path,est_info_proj.proj_name,'/');

%% 2. create model file [done]

est_info_proj.model.filename = generate_model_filename(est_info_proj);

save('ProjectInfo','est_info_proj')

function data_info = subfunc_data_info_generate(data,sheetname)


stimuli_info_tbl = get_stimuli_info_tbl(sheetname);
ligand_all_cat = unique(stimuli_info_tbl.Ligand,'stable');
ligand_all = categorical2cellstr(ligand_all_cat);


data_info = data;
data_info.stim = data.ligand_vec;

for ii = 1: length(data.ligand_vec)
    ligand_index = find(strcmp(data.ligand_vec{ii},ligand_all));
    
    dose_str_all = categorical2cellstr(stimuli_info_tbl.dose_fold(stimuli_info_tbl.Ligand==ligand_all_cat(ligand_index)));% {sti_vec.ligand_index_vec(sti_indx(1))};
dose_val_all = stimuli_info_tbl.dose_val(stimuli_info_tbl.Ligand==ligand_all_cat(ligand_index));
    data_info.dose_val{ii} = dose_val_all(data.dose_index_vec{ii});
    data_info.dose_str{ii} = dose_str_all{data.dose_index_vec{ii}};
    
end

end

function  filename = generate_model_filename(est_info_proj)

filename = '';

for ii = 1:length(est_info_proj.model.est_para)
    
    if isempty(est_info_proj.est.para_range_set{ii})
        filename = strcat(filename,est_info_proj.model.est_para{ii},'_');
    else
        filename = strcat(filename,est_info_proj.model.est_para{ii},...
            num2str(est_info_proj.est.para_range_set{ii}(1)),'to',...
            num2str(est_info_proj.est.para_range_set{ii}(2)),'_');
    end
    
end

filename = filename(1:end-1);

end
