% clear all
% tested 05/15/2024
xls_file_save = 1;
unstim_file_save = 1;
SAEM_save = 1;
unstim_ligand_metric_save = 1;

addpath('../NFkB_common/')
addpath('./lib/')

rescale_data = 1;


%% the path to save the data
data_save_path = '../raw_data2023/SAEM_format/';

params.save_path = data_save_path;

if ~exist(params.save_path,'dir')
    mkdir(params.parentpath,params.folder_name)
end
params.exp_original_data_file='Exp_data_Ade';
params.exp_info='exp_info_initialization_Ade';


%% rescaling the data and save

if rescale_data
    % the macrophage cell volume is 6 pl, measured by Stefanie Luecke
    NFkB_max_range = [0.25, 1.75]/6;
    
    % We assume cells in response to LPS 100ng can reach
    % the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
    % (Pam3CSK might not be able to reach the highest Nuc NFkB conc)
    
    ligand_index = 1; %TNF:1  CpG:2  LPS:3  Pam3CSK:4  PolyIC:5
    
    dose_index = [1,2,3]; %LPS 100ng
    
    ligand_index = 3; %TNF:1  CpG:2  LPS:3  Pam3CSK:4  PolyIC:5
    
    dose_index = [5]; %LPS 100ng
    %   notes:
    %     {'100pg','1ng','10ng'} %TNF
    %     {'10nM','33nM','100nM','333nM','1uM'} %CpG
    %     {'1ng','3ng','10ng','33ng','100ng'} %LPS
    %     {'10ng','100ng','1ug'} %P3CSK
    %     {'10ug','33ug','100ug'}  %PolyIC
    
    % 2022 version
    %     rescale_factor = data_rescaling_factor(NFkB_max_range,ligand_index,dose_index,params);
    % new version
    rescale_factor = data_rescaling_factor_2023(NFkB_max_range,ligand_index,dose_index,params);
    
    % run(params.exp_info)
    
    if xls_file_save
        for ligand_index=5 %1;5
            
            if isfield(params,'sheetname')
                stimuli_info_tbl = get_stimuli_info_tbl(params.sheetname);
            else
                stimuli_info_tbl = get_stimuli_info_tbl();
            end
            
            num_dose_all = sum((stimuli_info_tbl.ADM == ligand_index));
            
            
            for dose_index = 1:num_dose_all
                data_rescale(rescale_factor,ligand_index,dose_index,params);
                
            end
        end
        
    end
    
end

%% non-stim data
% data_rescale(rescale_factor,ligand_index,dose_index,params)
% updated on 11/24/2021
% rescale the data using the rescaling factor
if unstim_file_save
    load(params.exp_original_data_file);
    
    sti_info.ligand_name = 'none';
    sti_info.dose = 'x0';
    sti_info.data_name_rescale = 'non_stim.xls';
    %% recale all the data
    
    index = ((dataTbl.Ligand == sti_info.ligand_name) & (dataTbl.Dose == sti_info.dose));
    data = dataTbl.time_series(index,:);
    
    % for debugging
    %histogram(max(data,[],2));hold on
    
    data_rescale = data * rescale_factor;
    
    % for debugging
    % histogram(max(data_rescale,[],2));hold on
    
    
    if isfolder(params.save_path)
    else
        mkdir(params.save_path)
    end
    
    if isfile(strcat(params.save_path,sti_info.data_name_rescale))
        delete(strcat(params.save_path,sti_info.data_name_rescale));
    end
    
    writematrix(data_rescale,strcat(params.save_path,sti_info.data_name_rescale));
    
end


%% SAEM data format

if 1
    %params.TimePtsWeight
    
    % sti_vec.ligand_index_vec = [1,1,1];
    % sti_vec.dose_index_vec = [1,2,3];
    sti_ligand = 'TNF';
    data_info.ligand_vec = {sti_ligand ,sti_ligand ,sti_ligand ,sti_ligand ,sti_ligand  };
    data_info.ligand_vec = {sti_ligand ,sti_ligand ,sti_ligand };
    
    
    data_info.dose_index_vec = {1,2,3,4,5};
    data_info.dose_index_vec = {1,2,3};
    data_info.cell_num_vec = {50,50,50};% optional, default values is all the cells;
    % can only deal with the same amount or default number of cells
    data_info.SAEM_ending_time_mins_vec = {480,480,480,480,480}; %min  % optional, default values is 690 mins
    data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min  % optional, default values is 690 mins
    
    % 690 min is the data length. here we chopped to 8 hours.
    
    
    % sti_vec.ligand_index_vec = [1,1,2];
    % sti_vec.dose_index_vec = [2,3,1];
    % sti_vec.cell_num_vec = [30,30,30]; %or not define
    % sti_vec.SAEM_ending_time_mins_vec = [480,480,480];
    
    %sti_vec.
    
    
    sheetname = 'Ade_exp_data';
    
    params.write_SAEM_data_filename = generate_SAEM_data_filename(data_info,sheetname);
    
    params.write_SAEM_data_filename = strcat(params.save_path,params.write_SAEM_data_filename);
    
    params.end_time_mins = 690;%the default ending time
    params.SAEM_ending_time_mins = 480;
    params.time_each_frame = 5;
    
    if SAEM_save
        data_SAEM_format(data_info,params)
        
    end
    
end
% data_rescale = data_rescale(1:30,:);% only first 30 lines
%% visulization data

%         figure(liagand_index)
%         histogram(max(data_rescale,[],2));hold on
%     figure(liagand_index)
%     title(strcat(ligand_name,'(',string(unit(1)),')'))
%     legend(strrep(dose,'_','.'))
%
%     set(gca,'FontSize',20,'FontWeight','b')


%% non-stim data
% data_rescale(rescale_factor,ligand_index,dose_index,params)
% updated on 11/24/2021
% rescale the data using the rescaling factor
if unstim_ligand_metric_save
    load('Exp_data_Ade');
    clear data metric_structs collect_feature_vects
    data.exp = cell(2,1);
    data.info_ligand = cell(2,1);
    data.info_dose_str = cell(2,1);
    
    ligand_all = { 'TNF';
        'CpG';
        'LPS'
        'Pam3CSK';
        'polyIC';
        'none'};
    
    
    dose_all = {{'x0_1','x1','x10'} %TNF
        {'x33','x100','x333'} %CpG
        {'x1','x3_3','x10'} %LPS
        {'x10','x100','x1000'} %Pam3CSK4
        {'x10','x33','x100'} %polyIC
        {'x0'}};
    
    dose_str_all= {{'100pg/mL','1ng/mL','10ng/mL'} %TNF
        {'33nM','100nM','333nM'} %CpG
        {'1ng/mL','3.3ng/mL','10ng/mL'} %LPS
        {'10ng/mL','100ng/mL','1000ng/mL'} %P3CSK
        {'10ug','33ug','100ug'}
        {'0'}}; %polyIC
    
    % recale all the data
    i_data =1;
    for i_ligand = 1:length(ligand_all)
        for i_dose = 1:length(dose_all{i_ligand})
            index = ((dataTbl.Ligand == ligand_all{i_ligand}) & (dataTbl.Dose == dose_all{i_ligand}{i_dose}));
            data_sti = dataTbl.time_series(index,:);
            
            % for debugging
            %histogram(max(data,[],2));hold on
            data_sti(any(isnan(data_sti), 2), :) = [];
            sum(sum(isnan(data_sti)))
            data.exp{i_data} = data_sti * rescale_factor;
            data.info_ligand{i_data} = replace(ligand_all{i_ligand},'polyIC','PolyIC');
            data.info_dose_str{i_data} = dose_str_all{i_ligand}{i_dose};
            i_data = i_data+1;
        end
        % for debugging
        % histogram(max(data_rescale,[],2));hold on
        
    end
    vis_data_field = {'exp'};
    % data_label = {'exp'};%,'sample'};

    [collect_feature_vects,metric_structs] = calculate_codon(data,vis_data_field );%,data_label
    save(strcat(data_save_file_path,'Ade_all_stim_unstim_codon.mat'),'data','metric_structs','collect_feature_vects');
end
