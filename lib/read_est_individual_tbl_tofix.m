function [parameters,para_name,estimates] = read_est_individual_tbl(proj_num,data_parent_path)


% function data = read_monolix_data_to_matrix(data_proj_num_vec,data_parent_path,cell_num_data)
% input
% data path, how many doses and how many cells in each dose
% are included in the predictions.txt

% output
% {i_data} is the i_data-th dose
% data.pred_SAEM{i_data}
% data.pred_mode{i_data}
% data.pred_mean{i_data}
% data.exp{i_data}
% data.exp_SAEM_filter_nan{i_data}
% etc.

%% initialization

sim_data_tbl = SimDataTblInitialize();

i_data =1;

data.pred_SAEM=cell(2,1);
data.pred_mode=cell(2,1);
data.pred_mean=cell(2,1);
data.exp=cell(2,1);
data.time_points_num=cell(2,1);
data.exp_SAEM_filter_nan=cell(2,1);
data.pred_SAEM_filter_nan=cell(2,1);
data.exp_mode_filter_nan=cell(2,1);
data.pred_mode_filter_nan=cell(2,1);
data.exp_mean_filter_nan=cell(2,1);
data.pred_mean_filter_nan=cell(2,1);
data.order=cell(2,1);
data.info_ligand= cell(2,1);
data.info_dose_index= cell(2,1);
data.info_dose_str= cell(2,1);
data.info_num_cells = cell(2,1);

%% get the data
for i_proj = 1:length(data_proj_num_vec)
    
    proj_num = data_proj_num_vec(i_proj);
    
    [proj_path,proj_name] = subfunc_get_proj_path(proj_num,data_parent_path);
    % the proj_path is the relative path to Postdoct projects
    clear est_info_proj
    load(strcat(proj_path,'ProjectInfo.mat')) %est_info_proj
    exp_info_initialization
    
    clear info_dose_str_vec cell_num_vec
    for i_sti = 1:length(est_info_proj.data.ligand_vec)
        dose_str_all = get_stimuli_info_dose_str_all(est_info_proj.data.ligand_vec{i_sti});
        info_dose_str_vec{i_sti} = string( dose_str_all(est_info_proj.data.dose_index_vec{i_sti}));
        if ~isfield(est_info_proj.data,'cell_num_vec')
            cell_num_vec(i_sti) = get_stimuli_info_num_cells(est_info_proj.data.ligand_vec{i_sti},info_dose_str_vec{i_sti})
        else
            cell_num_vec(i_sti)  = est_info_proj.data.cell_num_vec{i_sti};
        end
    end
    
    results_filepath= strcat(proj_path,'results/');
    
    estimates=readtable(strcat(results_filepath,'predictions.txt'));
    
    id_end=cumsum(cell_num_vec);
    id_start=1+[0,id_end(1:end-1)];
    
    time_points=max(estimates.time)/5+1;
    
    
    for i_sti=1:length(est_info_proj.data.ligand_vec)
        sim_data_tbl.ligand(i_data,1) = est_info_proj.data.ligand_vec{i_sti};
        sim_data_tbl.dose_str
        sim_data_tbl.dose_val
        sim_data_tbl.species
        sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) = 0; 
        
        % data.info_ligand{i_data} = est_info_proj.data.ligand_vec{i_sti};
        data.info_dose_index{i_data} = est_info_proj.data.dose_index_vec{i_sti};
        data.info_dose_str{i_data} = info_dose_str_vec{i_sti} ;
        data.info_num_cells{i_data} = cell_num_vec(i_sti);
        
        cell_num=id_end(i_sti)-id_start(i_sti)+1;
        data.pred_SAEM{i_data}=NaN(cell_num,time_points);
        data.pred_mode{i_data}=NaN(cell_num,time_points);
        data.pred_mean{i_data}=NaN(cell_num,time_points);
        data.exp{i_data}=NaN(cell_num,time_points);
        
        i_cell=1;
        
        for i_cell_each_sti=id_start(i_sti):id_end(i_sti)
            %length_celli=sum(estimates.id==i);
            time_celli=estimates.time(estimates.id==i_cell_each_sti);
            data.pred_SAEM{i_data}(i_cell,time_celli/5+1)=estimates.indivPred_SAEM(estimates.id==i_cell_each_sti);
            data.pred_mode{i_data}(i_cell,time_celli/5+1)=estimates.indivPred_mode(estimates.id==i_cell_each_sti);
            data.pred_mean{i_data}(i_cell,time_celli/5+1)=estimates.indivPred_mean(estimates.id==i_cell_each_sti);
            data.exp{i_data}(i_cell,time_celli/5+1)=estimates.Y(estimates.id==i_cell_each_sti);
            i_cell=i_cell+1;
        end
        data.time_points_num{i_data}=max(estimates.time(estimates.id<=id_end(i_sti)&estimates.id>=id_start(i_sti)))/5+1;
        
        data.exp_SAEM_filter_nan{i_data}=data.exp{i_data}(~sum(isnan(data.pred_SAEM{i_data}),2),:);
        data.pred_SAEM_filter_nan{i_data}=data.pred_SAEM{i_data}(~sum(isnan(data.pred_SAEM{i_data}),2),:);
        data.exp_mode_filter_nan{i_data}=data.exp{i_data}(~sum(isnan(data.pred_mode{i_data}),2),:);
        data.pred_mode_filter_nan{i_data}=data.pred_mode{i_data}(~sum(isnan(data.pred_mode{i_data}),2),:);
        data.exp_mean_filter_nan{i_data}=data.exp{i_data}(~sum(isnan(data.pred_mean{i_data}),2),:);
        data.pred_mean_filter_nan{i_data}=data.pred_mean{i_data}(~sum(isnan(data.pred_mean{i_data}),2),:);
        
        [~, data.order{i_data}]=sort(max(data.exp{i_data},[],2),'descend');
        i_data = i_data +1;
    end
    
end

data_field_name={'pred_SAEM','pred_mode','pred_mean','exp',...
    'exp_SAEM_filter_nan','pred_SAEM_filter_nan',...
    'exp_mode_filter_nan','pred_mode_filter_nan',...
    'exp_mean_filter_nan','pred_mean_filter_nan'};%...


if nargin==3
    for i_data=1:length(data.exp)
        for i_field_name = 1:length(data_field_name)
            data.(data_field_name{i_field_name}){i_data} = data.(data_field_name{i_field_name}){i_data}(cell_num_data(1):cell_num_data(2),:);
        end
        [~, data.order{i_data}]=sort(max(data.exp{i_data},[],2),'descend');
    end
end


end
% fnames={'reading_table'} ;			% Name of folder
% estimates=readtable(sprintf('%s/predictions.txt',fnames{1}));%

