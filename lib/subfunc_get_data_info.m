function data_info = subfunc_get_data_info(data_info,sheetname)
% ./NFkB_common/generate_SAEM_data_filename
% ./NFkB_data/data_SAEM_format
% Ade's data: stimuli_info, 
% sheetname = 'Ade_exp_data';



for i_sti_num_vec = 1: length(data_info.ligand_vec)
    if isfield(data_info,'dose_str_vec')
    stim_info{i_sti_num_vec}.default_cell_num  =...
               get_stimuli_info_num_cells(data_info.ligand_vec{i_sti_num_vec},data_info.dose_str_vec{i_sti_num_vec},sheetname);
 % need to update other files!!
    elseif isfield(data_info,'dose_val')
       
        get_stimuli_info_num_cells(data_info.ligand_vec{i_sti_num_vec},data_info.dose_val{i_sti_num_vec},sheetname);
    else
        error('dose info is not specified!')
    end
        % num_cells_each_dose_all{strcmp(ligand_all,data_info.ligand_vec{i_sti_num_vec})}{data_info.dose_index_vec{i_sti_num_vec}};
    
end

for i_sti_num_vec = 1: length(data_info.ligand_vec)
    if isfield(data_info,'cell_num_vec') && (length(data_info.cell_num_vec)>=i_sti_num_vec) && ~isempty(data_info.cell_num_vec{i_sti_num_vec})
    else
        data_info.cell_num_vec{i_sti_num_vec} = stim_info{i_sti_num_vec}.default_cell_num;
    end
    if isfield(data_info,'SAEM_ending_time_mins_vec') && (length(data_info.SAEM_ending_time_mins_vec)>=i_sti_num_vec) && ~isempty(data_info.SAEM_ending_time_mins_vec{i_sti_num_vec})
    else
        data_info.SAEM_ending_time_mins_vec{i_sti_num_vec} = 690;
    end
end
end