function write_data_filename = generate_SAEM_data_filename_down_sampling_20241001(data,sheetname)

stimuli_info_tbl = get_stimuli_info_tbl(sheetname);
for i_ligand_vec = 1:length(data.ligand_vec)
    dose_vals = stimuli_info_tbl.dose_val...
        (stimuli_info_tbl.Ligand == data.ligand_vec{i_ligand_vec});
    
    data.dose_val{i_ligand_vec} = dose_vals(i_ligand_vec);
end

% data = subfunc_get_data_info(data,sheetname);



% defined ligand_all, dose_str_all, num_dose_all
sti_vec.ligand_index_vec = [];

%sti_vec.ligand = data.

for ii = 1:length(data.ligand_vec)
    switch data.ligand_vec{ii}
        case 'TNF'
            sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 1];
        case 'CpG'
            sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 2];
        case 'Pam3CSK'
            sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 4];
        case 'LPS'
            sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 3];
        case 'polyIC'
            sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 5];
    end
    
end

sti_vec.dose_index_vec = cell2mat( data.dose_index_vec);

if isfield(data,'cell_num_vec')
else
    stimuli_info_tbl = get_stimuli_info_tbl(sheetname);
    for sti_index = 1:length(sti_vec.dose_index_vec)
        index = (stimuli_info_tbl.Ligand == char(data.ligand_vec{sti_index})) & (stimuli_info_tbl.dose_index == data.dose_index_vec{sti_index});
        data.cell_num_vec{sti_index} = stimuli_info_tbl.num_cells(index);
    end
end
sti_vec.cell_num_vec = cell2mat(data.cell_num_vec ) ;
sti_vec.SAEM_ending_time_mins_vec  =cell2mat(data.SAEM_ending_time_mins_vec);

write_data_filename = subfunc_data_filename(sti_vec,sheetname);

write_data_filename = strcat('SAEM_data_',write_data_filename,'.txt');

end



function write_filename = subfunc_data_filename(sti_vec,sheetname)

stimuli_info_tbl = get_stimuli_info_tbl(sheetname);
ligand_all_cat = unique(stimuli_info_tbl.Ligand,'stable');
ligand_all = categorical2cellstr(ligand_all_cat);


write_filename = '';

for ligand_index=1:length(ligand_all)
    
    % find the corresponding ligand name
    sti_indx = find(sti_vec.ligand_index_vec==ligand_index);
    
    if ~isempty(sti_indx)
        
        % if all doses included, name the file only with ligand name
        % otherwise, specifiy the doses in the file name
        write_filename = strcat(write_filename,ligand_all{sti_vec.ligand_index_vec(sti_indx(1))});
        
        dose_str_ve = categorical2cellstr(stimuli_info_tbl.dose_fold(stimuli_info_tbl.Ligand==ligand_all_cat(ligand_index)));% {sti_vec.ligand_index_vec(sti_indx(1))};
        num_cell_ve = stimuli_info_tbl.num_cells(stimuli_info_tbl.Ligand==ligand_all_cat(ligand_index));%num_cells_each_dose_all{sti_vec.ligand_index_vec(sti_indx(1))};
        
        if length(sti_indx) ~= length(num_cell_ve)
            
            
            for i_sti_indx=1:length(sti_indx)
                
                
                write_filename = strcat(write_filename,dose_str_ve{sti_vec.dose_index_vec(sti_indx(i_sti_indx))});
                
                num_cells = num_cell_ve(sti_vec.dose_index_vec(sti_indx(i_sti_indx)));
                if isfield(sti_vec,'cell_num_vec') && (sti_vec.cell_num_vec(sti_indx(i_sti_indx)) < num_cells)
                    write_filename = strcat(write_filename,num2str(sti_vec.cell_num_vec(sti_indx(i_sti_indx))),'C');
                end
                
                if isfield(sti_vec,'SAEM_ending_time_mins_vec')
                    write_filename = strcat(write_filename,num2str(sti_vec.SAEM_ending_time_mins_vec(sti_indx(i_sti_indx))),'min');
                end
                
            end
        else
            num_cells = num_cell_ve(sti_vec.dose_index_vec(sti_indx(1)));
            if isfield(sti_vec,'cell_num_vec') && (sti_vec.cell_num_vec(sti_indx(1)) < num_cells)
                write_filename = strcat(write_filename,num2str(sti_vec.cell_num_vec(sti_indx(1))),'C');
            end
            
            if isfield(sti_vec,'SAEM_ending_time_mins_vec')
                write_filename = strcat(write_filename,num2str(sti_vec.SAEM_ending_time_mins_vec(sti_indx(1))),'min');
            end
        end
        
        
        
        write_filename = strcat(write_filename,'_');
    end
end

if write_filename(end) == '_'
    write_filename = write_filename(1:end-1);
end

end