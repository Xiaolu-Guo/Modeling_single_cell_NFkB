function cell_num = get_stimuli_info_num_cells(ligand,dose_info,sheetname)

if nargin<3
    stimuli_info_tbl = get_stimuli_info_tbl();
else
    stimuli_info_tbl = get_stimuli_info_tbl(sheetname);
end

if isnumeric(dose_info)
    cell_num = stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == ligand) & (stimuli_info_tbl.dose_val == dose_info));
else
    cell_num = stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == ligand) & (stimuli_info_tbl.dose_str == dose_info));
end