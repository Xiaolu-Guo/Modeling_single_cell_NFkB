function dose_str_all = get_stimuli_info_dose_str_all(ligand,sheetname)

if nargin<2
    stimuli_info_tbl = get_stimuli_info_tbl();
else
    stimuli_info_tbl = get_stimuli_info_tbl(sheetname);    
end

dose_str_all = stimuli_info_tbl.dose_str(stimuli_info_tbl.Ligand == ligand);
