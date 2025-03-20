function dose_scale = get_stimuli_info_dose_scale(ligand,sheetname)

if nargin<2
    stimuli_info_tbl = get_stimuli_info_tbl();
else
    stimuli_info_tbl = get_stimuli_info_tbl(sheetname);    
end

dose_scale = stimuli_info_tbl.dose_scale(stimuli_info_tbl.Ligand == ligand);
dose_scale = 1/dose_scale(1);