% called in rescaling factor and data_rescale_truncate
function sti_info = get_stimulus_info(ligand_index,dose_index,exp_info_func)

%load the experimental information
%see initialization.m

if nargin <3
    exp_info_initialization_Ade
else
    run(exp_info_func)
end
%% stimuli information corresponding to ligand
if ischar(ligand_index) || isstring(ligand_index)
    sti_info.ligand_name = ligand_index;
    ligand_index = find(strcmp(ligand_all,sti_info.ligand_name));
elseif isnumeric(ligand_index)
    sti_info.ligand_name = ligand_all{ligand_index};
end

sti_info.num_dose = num_dose_all{ligand_index};

sti_info.dose_scale = dose_scale_all{ligand_index};

sti_info.dose_ligand = dose_all{ligand_index};

sti_info.dose_str_ligand = dose_str_all{ligand_index};

sti_info.dose_val_ligand = dose_val_all{ligand_index};

sti_info.num_cells_ligand = num_cells_each_dose_all{ligand_index};

sti_info.ligand_index = ligand_index;

sti_info.ADM = sti_info.ligand_index;


%% stimuli information corresponding to speicific dose
sti_info.dose = sti_info.dose_ligand{dose_index};

sti_info.dose_str = sti_info.dose_str_ligand{dose_index};

sti_info.dose_val = sti_info.dose_val_ligand{dose_index};

sti_info.num_cells = sti_info.num_cells_ligand{dose_index};

sti_info.data_name_rescale = strcat(sti_info.ligand_name,'_',sti_info.dose_str,'.xls');
end