function rescale_factor = data_rescaling_factor(NFkB_max_range,ligand_index,dose_index,params)

% We assume cells in response to LPS 100ng can reach
% the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
% (Pam3CSK might not be able to reach the highest Nuc NFkB conc)

%% define the rescale factor: based on highest NFkB response to LPS

load(params.exp_original_data_file);

sti_info = get_stimulus_info(ligand_index,dose_index,params.exp_info);

% data_index_all = (dataTbl.Ligand == sti_info.ligand_name);

index = ((dataTbl.Ligand == sti_info.ligand_name) & (dataTbl.Dose == sti_info.dose));
data = dataTbl.time_series(index,:);

rescale_factor = NFkB_max_range(2) / prctile(max(data,[],2),90) ;% highest dose of LPS