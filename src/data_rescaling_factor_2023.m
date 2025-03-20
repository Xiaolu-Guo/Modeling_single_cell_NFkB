function rescale_factor = data_rescaling_factor_2023(NFkB_max_range,ligand_index,dose_index,params)

% We assume cells in response to LPS 100ng can reach
% the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
% (Pam3CSK might not be able to reach the highest Nuc NFkB conc)

%% define the rescale factor: based on highest NFkB response to LPS

load(params.exp_original_data_file);

% ligand_index = 4; dose_index = 5;
% sti_info = get_stimulus_info(ligand_index,dose_index,params.exp_info);
clear sti_info
        sti_info.ligand_name = 'LPS';
             sti_info.num_dose= 5;
           sti_info.dose_scale= 24000;
          sti_info.dose_ligand= {'x1'  'x3_3'  'x10'  'x33'  'x100'};
      sti_info.dose_str_ligand= {'1ng'  '3ng'  '10ng'  '33ng'  '100ng'};
      sti_info.dose_val_ligand= {[1]  [3]  [10]  [33]  [100]};
     sti_info.num_cells_ligand= {[176]  [245]  [437]  [310]  [327]};
         sti_info.ligand_index= 4;
                  sti_info.ADM= 4;
                 sti_info.dose= 'x100';
             sti_info.dose_str= '100ng';
             sti_info.dose_val= 100;
            sti_info.num_cells= 327;
    sti_info.data_name_rescale= 'LPS_100ng.xls';
    
% data_index_all = (dataTbl.Ligand == sti_info.ligand_name);

index = ((dataTbl.Ligand == sti_info.ligand_name) & (dataTbl.Dose == sti_info.dose));
data = dataTbl.time_series(index,:);

rescale_factor = NFkB_max_range(2) / prctile(max(data,[],2),90) ;% highest dose of LPS