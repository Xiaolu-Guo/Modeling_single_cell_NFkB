%% put all data in one file

% sti_vec.ligand_index_vec = [1,1,1];
% sti_vec.dose_index_vec = [1,2,3];
generate_all_data = 0;
generate_TNF_P3CSK_data = 1;
generate_single_data = 0;
addpath('../NFkB_common/')

if generate_all_data
    sti_ligand = 'TNF';
    data_info.ligand_vec = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec = {1,2,3};
    data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min  % optional, default values is 690 mins
    
    sti_ligand = 'LPS';
    data_info.ligand_vec(4:8) = {sti_ligand, sti_ligand, sti_ligand, sti_ligand, sti_ligand};
    data_info.dose_index_vec(4:8) = {1,2,3,4,5};
    data_info.SAEM_ending_time_mins_vec(4:8) = {480,480,480,480,480}; %min  % optional, default values is 690 mins
    
    sti_ligand = 'CpG';
    data_info.ligand_vec(9:13) = {sti_ligand, sti_ligand, sti_ligand, sti_ligand, sti_ligand};
    data_info.dose_index_vec(9:13) = {1,2,3,4,5};
    data_info.SAEM_ending_time_mins_vec(9:13) = {480,480,480,480,480}; %min  % optional, default values is 690 mins
    
    sti_ligand = 'polyIC';
    data_info.ligand_vec(14:16) = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec(14:16) = {1,2,3};
    data_info.SAEM_ending_time_mins_vec(14:16) = {480,480,480}; %min  % optional, default values is 690 mins
    
    sti_ligand = 'Pam3CSK';
    data_info.ligand_vec(17:19) = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec(17:19) = {1,2,3};
    data_info.SAEM_ending_time_mins_vec(17:19) = {480,480,480}; %min  % optional, default values is 690 mins
  
    sheetname = 'Ade_exp_data';
    params.write_SAEM_data_filename = 'SAEM_data_all480min.txt';

end

if generate_TNF_P3CSK_data
    sti_ligand = 'TNF';
    data_info.ligand_vec = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec = {1,2,3};
    data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min  % optional, default values is 690 mins
        
    sti_ligand = 'Pam3CSK';
    data_info.ligand_vec(4:6) = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec(4:6) = {1,2,3};
    data_info.SAEM_ending_time_mins_vec(4:6) = {480,480,480}; %min  % optional, default values is 690 mins
  
    sheetname = 'Ade_exp_data';
    params.write_SAEM_data_filename = 'SAEM_data_TNF_P3CSK_480min.txt';

end

if generate_single_data
    sti_ligand = 'TNF';
    data_info.ligand_vec = {sti_ligand, sti_ligand, sti_ligand };
    data_info.dose_index_vec = {1,2,3};
    data_info.SAEM_ending_time_mins_vec = {480,480,480}; %min  % optional, default values is 690 mins
% data_info.cell_num_vec = {50,50,50};% optional, default values is all the cells;
% can only deal with the same amount or default number of cells

    sheetname = 'Ade_exp_data';

params.write_SAEM_data_filename = generate_SAEM_data_filename(data_info,sheetname);

end



params.SAEM_ending_time_mins = 480;
params.time_each_frame = 5;
params.exp_info='../NFkB_common/exp_info_initialization_Ade';
params.sheetname = sheetname;
params.save_path = '../../NFkB_para_estm_project/Experiments/202202_rescaled_byXiaolu/';
% params.exp_original_data_file='Exp_data_Ade';

data_SAEM_format(data_info,params)