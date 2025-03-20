% TNFo_dual_para_sample_main
% clear all
function filename_list = para_sample_sim_namelist(data_save_file_path)
data_info.save_file_path = data_save_file_path;
Num_sample =10;
% the paramters that will be sampled

gene_info.gene_type = {'wt'};

proj_num_vec = cell(0);

% proj_num_vec = {18,25,26,27,28};% ,[18,26],[18,25],[26,27]}
% proj_ligand_vec = {{'TNF'},{'LPS'},{'CpG'},{'polyIC'},{'Pam3CSK'}};% ,{'TNF','CpG'},{'TNF','LPS'},{'CpG','polyIC'}}
% proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};%
% proj_dose_val_vec = {{10},{10},{100},{100000},{100}};

% proj_num_vec = {[18,26],[18,25],[26,27]};
% proj_ligand_vec = {{'TNF','CpG'},{'TNF','LPS'},{'CpG','polyIC'}};
% proj_dose_str_vec = {{'10ng/mL','100nM'},{'10ng/mL','10ng/mL'},{'100nM','100ug/mL'}};%
% proj_dose_val_vec = {{10,100},{10,10},{100,100000}};

proj_num_vec ={[5]
    [7]
    [7,3]
    [3]
    [3,5]
    [4]
    [6]
    [3,4]
    [3,6]
    [5]
    [7]
    [6]
    [6,5]
    [6,7]
    [6]
    [6,5]
    [6]
    [6,5]
    [6,4]
    [6,3,5]
    [6,3,7]
    [5]
    [7]
    [4]
    [4,5]
    [4,7]
    [7]};
proj_ligand_vec = {{'CpG'}
    {'Pam3CSK'}
    {'Pam3CSK','TNF'}
    {'TNF'}
    {'TNF','CpG'}
    {'LPS'}
    {'polyIC'}
    {'TNF','LPS'}
    {'TNF','polyIC'}
    {'CpG'}
    {'Pam3CSK'}
    {'polyIC'}
    {'polyIC','CpG'}
    {'polyIC','Pam3CSK'}
    {'polyIC'}
    {'polyIC','CpG'}
    {'polyIC'}
    {'polyIC','CpG'}
    {'polyIC','LPS'}
    {'polyIC','TNF','CpG'}
    {'polyIC','TNF','Pam3CSK'}
    {'CpG'}
    {'Pam3CSK'}
    {'LPS'}
    {'LPS','CpG'}
    {'LPS','Pam3CSK'}
    {'Pam3CSK'}};
proj_dose_str_vec = {{'100nM'}
    {'100ng'}
    {'100ng','10ng'}
    {'10ng'}
    {'10ng','100nM'}
    {'10ng'}
    {'100ug'}
    {'10ng','10ng'}
    {'10ng','100ug'}
    {'100nM'}
    {'100ng'}
    {'100ug'}
    {'100ug','100nM'}
    {'100ug','100ng'}
    {'100ug'}
    {'100ug','100nM'}
    {'100ug'}
    {'100ug','100nM'}
    {'100ug','10ng'}
    {'100ug','10ng','100nM'}
    {'100ug','10ng','100ng'}
    {'100nM'}
    {'100ng'}
    {'10ng'}
    {'10ng','100nM'}
    {'10ng','100ng'}
    {'100ng'}};
proj_dose_val_vec = {{100}
    {100}
    {100,10}
    {10}
    {10,100}
    {10}
    {100000}
    {10,10}
    {10,100000}
    {100}
    {100}
    {100000}
    {100000,100}
    {100000,100}
    {100000}
    {100000,100}
    {100000}
    {100000,100}
    {100000,10}
    {100000,10,100}
    {100000,10,100}
    {100}
    {100}
    {10}
    {10,100}
    {10,100}
    {100}};
%
% % TNF
% index_end = length(proj_num_vec);
% dose_val_ = 10.^linspace(-1,1,21);
% for i_dose =  1:length(dose_val_)
%     proj_num_vec{index_end+i_dose} = 18;
%     proj_ligand_vec{index_end+i_dose} = {'TNF'};
%     proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
%     proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
% end
%
% % LPS
% index_end = length(proj_num_vec);
% dose_val_ = 10.^linspace(0,2,21);
% for i_dose =  1:length(dose_val_)
%     proj_num_vec{index_end+i_dose} = 25;
%     proj_ligand_vec{index_end+i_dose} = {'LPS'};
%     proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
%     proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
% end
%
% % CpG
% index_end = length(proj_num_vec);
% dose_val_ = 10.^linspace(1,3,21);
% for i_dose =  1:length(dose_val_)
%     proj_num_vec{index_end+i_dose} = 26;
%     proj_ligand_vec{index_end+i_dose} = {'CpG'};
%     proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'nM')};
%     proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
% end
%
% % polyIC
% index_end = length(proj_num_vec);
% dose_val_ = 10.^linspace(4,5,21);
% for i_dose =  1:length(dose_val_)
%     proj_num_vec{index_end+i_dose} = 27;
%     proj_ligand_vec{index_end+i_dose} = {'CpG'};
%     proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)/1000),'ug/mL')};
%     proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
% end
%
% % Pam3CSK
% index_end = length(proj_num_vec);
% dose_val_ = 10.^linspace(1,3,21);
% for i_dose =  1:length(dose_val_)
%     proj_num_vec{index_end+i_dose} = 28;
%     proj_ligand_vec{index_end+i_dose} = {'Pam3CSK'};
%     proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
%     proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
% end


gene_info.gene_parameter_value_vec_genotype = cell(0);
for i_proj_num_vec = 1:length(proj_num_vec)
    
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            [~,est]  = sample_est_parameters_2(proj_num_vec{i_proj_num_vec},Num_sample,'../SAEM_proj_2022_2/');
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(non_NFkB_index)};
        case 2
            % sample bi-stimulation
            [~,estimates_2ligand] = sample_est_parameters_2receptor_2(proj_num_vec{i_proj_num_vec},Num_sample*2,'../SAEM_proj_2022_2/');
            gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
    end
    
    % stimili info
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    save_filename = 'para_sampled_NFkBinit';
    
    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end
    
    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'}};
    data_info.save_file_name = save_filename; % only beginning, no .mat
    
    
    for i_para_vec = 1:length(gene_info.parameter_name_vec)
        
        sim_info.parameter_name = gene_info.parameter_name_vec{i_para_vec};
        
        %para_val = gene_info.para_val_vec{i_para_vec};
        
        
        if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
            save_file_name = data_info.save_file_name;
            for i_para_name = 1:length(sim_info.parameter_name)
                save_file_name = strcat(save_file_name,'_',sim_info.parameter_name{i_para_name});
            end
            filename_list{i_proj_num_vec} = strcat(save_file_name,'.mat');
        end
        
    end
    
    
    
end

end


