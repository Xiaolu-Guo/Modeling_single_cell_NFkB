% TNFo_dual_para_sample_main
% clear all
function [] = para_sample_sim_sti_doses(data_save_file_path,Num_sample,ligand_str)
data_info.save_file_path = data_save_file_path;

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

switch ligand_str
    case 'TNF'
        index_end = length(proj_num_vec);
        dose_val_ = 10.^linspace(-1,1,21);
        for i_dose =  1:length(dose_val_)
            proj_num_vec{index_end+i_dose} = 18;
            proj_ligand_vec{index_end+i_dose} = {'TNF'};
            proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
            proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
        end
        
    case 'LPS'
        index_end = length(proj_num_vec);
        dose_val_ = 10.^linspace(0,2,21);
        for i_dose =  1:length(dose_val_)
            proj_num_vec{index_end+i_dose} = 25;
            proj_ligand_vec{index_end+i_dose} = {'LPS'};
            proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
            proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
        end
        
    case 'CpG'
        index_end = length(proj_num_vec);
        dose_val_ = 10.^linspace(1,3,21);
        for i_dose =  1:length(dose_val_)
            proj_num_vec{index_end+i_dose} = 26;
            proj_ligand_vec{index_end+i_dose} = {'CpG'};
            proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'nM')};
            proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
        end
        
    case 'polyIC'
        index_end = length(proj_num_vec);
        dose_val_ = 10.^linspace(4,5,21);
        for i_dose =  1:length(dose_val_)
            proj_num_vec{index_end+i_dose} = 27;
            proj_ligand_vec{index_end+i_dose} = {'polyIC'};
            proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)/1000),'ug/mL')};
            proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
        end
        
    case 'Pam3CSK'
        index_end = length(proj_num_vec);
        dose_val_ = 10.^linspace(1,3,21);
        for i_dose =  1:length(dose_val_)
            proj_num_vec{index_end+i_dose} = 28;
            proj_ligand_vec{index_end+i_dose} = {'Pam3CSK'};
            proj_dose_str_vec{index_end+i_dose} = {strcat(num2str(dose_val_(i_dose)),'ng/mL')};
            proj_dose_val_vec{index_end+i_dose} = {dose_val_(i_dose)};
        end
end


gene_info.gene_parameter_value_vec_genotype = cell(0);
for i_proj_num_vec = 1:length(proj_num_vec)
    
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            [para_val,est]  = sample_est_parameters(proj_num_vec{i_proj_num_vec},Num_sample,'../SAEM_proj_2022/');
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';
        case 2
            % sample bi-stimulation
            [para_sample_2ligand,estimates_2ligand] = sample_est_parameters_2receptor(proj_num_vec{i_proj_num_vec},Num_sample*2,'../SAEM_proj_2022/');
            gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_sample_2ligand(:,estimates_2ligand.non_NFkB_index)';
            gene_info.species_name_vec= {{'NFkB'}};
            gene_info.species_value_vec_genotype{1}{1} = para_sample_2ligand(:,estimates_2ligand.NFkB_index)';
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
    
    genotype_sim_save(sim_info,data_info,gene_info);
    
end


