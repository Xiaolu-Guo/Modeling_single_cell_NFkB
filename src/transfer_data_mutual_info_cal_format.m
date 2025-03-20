function [] = transfer_data_mutual_info_cal_format(data_save_file_path,MI_file_save_path)
% For Guo et al. Figure 5A, extended doses of each ligand simulated heatmaps
% 
% tested 05/15/2024, Matlab 2020a

load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
load(strcat(data_save_file_path,'All_codon_dual_ligand.mat'))

index_supriya = 1:5;
index_ade = 6:10;
index_fitting = 11:15;
index_sampling = 16:20;
index_supriya_dual = 21:29;
index_sampling_dual = [30:33,35:39];


if 1   
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    data_name = {'Ade','Supriya','Fitting','Sampling'};
    index_vec = {index_ade,index_supriya,index_fitting,index_sampling};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            
            for i_codon =1:length(codon_list)
                nfkb(i_ligand).id = ligand_all{i_ligand};
                nfkb(i_ligand).ids = ligand_all;
                nfkb(i_ligand).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){index_data(i_ligand)};
            end
        end
        save(strcat(MI_file_save_path,'mutual_info_format_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end

