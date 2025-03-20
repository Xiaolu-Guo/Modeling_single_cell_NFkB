clear nfkb collect_feature_vects
load('./example_data_format/mutual_info_cal_data_example.mat')
nfkb = nfkb(1:3);

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

load(strcat(data_save_file_path_1,'Ade_all_stim_unstim_codon.mat'));%All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));




if 1
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    
    index_vec = {[1,2,3,16],...%TNF
        [10:12,16],... %Pam3CSK
        [4:6,16],... % CpG
        [7:9,16],... % LPS
        [13:15,16]};%PolyIC
    dose_all = {'Low','Medium','High','Unstim'};
    for i_ligand = 1:length(ligand_all)
        index_data = index_vec{i_ligand};

        for i_index_data = 1:length(index_data)
            nfkb(i_index_data).sc_metrics = struct();
            for i_codon =1:length(codon_list)
                nfkb(i_index_data).id = dose_all{i_index_data};
                nfkb(i_index_data).ids = dose_all;
                nfkb(i_index_data).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){index_data(i_index_data)};
            end
        end
        save(strcat('mutual_info_format_ade_exp_',ligand_all{i_ligand},'.mat'),'nfkb')
        
    end
end


