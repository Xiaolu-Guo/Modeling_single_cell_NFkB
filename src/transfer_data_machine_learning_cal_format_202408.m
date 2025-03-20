 clear all
load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
load(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'))

index_supriya = 1:5;
index_ade = 6:10;
index_fitting = 11:15;
index_sampling = 16:20;
index_supriya_dual = 21:29;
index_sampling_dual = [30:33,35:39];

data_save_path = './data_signaling_codons_Machine_learning_format/';

% Check if the directory exists
if ~exist(data_save_path, 'dir')
    % If the directory does not exist, create it
    mkdir(data_save_path);
end

if 1   
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    data_name = {'Ade','Fitting','Sampling'};
    index_vec = {index_ade,index_fitting,index_sampling};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        nfkb_codon_all = [];
        nfkb_id_all = [];
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            nfkb_codon = [];
            nfkb_id = [];
            for i_codon =1:length(codon_list)
                nfkb_codon = [nfkb_codon,collect_feature_vects.(codon_list{i_codon}){index_data(i_ligand)}];
            end
            nfkb_id = (i_ligand-1)*ones(size(nfkb_codon,1),1);
            nfkb_codon_all = [nfkb_codon_all;nfkb_codon];
            nfkb_id_all = [nfkb_id_all;nfkb_id];
           [~, order_data] = sort(max(metrics{index_data(i_ligand)}.time_series,[],2),'descend')
            % h=heatmap(metrics{index_data(i_ligand)}.time_series(order_data,:) ,'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF

        end
        
         writematrix(nfkb_codon_all,strcat(data_save_path,'X_codon_stim_',data_name{i_data_set},'.csv'));
         writematrix(nfkb_id_all,strcat(data_save_path,'y_codon_stim_',data_name{i_data_set},'.csv'));

    end
end


