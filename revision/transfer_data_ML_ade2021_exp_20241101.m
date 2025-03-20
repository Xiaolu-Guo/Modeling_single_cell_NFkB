clear all
load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
data_savepath = '../raw_data2023/2024/';
load(strcat(data_savepath,'Ade2021_IkBamut_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');

%index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
% need to be revised if change
index_IkBamut = [1:3];
index_WT = [4:6];


data_save_path = './data_signaling_codons_Machine_learning_format/';

% Check if the directory exists
if ~exist(data_save_path, 'dir')
    % If the directory does not exist, create it
    mkdir(data_save_path);
end

if 1   
    % ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
    
    data_name = {'IkBamut','WT'};
    index_vec = {index_IkBamut,index_WT};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        nfkb_codon_all = [];
        nfkb_id_all = [];
        for i_ligand = 1:length(index_data)% in total length(index_data) conditions 
            nfkb(i_ligand).sc_metrics = struct();
            nfkb_codon = [];
            nfkb_id = [];
            for i_codon =1:length(codon_list)
                nfkb_codon = [nfkb_codon,collect_feature_vects.(codon_list{i_codon}){index_data(i_ligand)}];
            end
            nfkb_id = (i_ligand-1)*ones(size(nfkb_codon,1),1);
            nfkb_codon_all = [nfkb_codon_all;nfkb_codon];
            nfkb_id_all = [nfkb_id_all;nfkb_id];
           %[~, order_data] = sort(max(metrics{index_data(i_ligand)}.time_series,[],2),'descend')
            % h=heatmap(metrics{index_data(i_ligand)}.time_series(order_data,:) ,'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF

        end
        
         writematrix(nfkb_codon_all,strcat(data_save_path,'Ade2021exp_X_codon_stim_',data_name{i_data_set},'.csv'));
         writematrix(nfkb_id_all,strcat(data_save_path,'Ade2021exp_y_codon_stim_',data_name{i_data_set},'.csv'));

    end
end


