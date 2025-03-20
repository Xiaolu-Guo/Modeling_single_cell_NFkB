clear all
load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
data_savepath = '../raw_data2023/2024/';
load(strcat(data_savepath,'Sim18_wt_IkBo_codon_metric_wt_revision_1115.mat'),'data','metrics','collect_feature_vects');

collect_feature_vects_wt = collect_feature_vects;
data_wt = data;
metrics_wt = metrics;

load(strcat(data_savepath,'Sim18_wt_IkBo_codon_metric_ikbamut_revision_1115.mat'),'data','metrics','collect_feature_vects');

collect_feature_vects_mut = collect_feature_vects;
data_mut = data;
metrics_mut = metrics;
%index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
% need to be revised if change
index_IkBamut = [1,2,4];
index_WT = [1,2,4];

metrics_cal = cell(1,length(metrics));
metric_names = fieldnames(metrics{1});

i_data = 1;
for i_index_wt = 1:length(index_WT)
    data_info.info_ligand{i_data} = data_wt.info_ligand{index_WT(i_index_wt)};
    data_info.info_dose_str{i_data} = data_wt.info_dose_str{index_WT(i_index_wt)};
    data_info.data_label{i_data} = 'wt';
    % metrics_cal{i_data} = metrics_wt{index_WT(i_index_wt)};

    for i_metric_name = 1:length(metric_names)
        metrics_cal{i_data}.(metric_names{i_metric_name}) = metrics_wt{index_WT(i_index_wt)}.(metric_names{i_metric_name})(1:9:end,:);
    end

    i_data = i_data+1;
end

for i_index_mut = 1:length(index_IkBamut)
    data_info.info_ligand{i_data} = data_mut.info_ligand{index_IkBamut(i_index_mut)};
    data_info.info_dose_str{i_data} = data_mut.info_dose_str{index_IkBamut(i_index_mut)};
    data_info.data_label{i_data} = 'IkBamut';
    % metrics_cal{i_data} = metrics_mut{index_IkBamut(i_index_mut)};
    for i_metric_name = 1:length(metric_names)
        metrics_cal{i_data}.(metric_names{i_metric_name}) = metrics_mut{index_IkBamut(i_index_mut)}.(metric_names{i_metric_name})(1:9:end,:);
    end
    i_data = i_data+1;
end

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled


data_save_path = './data_signaling_codons_Machine_learning_format/';

% Check if the directory exists
if ~exist(data_save_path, 'dir')
    % If the directory does not exist, create it
    mkdir(data_save_path);
end


index_WT = [1,2,3];
index_IkBamut = [4,5,6];

if 1   
    % ligand_all = {'TNF','LPS','PolyIC'}; this is the ligand order
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
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
        
         writematrix(nfkb_codon_all,strcat(data_save_path,'SimIkBamutClass_X_codon_stim_',data_name{i_data_set},'.csv'));
         writematrix(nfkb_id_all,strcat(data_save_path,'SimIkBamutClass_y_codon_stim_',data_name{i_data_set},'.csv'));

    end
end


