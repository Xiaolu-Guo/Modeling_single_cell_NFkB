%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

% vers = '_r2';
% vers_savefig = strcat(vers,'');%_IkBao _IkBao_onlyPeak
file_names = {'Sim16_IkBao_5_signle_ligand_codon_metric_r1','Sim15_5_signle_ligand_codon_metric_r2'};

% Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
ligand_all = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

data_name = {'wt','IkBao'};

for i_data_set = 1:length(data_name)
    load(strcat('../raw_data2023/',file_names{i_data_set},'.mat'))
load('mutual_info_cal_data_example.mat')

    for i_data = 1:length(data.model_sim)
        
        data_info.info_ligand{i_data} = data.info_ligand{i_data};
        data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
        data_info.data_label{i_data} = 'SamplingComb';
        
    end
    metrics_cal = cell(1,length(metrics));
    metric_names = fieldnames(metrics{1});
    
    for i_metric_index = 1:length(metrics)
        for i_metric_name = 1:length(metric_names)
            metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
        end
    end
    
    %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
    [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
    
    
    for i_ligand = 1:length(ligand_all)
        nfkb(i_ligand).sc_metrics = struct();
        
        for i_codon =1:length(codon_list)
            nfkb(i_ligand).id = ligand_all{i_ligand};
            nfkb(i_ligand).ids = ligand_all;
            nfkb(i_ligand).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){i_ligand};
        end
    end
    save(strcat('../raw_data2023/','mutual_info_format_dual_ligand_',data_name{i_data_set},'.mat'),'nfkb')
    
end
