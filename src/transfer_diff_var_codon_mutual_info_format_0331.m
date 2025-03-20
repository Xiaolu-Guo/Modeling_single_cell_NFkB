function [] = transfer_diff_var_codon_mutual_info_format_0331(data_save_file_path,MI_file_save_path)
% For Guo et al. simulation data for Figure ??
% run simulation for unstimulated case: to get Sim_unstim_fitting_alldose_r2_codon_metric.mat
% 
% tested 05/15/2024, Matlab 2020a


%% NFkB red var, 10 groups, vs wt, vs IkB-/-.
if 1
    file_list = {'Sim17_reduce_core_heterogeneity_codon_metric1.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric2.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric3.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric4.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric5.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric6.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric7.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric8.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric9.mat'
        'Sim17_reduce_core_heterogeneity_codon_metric10.mat'
        'Sim18_wt_IkBo_codon_metric1.mat'
        'Sim18_wt_IkBo_codon_metric2.mat'
        'Sim18_wt_IkBo_codon_metric3.mat'
        'Sim5_SS_codon_metric_p25x_r1.mat'
        'Sim19_reduce_rcpt_heterogeneity_codon_metric1.mat'
        'Sim20_reduce_TAK1ac_heterogeneity_codon_metric1.mat'
        'Sim21_reduce_TAK1ac_NFkB_heterogeneity_codon_metric1.mat'
        'Sim21_reduce_TAK1ac_NFkB_heterogeneity_codon_metric2.mat'
        'Sim22_TAK1_heterogeneity_codon_metric1.mat'
        'Sim23_NFkB_heterogeneity_codon_metric1.mat'
        'Sim24_no_heterogeneity_codon_metric1.mat'
        'Sim19_IkBas_reduce_rcpt_heterogeneity_codon_metric1.mat'};
    
    data_labels = {'NFkB_var_red1',...
        'NFkB_var_red2',...
        'NFkB_var_red3',...
        'NFkB_var_red4',...
        'NFkB_var_red5',...
        'NFkB_var_red6',...
        'NFkB_var_red7',...
        'NFkB_var_red8',...
        'NFkB_var_red9',...
        'NFkB_var_red10',...
        'wt1',...
        'IkBo2',...
        'IkBo3',...
        'IkBop25',...
        'rcpt_red',...
        'TAKac_red',...
        'TAKac_NFkB_red1',...
        'TAKac_NFkB_red2',...
        'TAK_noise',...
        'NFkB_noise',...
        'no_noise',...
        'IkBo_rcpt_red'};
    
    %% load data
    % 1 reduce NFkB var to representative cell paras
    for i_data = 1:length(file_list)
        load(strcat(data_save_file_path,file_list{i_data}))
        metrics_all{i_data} = metrics;
        data_all{i_data} = data;
    end
    
    %% cal codon and sim vs exp diff
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    i_index = 1;
    clear metrics_cal data_info
    field_list = fieldnames(metrics_all{1}{1});
    for i_cond = 1:length(metrics_all{1})
        
        for i_data = 1:length(file_list)
            for i_field = 1:length(field_list)
                
                metrics_cal{i_index}.(field_list{i_field}) = metrics_all{i_data}{i_cond}.(field_list{i_field})(1:9:end,:);
            end
            
            data_info.info_ligand(i_index) = [data_all{i_data}.info_ligand(i_cond)];
            data_info.info_dose_str(i_index) =[data_all{i_data}.info_dose_str(i_cond)];
            data_info.data_label(i_index) = data_labels(i_data);
            
            i_index = i_index+1;
        end
        
    end
    
    [collect_feature_vects_all_minmax_scaled,~] = calculate_codon_from_metric2023(data_info,metrics_cal); %,  parameter
    
    data_name = data_labels;
    for i_data = 1:length(data_labels)
        
        ligand_all = {'TNF', 'LPS', 'CpG', 'PolyIC', 'Pam3CSK'};
        codon_list = {'TotalActivity', 'Duration', 'EarlyVsLate', 'Speed', 'PeakAmplitude', 'OscVsNonOsc'};
        
        for i_ligand = 1:length(ligand_all)
            
            nfkb(i_ligand).sc_metrics = struct();
                        
            for i_codon =1:length(codon_list)
                nfkb(i_ligand).id = ligand_all{i_ligand};
                nfkb(i_ligand).ids = ligand_all;
                nfkb(i_ligand).sc_metrics.(codon_list{i_codon}) = collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(i_ligand-1)*length(data_labels)+i_data};
                
            end
            
        end
        
        save(strcat(MI_file_save_path,'mutual_info_format_alldata_sc_',data_name{i_data},'_min_max_rescale_0331.mat'),'nfkb')
        
    end
end
end
