clear all
load('mutual_info_cal_data_example.mat')


data_filename = 'Sim5_codon_all5dose_metric.mat';
load(strcat('../raw_data2023/',data_filename))
% TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
% Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
% CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
% LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
% PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15
sample_index = [13;12;11;14;15];

%% % transfer NFkB conds
if 0 % transfer time points: 0 bits, not working
    
    data_info.info_ligand = data.info_ligand(sample_index);
    data_info.info_dose_str = data.info_dose_str(sample_index);
    data_info.data_label = {'sample_sim','sample_sim','sample_sim','sample_sim','sample_sim'};
    metric_fields = fieldnames(metrics{1});
    for i_metric_index = 1:length(sample_index)
        i_metric = sample_index(i_metric_index);
        for i_fields = 1:length(metric_fields)
            metrics_cal{i_metric_index}.(metric_fields{i_fields}) = metrics{i_metric}.(metric_fields{i_fields})(:,:);
        end
    end
    
    [collect_feature_vects,~] = calculate_codon_from_metric2023(data_info,metrics_cal); %,  parameter
    
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    timept_list ={};
    for i_tpt = 1:size(data.model_sim{1},2)
        timept_list{i_tpt} = strcat('timept_',num2str((i_tpt-1)*5),'min');
    end
    data_name = {'Sampling'};
    index_vec = {sample_index};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            
            for i_codon =1:length(codon_list)
                nfkb(i_ligand).id = ligand_all{i_ligand};
                nfkb(i_ligand).ids = ligand_all;
                nfkb(i_ligand).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){i_ligand}(1:9:end);
            end
        end
        save(strcat('mutual_info_format_NFkB_codon_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end


%% % transfer NFKB time points: 0 bits, not working
if 0 % transfer time points: 0 bits, not working
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    timept_list ={};
    for i_tpt = 1:size(data.model_sim{1},2)
        timept_list{i_tpt} = strcat('timept_',num2str((i_tpt-1)*5),'min');
    end
    data_name = {'Sampling'};
    index_vec = {sample_index};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            
            for i_time_point =1:length(timept_list)
                nfkb(i_ligand).id = ligand_all{i_ligand};
                nfkb(i_ligand).ids = ligand_all;
                nfkb(i_ligand).sc_metrics.(timept_list{i_time_point}) = data.model_sim{index_data(i_ligand)}(1:9:end,i_time_point);
            end
        end
        save(strcat('mutual_info_format_NFKB_time_traj_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end

%% % transfer IKK TAK time points: 0 bits, not working
if 0 % transfer time points: 0 bits, not working
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    timept_list ={};
    for i_tpt = 1:size(data.model_sim{1},2)
        timept_list{i_tpt} = strcat('timept_',num2str((i_tpt-1)*5),'min');
    end
    data_name = {'Sampling'};
    index_vec = {sample_index};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            
            for i_time_point =1:length(timept_list)
                nfkb(i_ligand).id = ligand_all{i_ligand};
                nfkb(i_ligand).ids = ligand_all;
                nfkb(i_ligand).sc_metrics.(timept_list{i_time_point}) = data.model_sim{index_data(i_ligand)}(7:9:end,i_time_point);
            end
        end
        save(strcat('mutual_info_format_IKK_time_traj_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end

%% transfer into TAK1 codon, metrics, etc
if 0 % transfer into TAK1 codon, metrics, etc
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    timept_list ={};
    threshold_duration = 0.00005:0.00005:0.001;
    
    dynamic_feature = cell(1);
    for i_ligand = 1:length(sample_index)
        i_sti = sample_index(i_ligand);
        dynamic_feature{i_ligand}.integrals = diff([zeros(size(metrics{i_sti}.integrals,1)/9,1),metrics{i_sti}.integrals(8:9:end,7:6:end)],1,2);
        for i_duration = 1:length(threshold_duration)
            dynamic_feature{i_ligand}.duration(:,i_duration) = sum(metrics{i_sti}.time_series(8:9:end,:)>threshold_duration(i_duration),2);% unit is 5min
        end
        %dynamic_feature{i_sti}.first_half_int() = ;
        %dynamic_feature{i_sti}.second_half_int() = ;
        dynamic_feature{i_ligand}.max_val = metrics{i_sti}.max_amplitude(8:9:end);% 8 is TAK1, 7 is IKK
        dynamic_feature{i_ligand}.speed = metrics{i_sti}.time2HalfMaxIntegral(8:9:end);
    end
    
    for i_integral = 1:size(dynamic_feature{1}.integrals,2)
        integral_list{i_integral} = strcat('integral_',num2str((i_integral)*30),'min');
    end
    
    for i_duration = 1:size(dynamic_feature{1}.duration,2)
        duration_list{i_duration} = strcat('duration_',num2str((i_duration)),'th');
    end
    
    feature_list ={'max_val','speed'};
    
    
    for i_integral =1:length(integral_list)
        data_total = [];
        for i_ligand = 1:length(dynamic_feature)
            data_total = [data_total;dynamic_feature{i_ligand}.integrals(:,i_integral)];
        end
        dynamic_zscore.mean.integrals(:,i_integral) = nanmean(data_total);
        dynamic_zscore.std.integrals(:,i_integral) = nanstd(data_total);
    end
    
    for i_duration = 1:length(duration_list)
        data_total = [];
        
        for i_ligand = 1:length(dynamic_feature)
            
            data_total = [data_total;dynamic_feature{i_ligand}.duration(:,i_duration)];
        end
        dynamic_zscore.mean.duration(:,i_duration) = nanmean(data_total);
        dynamic_zscore.std.duration(:,i_duration) = nanstd(data_total);
        
    end
    
    for i_feature = 1:length(feature_list)
        data_total = [];
        for i_ligand = 1:length(dynamic_feature)
            
            data_total = [data_total;dynamic_feature{i_ligand}.(feature_list{i_feature})];
        end
        dynamic_zscore.mean.(feature_list{i_feature}) = nanmean(data_total);
        dynamic_zscore.std.(feature_list{i_feature}) = nanstd(data_total);
    end
    
    data_name = {'Sampling'};
    index_vec = {sample_index};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            nfkb(i_ligand).id = ligand_all{i_ligand};
            nfkb(i_ligand).ids = ligand_all;
            % normalize or not?
            for i_integral =1:length(integral_list)
                nfkb(i_ligand).sc_metrics.(integral_list{i_integral}) =...
                    (dynamic_feature{i_ligand}.integrals(:,i_integral)...
                    - dynamic_zscore.mean.integrals(:,i_integral))...
                    /dynamic_zscore.std.integrals(:,i_integral);
            end
            
            for i_duration = 1:length(duration_list)
                nfkb(i_ligand).sc_metrics.(duration_list{i_duration}) =...
                    (dynamic_feature{i_ligand}.duration(:,i_duration)...
                    -dynamic_zscore.mean.duration(:,i_duration))...
                    /  dynamic_zscore.std.duration(:,i_duration);
                
            end
            
            for i_feature = 1:length(feature_list)
                nfkb(i_ligand).sc_metrics.(feature_list{i_feature}) =...
                    (dynamic_feature{i_ligand}.(feature_list{i_feature})...
                    -dynamic_zscore.mean.(feature_list{i_feature}))...
                    /dynamic_zscore.std.(feature_list{i_feature});
                
            end
        end
        save(strcat('mutual_info_format_TAK1_DynFeature_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end

%% transfer into IKK codon, metrics, etc
if 0 % transfer into IKK codon, metrics, etc
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    timept_list ={};
    threshold_duration = 0.001:0.001:0.03;
    
    dynamic_feature = cell(1);
    for i_ligand = 1:length(sample_index)
        i_sti = sample_index(i_ligand);
        dynamic_feature{i_ligand}.integrals = diff([zeros(size(metrics{i_sti}.integrals,1)/9,1),metrics{i_sti}.integrals(7:9:end,7:6:end)],1,2);
        for i_duration = 1:length(threshold_duration)
            dynamic_feature{i_ligand}.duration(:,i_duration) = sum(metrics{i_sti}.time_series(7:9:end,:)>threshold_duration(i_duration),2);% unit is 5min
        end
        %dynamic_feature{i_sti}.first_half_int() = ;
        %dynamic_feature{i_sti}.second_half_int() = ;
        dynamic_feature{i_ligand}.max_val = metrics{i_sti}.max_amplitude(7:9:end);% 8 is TAK1, 7 is IKK
        dynamic_feature{i_ligand}.speed = metrics{i_sti}.time2HalfMaxIntegral(7:9:end);
    end
    
    for i_integral = 1:size(dynamic_feature{1}.integrals,2)
        integral_list{i_integral} = strcat('integral_',num2str((i_integral)*30),'min');
    end
    
    for i_duration = 1:size(dynamic_feature{1}.duration,2)
        duration_list{i_duration} = strcat('duration_',num2str((i_duration)),'th');
    end
    
    feature_list ={'max_val','speed'};
    
    
    for i_integral =1:length(integral_list)
        data_total = [];
        for i_ligand = 1:length(dynamic_feature)
            data_total = [data_total;dynamic_feature{i_ligand}.integrals(:,i_integral)];
        end
        dynamic_zscore.mean.integrals(:,i_integral) = nanmean(data_total);
        dynamic_zscore.std.integrals(:,i_integral) = nanstd(data_total);
    end
    
    for i_duration = 1:length(duration_list)
        data_total = [];
        
        for i_ligand = 1:length(dynamic_feature)
            
            data_total = [data_total;dynamic_feature{i_ligand}.duration(:,i_duration)];
        end
        dynamic_zscore.mean.duration(:,i_duration) = nanmean(data_total);
        dynamic_zscore.std.duration(:,i_duration) = nanstd(data_total);
        
    end
    
    for i_feature = 1:length(feature_list)
        data_total = [];
        for i_ligand = 1:length(dynamic_feature)
            
            data_total = [data_total;dynamic_feature{i_ligand}.(feature_list{i_feature})];
        end
        dynamic_zscore.mean.(feature_list{i_feature}) = nanmean(data_total);
        dynamic_zscore.std.(feature_list{i_feature}) = nanstd(data_total);
    end
    
    data_name = {'Sampling'};
    index_vec = {sample_index};
    
    for i_data_set = 1:length(data_name)
        index_data = index_vec{i_data_set};
        
        for i_ligand = 1:length(ligand_all)
            nfkb(i_ligand).sc_metrics = struct();
            nfkb(i_ligand).id = ligand_all{i_ligand};
            nfkb(i_ligand).ids = ligand_all;
            % normalize or not?
            for i_integral =1:length(integral_list)
                nfkb(i_ligand).sc_metrics.(integral_list{i_integral}) =...
                    (dynamic_feature{i_ligand}.integrals(:,i_integral)...
                    - dynamic_zscore.mean.integrals(:,i_integral))...
                    /dynamic_zscore.std.integrals(:,i_integral);
            end
            
            for i_duration = 1:length(duration_list)
                nfkb(i_ligand).sc_metrics.(duration_list{i_duration}) =...
                    (dynamic_feature{i_ligand}.duration(:,i_duration)...
                    -dynamic_zscore.mean.duration(:,i_duration))...
                    /  dynamic_zscore.std.duration(:,i_duration);
                
            end
            
            for i_feature = 1:length(feature_list)
                nfkb(i_ligand).sc_metrics.(feature_list{i_feature}) =...
                    (dynamic_feature{i_ligand}.(feature_list{i_feature})...
                    -dynamic_zscore.mean.(feature_list{i_feature}))...
                    /dynamic_zscore.std.(feature_list{i_feature});
                
            end
        end
        save(strcat('mutual_info_format_IKK_DynFeature_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
        
    end
end
