function [] = draw_all_ligand_codon_distrib_2023_03(data_save_file_path,fig_save_path)
% to be deleted 
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

fig_opt.paper_opt.paperpos=[0,0,150,180]*3;
fig_opt.paper_opt.papersize=[150 180]*3;

    % data_proj_num_vec = data_proj_nums{i_ligand};    

    load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
       
for i_sti = 1:length(data.pred_mode_filter_nan)
    
    data.pred_mode_cv{i_sti} = std(data.pred_mode_filter_nan{i_sti},[],2)./mean(data.pred_mode_filter_nan{i_sti},2);
    [~,data.pred_mode_cv_order{i_sti}] = sort(data.pred_mode_cv{i_sti},'descend');
    [~,data.pred_mode_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2},'descend');
    [~,data.exp_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2-1},'descend');

end

% input
    thresh_TNF=0.33;
    thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
for i_data = 1:3
    index = data.(thresh_field{1}){i_data} (1: ceil(length(data.(thresh_field{1}){i_data}) * thresh_TNF));
    data.pred_mode{i_data} = data.pred_mode_filter_nan{i_data}(index ,:);
    data.exp{i_data} = data.exp_mode_filter_nan{i_data}(index,:);
    [~, data.order{i_data}]=sort(max(data.exp{i_data},[],2),'descend');
end

% codon_list = {'Duration','EarlyVsLate','OscVsNonOsc','PeakAmplitude','Speed','TotalActivity'};



metrics_list = fieldnames(metrics{1});

for i_metrics =1:length(metrics_list)
    for i_data = 1:3
        index = data.(thresh_field{1}){i_data} (1: ceil(length(data.(thresh_field{1}){i_data}) * thresh_TNF));
        
        metrics{i_data*2-1}.(metrics_list{i_metrics}) = metrics{i_data*2-1}.(metrics_list{i_metrics})( index,:);
        metrics{i_data*2}.(metrics_list{i_metrics}) = metrics{i_data*2}.(metrics_list{i_metrics})( index,:);
    end
end

ids = 1:length(collect_feature_vects.Duration);

codewords = {'Duration'     ;
    'EarlyVsLate'  ;
    'OscVsNonOsc'  ;
    'PeakAmplitude';
    'Speed'        ;
    'TotalActivity'};


CodewordList = {'Duration';
    'EarlyVsLate'  ;
    'OscVsNonOsc'  ;
    'PeakAmplitude';
    'PeakAmplitude';
    'Speed'        ;
    'Speed'        ;
    'Speed'        ;
    'TotalActivity'};

FeatureSpecifiers = {[ 2];
    [-1];
    [ 0];
    [ 0];
    [ 0];
    [ 0];
    [-1];
    [ 2];
    [ 0]};

FeatureList = {'duration'    ;
    'time2HalfMaxPosIntegral';
    'oscpower'               ;
    'max_value'              ;
    'pos_pk1_amp'            ;
    'max_pos_pk1_speed'      ;
    'pos_pk1_time'           ;
    'derivatives'            ;
    'max_pos_integral'       };

for i = 1:length(codewords)
    vects = cell(1, length(ids));
    for j= 1:length(FeatureList)
        if strcmp(CodewordList{j}, codewords{i})
            for k=1:length(ids)
                metric_struct = metrics{k};
                feature = metric_struct.(FeatureList{j});
                if abs(FeatureSpecifiers{j}) ~= 0
                    feature = feature(:, abs(FeatureSpecifiers{j}));
                end
                if FeatureSpecifiers{j} < 0
                    feature = feature*-1;
                end
                if parameters.calcResponders
                    feature = feature(logical(metric_struct.responder_index));
                end
                vects{k} = [vects{k} feature];
            end
        end
    end
    if size(vects{1}, 2) > 1  %%% for codewords that are composed of more than one metric-->take z score and then take average
        vects_zscore=vects(1:2:end);
        all_data_zscore = cell2mat(vects_zscore(:));
        mean_zscore= nanmean(all_data_zscore,1);
        std_zscore=nanstd(all_data_zscore);
        all_data = cell2mat(vects(:));
        all_data = (all_data-mean_zscore)./std_zscore;
        
        count = 1;
        for k = 1:length(ids)
            vects{k} = nanmean(all_data(count:(size(vects{k}, 1)+count-1), :), 2);
            count = count + size(vects{k}, 1);
        end
    end
    all_data = cell2mat(vects(:)); %%% normalize each codeword
    all_data = (all_data-prctile(all_data,0.5))/(prctile(all_data,99.5)-prctile(all_data,0.5)); %normalized except those extrema
    count = 1;
    for k = 1:length(ids)
        vects{k} = all_data(count:(size(vects{k}, 1)+count-1), :);
        count = count + size(vects{k}, 1);
    end
    collect_feature_vects.(codewords{i}) = vects';
end

save(  strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'),'metrics','collect_feature_vects','data')% data,collect_feature_vects,metrics));


% change to metric




        %% codon distribution
    % violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
    % codon and sti seperately, and save
    fig_opt.save_file = strcat(fig_save_path,'All_ligand_codon_2023_distrib');
%    violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
   violin_plot_codon_2023_03(collect_feature_vects,fig_opt) % draw all codon and sti in one

 % violin_plot_codon_2023_CRI(collect_feature_vects,fig_opt)



