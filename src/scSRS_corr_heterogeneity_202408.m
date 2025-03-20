function [] = scSRS_corr_heterogeneity_202408()
% For Guo et al. Figure 7 S7: single-cell stimulus response specificity
%
% tested 05/15/2024, Matlab 2020a

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';
savepath = fig_save_path;
addpath('./lib/')
addpath('./src/')
addpath('./bin/')
plot_vers = '_0812';%'_0705';

%% initialization and data preparation

if 1
    vers = '_r1';
    vers_savefig = strcat(vers,'_matching_SS_p25x');%_IkBao _IkBao_onlyPeak ,'_SS_0x' _p1x
    
    load(strcat('../raw_data2023/Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x',vers,'.mat'))
    data_IkBo = data;
    para_mat_IkBo = sim_data_tbl.parameter_value(1:9:end,:);
    ligand_index = [1,5,3,2,4]; % reoder the stim
  
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
    
    length_data = length(data.model_sim);
    length_metrics = length(metrics);
    
    load(strcat('../raw_data2023/Sim8_5_signle_ligand_codon_metric_r3','.mat'))
    para_mat_wt = sim_data_tbl.parameter_value(1:9:end,:);
    
    for i_data = 1:length(data.model_sim)
        
        data_info.info_ligand{i_data+length_data} = data.info_ligand{i_data};
        data_info.info_dose_str{i_data+length_data} = data.info_dose_str{i_data};
        data_info.data_label{i_data+length_data} = 'SamplingComb';
        
    end
    
    for i_metric_index = 1:length(metrics)
        for i_metric_name = 1:length(metric_names)
            metrics_cal{i_metric_index+length_metrics}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
        end
    end
    
    
    %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
    [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
    
    codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
    weights_codon = [1,1,1,1,1,1];
    
    
    clear codon_single_cell_wt codon_single_cell_SS
    codon_single_cell_wt=cell(2,1);
    codon_single_cell_IkBo=cell(2,1);
    
    for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
        for i_ligand = 1:length_data
            for i_codon = 1:length(codon_list)
                codon_single_cell_IkBo{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
            end
        end
    end
    
    for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
        for i_ligand = 1:length_data
            for i_codon = 1:length(codon_list)
                codon_single_cell_wt{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand+length_data}(i_cell,:);
            end
        end
    end
    
    
    dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
    dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
    
    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
    
    dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
end

%% [tested] Figure 6A, S6A-B: wild type binary category map
if 1
    
    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
    dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
    
    
    conf_pairs = {'TNF-LPS', 'TNF-CpG', 'TNF-P3C', 'TNF-PIC',... [1,2,4,3]
        'LPS-CpG', 'LPS-P3C', 'CpG-P3C',... [5,7,9]
        'LPS-PIC', 'CpG-PIC', 'P3C-PIC'};% [6,8,10]
    
    
    dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
    color_limit = [0,50;0,20;0.5,1];
    threshold_vec = {1}; % {0.5,0.8,1,1.2,1.5};
    
    
    
    if 1 %Figure 6A S6B WT: tested 08/12/2024
        % to get Figure S6B different threshold results, set
        % threshold_new to different values:
        % 0.5,0.8,1.0,1.2,1.5
        colormap_mat = [0,0,0
            1,1,1];
        color_limit_6A = [0,1];
        threshold_new = 1;
        [fig1] = hierarchical_row_column_plot_kw_0411(dist_mat_wt,threshold_new,color_limit_6A,colormap_mat,conf_pairs)
        
        fig1
        saveas(gcf,strcat(savepath,'fig6_wt_conf_mat_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
        close
        
    end
    
    
    if 1 % Figure S6A: tested 08/12/2024
        for i_thresh = 1:length(threshold_vec)
            threshold_new = threshold_vec{i_thresh}; %0.8 1.2
            
            conf_mat_thresh = dist_mat_wt<=threshold_new;% confusion
            [fig1,fig2,fig3] = calculate_draw_corr_mat(conf_mat_thresh,color_limit);
            
            fig3
            saveas(gcf,strcat(savepath,'fig6_jaccard_dist_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
            close
            
            fig2
            %saveas(gcf,strcat(savepath,'fig6_jaccard_index_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
            close
            
            fig1
            %saveas(gcf,strcat(savepath,'fig6_cell_number_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
            close
            
            
        end
        
        
    end
    
end

%% [tested] figure 6D: barplot of SRS prob different threshold vs cell numbers

if 1 % Figure 6D : tested on 08/12/2024
    % epsilon = 2; % Define the epsilon value
    threshold_vec = {0.5,0.8,1,1.2,1.5};
    
    for i_threshold = 1:length(threshold_vec)
        threshold_new = threshold_vec{i_threshold};
        conf_mat_wt = cell2mat(cellfun(@(x) confusion_mat(x,threshold_new), codon_single_cell_wt, 'UniformOutput', false));
        
        % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
        %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
        %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
        %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
        
        %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
        %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
        %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
        
        conf_mat_wt = conf_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
        
        conf_score_wt = sum(conf_mat_wt,2);
        
        reduce_heterogeneity_conf_mat = zeros(size(conf_mat_wt));
        
        cell_num = 0;
        for cluster_num =0:10
            barplot(cluster_num+1,1) = sum(conf_score_wt == cluster_num);
            reduce_heterogeneity_conf_mat(cell_num+1:cell_num+barplot(cluster_num+1,1),1:cluster_num) = 1;
            cell_num = cell_num+barplot(cluster_num+1,1);
        end
        
        
        % [fig1] = prob_cell_num_SRS(dist_mat_wt,threshold_new);
        cell_num = 10;
        
        
        distinguish_mat = 1-reduce_heterogeneity_conf_mat;
        [prob_vec] = prob_cell_num_SRS(distinguish_mat)  ;
        
        cell_num_95 = find(prob_vec>=0.95,1);
        cell_num_99 = find(prob_vec>=0.99,1);
        
        if isempty(cell_num_95)
            bar_plot_red_heterogeneity_cell_95prob(i_threshold) = 12;
        else
            bar_plot_red_heterogeneity_cell_95prob(i_threshold) = cell_num_95;
        end
        
        if isempty(cell_num_99)
            bar_plot_red_heterogeneity_cell_99prob(i_threshold) = 12;
        else
            bar_plot_red_heterogeneity_cell_99prob(i_threshold) = cell_num_99;
        end
        
        
    end
    
    
    %{'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
    dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
    
    
    conf_pairs = {'TNF-LPS', 'TNF-CpG', 'TNF-P3C', 'TNF-PIC',... [1,2,4,3]
        'LPS-CpG', 'LPS-P3C', 'CpG-P3C',... [5,7,9]
        'LPS-PIC', 'CpG-PIC', 'P3C-PIC'};% [6,8,10]
    
    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
    
    dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
    color_limit = [0,50;0,20;0.5,1];
    threshold_vec = {0.5,0.8,1,1.2,1.5}; % {0.5,0.8,1,1.2,1.5};
    
    for i_threshold = 1:length(threshold_vec)
        threshold_new = threshold_vec{i_threshold}; %0.8 1.2
        
        % [fig1] = prob_cell_num_SRS(dist_mat_wt,threshold_new);
        cell_num = 10;
        
        paperpos=[0,0,100,100];
        
        papersize=[100,100];
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        distinguish_mat = dist_mat_wt > threshold_new;
        [prob_vec] = prob_cell_num_SRS(distinguish_mat) ;
        
        
        cell_num_95 = find(prob_vec>=0.95,1);
        cell_num_99 = find(prob_vec>=0.99,1);
        
        if isempty(cell_num_95)
            bar_plot_wt_cell_95prob(i_threshold) = 12;
        else
            bar_plot_wt_cell_95prob(i_threshold) = cell_num_95;
        end
        
        if isempty(cell_num_99)
            bar_plot_wt_cell_99prob(i_threshold) = 12;
        else
            bar_plot_wt_cell_99prob(i_threshold) = cell_num_99;
        end
        
    end
    
    
    barplot_all_95 = [bar_plot_wt_cell_95prob;bar_plot_red_heterogeneity_cell_95prob];
    barplot_all_99 = [bar_plot_wt_cell_99prob;bar_plot_red_heterogeneity_cell_99prob];
    
    figure(1)
    paperpos=[0,0,100,100];
    
    papersize=[100,100];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
    bar(barplot_all_95')
    ylim([0,10])
    set(gca,'XTickLabels',{},'YTickLabels',{})
    saveas(gcf,strcat(savepath,'fig6_probSRS95_cell_num_0705'),'epsc');%_20cells
    
    close()
    
    figure(1)
    paperpos=[0,0,100,100];
    
    papersize=[100,100];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
    bar(barplot_all_99')
    ylim([0,10])
    
    set(gca,'XTickLabels',{},'YTickLabels',{})
    saveas(gcf,strcat(savepath,'fig6_probSRS99_cell_num_0705'),'epsc');%_20cells
    
    close()
    
    
end

%% [tested] figure 6C: plot probSRS vs cell numbers

if 1 % Figure 6C: tested on 08/12/2024
    
    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
    dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
    
    
    conf_pairs = {'TNF-LPS', 'TNF-CpG', 'TNF-P3C', 'TNF-PIC',... [1,2,4,3]
        'LPS-CpG', 'LPS-P3C', 'CpG-P3C',... [5,7,9]
        'LPS-PIC', 'CpG-PIC', 'P3C-PIC'};% [6,8,10]
    
    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
    
    dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
    color_limit = [0,50;0,20;0.5,1];
    threshold_vec = {1}; % {0.5,0.8,1,1.2,1.5};
    
    threshold_new = [1]; %0.8 1.2
    
    % [fig1] = prob_cell_num_SRS(dist_mat_wt,threshold_new);
    cell_num = 10;
    
    paperpos=[0,0,150,100];
    
    papersize=[150,100];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
    
    distinguish_mat = dist_mat_wt > threshold_new;
    [prob_vec] = prob_cell_num_SRS(distinguish_mat)
    figure(1)
    plot(1:cell_num,prob_vec,'LineWidth',1.5);hold on
    
    conf_mat_wt = cell2mat(cellfun(@(x) confusion_mat(x,threshold_new), codon_single_cell_wt, 'UniformOutput', false));
    
    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
    
    conf_mat_wt = conf_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    
    conf_score_wt = sum(conf_mat_wt,2);
    
    reduce_heterogeneity_conf_mat = zeros(size(conf_mat_wt));
    
    cell_num_count = 0;
    for cluster_num =0:10
        barplot(cluster_num+1,1) = sum(conf_score_wt == cluster_num);
        reduce_heterogeneity_conf_mat(cell_num_count+1:cell_num_count+barplot(cluster_num+1,1),1:cluster_num) = 1;
        cell_num_count = cell_num_count+barplot(cluster_num+1,1);
    end
    
    
    distinguish_mat = 1-reduce_heterogeneity_conf_mat;
    [prob_vec] = prob_cell_num_SRS(distinguish_mat)
    figure(1)
    plot(1:cell_num,prob_vec,'LineWidth',1.5);hold on
    
    ylim([0,1])
    set(gca,'XTick',0:2:10,'YTick',0:0.2:1,'XTickLabels',{},'YTickLabels',{})
    saveas(gcf,strcat(savepath,'fig6_probSRS_cell_num_wt_vs_rh',plot_vers),'epsc');%_20cells
    close
    
end

%% [tested] figure 6B S6B: low diversity


if 1
    % to get different threshold set threshold_new to
    % 0.5,0.8,1.0,1.2,1.5
    threshold_vec = {1};
    threshold_new = 1; %0.5,0.8,1.0,1.2,1.5
    conf_mat_wt = cell2mat(cellfun(@(x) confusion_mat(x,threshold_new), codon_single_cell_wt, 'UniformOutput', false));
    
    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
    
    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
    
    conf_mat_wt = conf_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
    
    conf_score_wt = sum(conf_mat_wt,2);
    
    reduce_heterogeneity_conf_mat = zeros(size(conf_mat_wt));
    
    cell_num = 0;
    for cluster_num =0:10
        barplot(cluster_num+1,1) = sum(conf_score_wt == cluster_num);
        reduce_heterogeneity_conf_mat(cell_num+1:cell_num+barplot(cluster_num+1,1),1:cluster_num) = 1;
        cell_num = cell_num+barplot(cluster_num+1,1);
    end
    
    if 1 % Figure 6B & S6B: tested 08/12/2024
        figure(1)
        
        colormap_mat = [0,0,0
            1,1,1];
        color_limit = [0,1];
        reordered_traj_mat = 1-reduce_heterogeneity_conf_mat;
        paperpos=[0,0,100,200]*3;
        papersize=[100,200]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',colormap_mat,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
        
        XLabels = 1:size(reordered_traj_mat,2);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/1);%conf_pairs
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:size(reordered_traj_mat,1);
        % Convert each number in the array into a string
        YCustomXLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        YCustomXLabels(:) = " ";
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.YDisplayLabels = YCustomXLabels;
        
        % xlabel('Time (hours)');
        % ylabel(vis_data_field{i_data_field});
        % clb=colorbar;
        % clb.Label.String = 'NFkB(A.U.)';
        colorbar('off')
        
        set(gca,'fontsize',14,'fontname','Arial');
        
        saveas(gcf,strcat(savepath,'fig6_red_hetero_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
        close
        
    end
    
    if 1 %Figure S6A low diversity: tested 08/12/2024
        color_limit = [0,50;0,20;0.5,1];
        [fig1,fig2,fig3] = calculate_draw_corr_mat(reduce_heterogeneity_conf_mat,color_limit);
        
        fig3
        saveas(gcf,strcat(savepath,'fig6_jaccard_dist_red_hetero_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
        close
        
        fig2
        %saveas(gcf,strcat(savepath,'fig6_jaccard_index_red_hetero_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
        close
        
        fig1
        %saveas(gcf,strcat(savepath,'fig6_cell_number_red_hetero_conf_pairs_',replace(num2str(threshold_new),'.','p'),plot_vers),'epsc');%_20cells
        close
        
    end
end

end


%% calculated the confusion matrix
function conf_mat = confusion_mat(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3
conf_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        conf_mat(i_col) = (sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2))<epsilon);
        i_col = i_col + 1;
    end
end

end


%% calculate the distance matrix of stimulus pairs
function dist_mat = dist_mat(matrix)

% epsilon = 2; % Set the epsilon value 1,2,3
dist_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        dist_mat(i_col) = sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2));
        i_col = i_col + 1;
    end
end

end


%% hierachical plot row cluster
function [fig3,Outperm_col] = hierarchical_row_column_plot_kw_0411(traj_mat,threshold,color_limit,color_map,conf_pairs,Outperm_col)

traj_mat_all = traj_mat > threshold;
index_0pair = find(sum(traj_mat_all,2) == 10);
traj_mat_0pair = traj_mat_all(index_0pair,:);

index_1pair = find(sum(traj_mat_all,2) == 9);
traj_mat_1pair = traj_mat_all(index_1pair,:);

index_2pair_more = find(sum(traj_mat_all,2) <= 8);
traj_mat_2pair_more = traj_mat_all(index_2pair_more,:);

index_all = index_0pair;
traj_mat_1pair_ordered = [];

for i_pair = 1:size(traj_mat,2)
    index_pair = traj_mat_1pair(:,i_pair) == 0;
    traj_mat_1pair_ordered = [traj_mat_1pair_ordered;traj_mat_1pair(index_pair,:)];
end



figure(1)

% Step 1: Perform hierarchical clustering on the rows
Y = pdist(traj_mat_2pair_more, 'euclidean'); % Compute the pairwise distances between rows
Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(1)
paperpos=[0,0,55,70];
papersize=[55 70];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
close();

% Step 3: Reorder the matrix based on the clustering result
traj_mat_2pair_more_ordered = traj_mat_2pair_more(Outperm, :);

traj_mat_plot = [traj_mat_0pair;traj_mat_1pair_ordered;traj_mat_2pair_more_ordered];

if nargin < 6
    Outperm_col = 1:size(traj_mat,2);
end


figure(3)

paperpos=[0,0,100,200]*3;
papersize=[100,200]*3;

set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

% subplot(1,length(vis_data_field),i_data_field)
if nargin <3
    h=heatmap(double(traj_mat_plot(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,1]);%[-0.001,0.2] for TNF
elseif nargin <4
    h=heatmap(double(traj_mat_plot(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
else
    h=heatmap(double(traj_mat_plot(:,:)),'ColorMap',color_map,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
    
end
%
XLabels = 1:size(traj_mat,2);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/1);%conf_pairs
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

if nargin ==6
    conf_pairs = conf_pairs(Outperm_col);
    
    
    CustomXLabels = string(conf_pairs);
end

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_mat,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig3 = gcf;


end


%% calculate draw correlation matrix
function [fig1,fig2,fig3] = calculate_draw_corr_mat(conf_mat,color_limit)


cn_confusion = -ones(size(conf_mat,2)); % cell number for each confusion pair

for i_cn_confusion = 1:size(cn_confusion,1)
    for j_cn_confusion = 1:size(cn_confusion,1)
        cn_confusion(i_cn_confusion, j_cn_confusion)  = sum(conf_mat(:,i_cn_confusion)&conf_mat(:,j_cn_confusion));
    end
end

jaccard_index_mat =  -ones(size(conf_mat,2));
jaccard_dist_mat = -ones(size(conf_mat,2));

for i_jaccard_index_mat = 1:size(jaccard_index_mat,1)
    for j_jaccard_index_mat = 1:size(jaccard_index_mat,1)
        if i_jaccard_index_mat == j_jaccard_index_mat
            jaccard_index_mat(i_jaccard_index_mat, j_jaccard_index_mat)  = cn_confusion(i_jaccard_index_mat, j_jaccard_index_mat) / size(conf_mat,1);
            jaccard_dist_mat(i_jaccard_index_mat, j_jaccard_index_mat) = 0;
        else
            jaccard_index_mat(i_jaccard_index_mat, j_jaccard_index_mat)  = ...
                cn_confusion(i_jaccard_index_mat, j_jaccard_index_mat)/...
                (cn_confusion(i_jaccard_index_mat, i_jaccard_index_mat) + cn_confusion(j_jaccard_index_mat, j_jaccard_index_mat) -cn_confusion(i_jaccard_index_mat, j_jaccard_index_mat));
            jaccard_dist_mat(i_jaccard_index_mat, j_jaccard_index_mat) = 1-jaccard_index_mat(i_jaccard_index_mat, j_jaccard_index_mat);
            
        end
    end
end

cn_confusion_plot = cn_confusion;
jaccard_index_mat_plot = jaccard_index_mat * 100;
jaccard_dist_mat_plot = jaccard_dist_mat;


figure(1)
paperpos=[0,0,130,130]*3;

papersize=[130,130]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

h = heatmap(cn_confusion_plot,'Colormap',parula,'GridVisible','off','MissingDataColor',[1 1 1]);
XLabels = 1:size(jaccard_index_mat,1);
CustomXLabels = string(XLabels);
CustomXLabels(:) = " ";
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
if nargin <2
    caxis([0,200]);
else
    caxis(color_limit(1,:));
end

colorbar off
fig1 = gcf;


figure(2)
paperpos=[0,0,130,130]*3;

papersize=[130,130]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

h = heatmap(jaccard_index_mat_plot,'Colormap',parula,'GridVisible','off','MissingDataColor',[1 1 1], 'CellLabelFormat', '%.1f%%');
XLabels = 1:size(jaccard_index_mat,1);
CustomXLabels = string(XLabels);
CustomXLabels(:) = " ";
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
if nargin <2
    caxis([0,100]);
else
    caxis(color_limit(2,:));
end

colorbar off
fig2 = gcf;


figure(3)
paperpos=[0,0,130,130]*3;

papersize=[130,130]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

h = heatmap(jaccard_dist_mat_plot,'Colormap',parula,'GridVisible','off','MissingDataColor',[1 1 1]);
XLabels = 1:size(jaccard_index_mat,1);
CustomXLabels = string(XLabels);
CustomXLabels(:) = " ";
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
if nargin <2
    caxis([0,1]);
else
    caxis(color_limit(3,:));
end

colorbar off
fig3 = gcf;

end

%% calculate SRS probablity vs cell number
function [prob_vec] = prob_cell_num_SRS(distinguish_mat)

maximal_cell_group_size = 10;
num_stim_pairs = 10;
N_samples = 100000;

prob_vec = -ones(1,maximal_cell_group_size);
for i_threshold = 1
    % distinguish_mat = traj_mat >= threshold(i_threshold);% confusion
    
    for i_cell_group_size = 1:maximal_cell_group_size
        
        srs_or_not = -ones(N_samples,1);
        for i_sample = 1: N_samples
            cell_ids = randperm(size(distinguish_mat,1),i_cell_group_size);
            distinguish_pairs_or_not = sum(distinguish_mat(cell_ids,:),1);
            srs_or_not(i_sample) = sum(distinguish_pairs_or_not>0) == num_stim_pairs;
        end
        
        prob_vec(i_threshold,i_cell_group_size) = sum(srs_or_not) /N_samples;
    end
    
    % plot(1:total_cell_num,prob(i_threshold,:)); hold on
    
end
end
