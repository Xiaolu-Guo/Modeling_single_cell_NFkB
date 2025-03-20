clear all

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';
savepath = fig_save_path;
addpath('./lib/')
addpath('./src/')
addpath('./bin/')

%% draw 5 ligand single ligand stim

vers = '_r1';
vers_savefig = strcat(vers,'_matching_SS_p01x');%_IkBao _IkBao_onlyPeak ,'_SS_0x' _p1x
if 1 %draw 5 ligand single ligand stim
    % Sim16_IkBao_5_signle_ligand_codon_metric _r1?
    % Sim15_5_signle_ligand_codon_metric _r2
    % Sim16_IkBao_5_signle_ligand_codon_metric_0x r1
    % Sim16_IkBao_5_signle_ligand_codon_metric_p01x
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_0x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p1x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p01x _r1
    load(strcat('../raw_data2023/Sim16_IkBao_matching_5_signle_ligand_codon_metric_p01x',vers,'.mat'))
    data_IkBo = data;
    para_mat_IkBo = sim_data_tbl.parameter_value(1:9:end,:);
    %% draw codons
    if 1
        % codon_mat = double;
        single_cell_specificity_plot = 1;
        
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
        %codon_list = {'TotalActivity'};
        
        %% draw single cells
        if single_cell_specificity_plot
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
            
            %%
            if 0
                epsilon_vec = 0.1:0.1:5;
                for i_epsilon = 1:length(epsilon_vec)
                    epsilon = epsilon_vec(i_epsilon);
                    % epsilon = 2; % Define the epsilon value
                    
                    clusterResults_wt = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_wt, 'UniformOutput', false));
                    clusterResults_IkBo = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_IkBo, 'UniformOutput', false));
                    
                    for cluster_num =1:5
                        SRS_wt(cluster_num,i_epsilon) = sum(clusterResults_wt == cluster_num) ;
                        SRS_IkBm(cluster_num,i_epsilon) = sum(clusterResults_IkBo == cluster_num) ;
                    end
                    
                end
                
                figure(1)
                paperpos = [0,0,70,50]*1.8;
                papersize = [70,50]*1.8;
                draw_pos = [10,10,50,30]*1.8;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                
                Category_percentage = SRS_wt ./ sum(SRS_wt); % Normalize to sum to 1 (or 100%)
                h = area(epsilon_vec, Category_percentage'*100); % Multiply by 100 to convert to percentages
                %                 xlabel('epsilon'); % Label for x-axis
                %                 ylabel('Percentage'); % Label for y-axis
                %                 title('Category Percentage Over Time'); % Title for the plot
                
                % Customize the colors if needed
                set(h(1), 'FaceColor', [0,0,1]); % Category 1 color
                set(h(2), 'FaceColor', [0,1,0]); % Category 2 color
                set(h(3), 'FaceColor', [1,1,0]); % Category 3 color
                set(h(4), 'FaceColor', [1,0.5,0]); % Category 4 color
                set(h(5), 'FaceColor', [1,0,0]); % Category 5 color
                
                %legend('SRS=5', 'SRS=4', 'SRS=3', 'SRS=2', 'SRS=1'); % Add legend
                
                % Ensure the y-axis goes from 0 to 100%
                ylim([0 100]);
                set(gca,'XTick',0:1:5)
                
                % Set x-axis tick labels to none
                xticklabels({});
                
                % Set y-axis tick labels to none
                yticklabels({});
                
                saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_perct_wt_cmp',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                close
                a = 1;
                
                figure(1)
                paperpos = [0,0,70,50]*1.8;
                papersize = [70,50]*1.8;
                draw_pos = [10,10,50,30]*1.8;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                
                Category_percentage = SRS_IkBm ./ sum(SRS_IkBm); % Normalize to sum to 1 (or 100%)
                h = area(epsilon_vec, Category_percentage'*100); % Multiply by 100 to convert to percentages
                %                 xlabel('epsilon'); % Label for x-axis
                %                 ylabel('Percentage'); % Label for y-axis
                %                 title('Category Percentage Over Time'); % Title for the plot
                
                % Customize the colors if needed
                set(h(1), 'FaceColor', [0,0,1]); % Category 1 color
                set(h(2), 'FaceColor', [0,1,0]); % Category 2 color
                set(h(3), 'FaceColor', [1,1,0]); % Category 3 color
                set(h(4), 'FaceColor', [1,0.5,0]); % Category 4 color
                set(h(5), 'FaceColor', [1,0,0]); % Category 5 color
                
                %legend('SRS=5', 'SRS=4', 'SRS=3', 'SRS=2', 'SRS=1'); % Add legend
                
                % Ensure the y-axis goes from 0 to 100%
                ylim([0 100]);
                set(gca,'XTick',0:1:5)
                
                % Set x-axis tick labels to none
                xticklabels({});
                
                % Set y-axis tick labels to none
                yticklabels({});
                
                saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_perct_IkBm_cmp',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                close
                a = 1;
                
            end
            
            %%
            if 0
                for epsilon = 1:3
                    % epsilon = 2; % Define the epsilon value
                    
                    clusterResults_wt = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_wt, 'UniformOutput', false));
                    clusterResults_IkBo = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_IkBo, 'UniformOutput', false));
                    
                    for cluster_num =1:5
                        barplot(cluster_num,1) = sum(clusterResults_wt == cluster_num) ;
                        barplot(cluster_num,2) = sum(clusterResults_IkBo == cluster_num) ;
                    end
                    
                    epsilon
                    barplot
                    
                    if 0
                        
                        figure(1)
                        paperpos = [0,0,40,50]*1.5;
                        papersize = [40,50]*1.5;
                        draw_pos = [10,10,20,30]*1.5;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        bar(barplot(:,1))
                        
                        
                        
                        ylim([0,800])
                        set(gca,'YTick',[0,400,800])
                        % Set x-axis tick labels to none
                        xticklabels({});
                        
                        % Set y-axis tick labels to none
                        yticklabels({});
                        
                        saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_bar_wt',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                    
                    if 0
                        
                        figure(1)
                        paperpos = [0,0,80,50]*1.5;
                        papersize = [80,50]*1.5;
                        draw_pos = [10,10,60,30]*1.5;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        bar(barplot)
                        
                        
                        
                        ylim([0,800])
                        set(gca,'YTick',[0,400,800])
                        % Set x-axis tick labels to none
                        xticklabels({});
                        
                        % Set y-axis tick labels to none
                        yticklabels({});
                        
                        saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_bar_wt_IkBm',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                    
                    a = 1;
                    if 0
                        figure(1)
                        paperpos = [0,0,70,50]*1.8;
                        papersize = [70,50]*1.8;
                        draw_pos = [10,10,50,30]*1.8;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        % Calculating the mean
                        meanValue = mean(clusterResults_wt)
                        histogram(clusterResults_wt,0.5:1:5.5)
                        set(gca,'XTick',[1,2,3,4,5])
                        
                        % Adding a line at the mean value
                        hold on; % Retain the current plot when adding new plots
                        xline(meanValue, 'r', 'LineWidth', 2); % Red line for the mean value
                        xlabel('cluster nums')
                        ylabel('counts')
                        saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_hist_wt_cmp',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                        figure(1)
                        paperpos = [0,0,70,50]*1.8;
                        papersize = [70,50]*1.8;
                        draw_pos = [10,10,50,30]*1.8;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        % Calculating the mean
                        meanValue = mean(clusterResults_IkBo)
                        histogram(clusterResults_IkBo,0.5:1:5.5)
                        set(gca,'XTick',[1,2,3,4,5])
                        
                        % Adding a line at the mean value
                        hold on; % Retain the current plot when adding new plots
                        xline(meanValue, 'r', 'LineWidth', 2); % Red line for the mean value
                        xlabel('cluster nums')
                        ylabel('counts')
                        saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_hist',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                end
                
            end
            
            
            if 0
                
                epsilon = 2; % Define the epsilon value
                
                clusterResults_wt = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_wt, 'UniformOutput', false));
                clusterResults_IkBo = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_IkBo, 'UniformOutput', false));
                
                for cluster_num =1:5
                    barplot(cluster_num,1) = sum(clusterResults_wt == cluster_num) ;
                    barplot(cluster_num,2) = sum(clusterResults_IkBo == cluster_num) ;
                end
                
                clear clusters_IkBo clusters_wt
                %                 clusters_IkBo = cell(size(codon_single_cell_IkBo,1),1);
                %                 clusters_wt = cell(size(codon_single_cell_IkBo,1),1);
                
                for i_cell = 1:size(codon_single_cell_IkBo,1)
                    
                    clusters_IkBo{i_cell} = applyDBSCAN_clusters(codon_single_cell_IkBo{i_cell},epsilon);
                    clusters_wt{i_cell} = applyDBSCAN_clusters(codon_single_cell_wt{i_cell},epsilon);
                    
                end
                ligand_id = [1,2,3,5,4];%TNF, LPS, CPG, Pam, PIC
                
                all_2combine = {[ligand_id(1),ligand_id(2)],[ligand_id(1),ligand_id(3)],...
                    [ligand_id(1),ligand_id(4)],[ligand_id(1),ligand_id(5)],...
                    [ligand_id(2),ligand_id(3)],[ligand_id(2),ligand_id(4)],...
                    [ligand_id(3),ligand_id(4)],[ligand_id(2),ligand_id(5)],...
                    [ligand_id(3),ligand_id(5)],[ligand_id(4),ligand_id(5)]};
                
                all_3combine = {[ligand_id(1),ligand_id(2),ligand_id(3)],[ligand_id(1),ligand_id(2),ligand_id(4)],...
                    [ligand_id(1),ligand_id(2),ligand_id(5)],[ligand_id(1),ligand_id(3),ligand_id(4)],...
                    [ligand_id(1),ligand_id(3),ligand_id(5)],[ligand_id(1),ligand_id(4),ligand_id(5)],...
                    [ligand_id(2),ligand_id(3),ligand_id(4)],[ligand_id(2),ligand_id(3),ligand_id(5)],...
                    [ligand_id(2),ligand_id(4),ligand_id(5)],[ligand_id(3),ligand_id(4),ligand_id(5)]};
                
                all_4combine = {[ligand_id(1),ligand_id(2),ligand_id(3),ligand_id(4)],...
                    [ligand_id(1),ligand_id(2),ligand_id(3),ligand_id(5)],...
                    [ligand_id(1),ligand_id(2),ligand_id(4),ligand_id(5)],...
                    [ligand_id(1),ligand_id(3),ligand_id(4),ligand_id(5)],...
                    [ligand_id(2),ligand_id(3),ligand_id(4),ligand_id(5)]};
                
                
                confusion_mat_IkBo_2conf = confusion_matrix_cal(clusters_IkBo,all_2combine);
                confusion_mat_IkBo_3conf = confusion_matrix_cal(clusters_IkBo,all_3combine);
                confusion_mat_IkBo_4conf = confusion_matrix_cal(clusters_IkBo,all_4combine);
                confusion_mat_IkBo_all2conf = confusion_all_matrix_cal(clusters_IkBo,all_2combine);
                
                
                confusion_mat_wt_2conf = confusion_matrix_cal(clusters_wt,all_2combine);
                confusion_mat_wt_3conf = confusion_matrix_cal(clusters_wt,all_3combine);
                confusion_mat_wt_4conf = confusion_matrix_cal(clusters_wt,all_4combine);
                confusion_mat_wt_all2conf = confusion_all_matrix_cal(clusters_wt,all_2combine);
                
                if 0 % plot wt
                    'wt'
                    total_cell = sum(clusterResults_wt == 4) + sum(clusterResults_wt == 3) + sum(clusterResults_wt == 2);
                    
                    '2conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_wt_2conf,total_cell,[]);
                    
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_wt_pairs',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '2confall, all'
                    plot_bar_conf(confusion_mat_wt_all2conf,1000,[]);
                    
                    ylim([0 1]);
                    set(gca,'YTick',0:0.5:1)%'XTick',0:1:5,
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_wt_pairs_confall',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '3conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_wt_3conf,total_cell,[]);
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_wt_triples',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '4conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_wt_4conf,total_cell,[]);
                    
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_wt_4ligands',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '2conf, 4'
                    [ax,fig] = plot_bar_conf(confusion_mat_wt_2conf,sum(clusterResults_wt == 4),clusterResults_wt == 4);
                    ylim([0 0.3]);
                    set(gca,'YTick',0:0.1:0.3)%'XTick',0:1:5,
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_wt_pairs_SRS4',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                end
                
                if 0 % plot IkBo
                    
                    'IkBo'
                    total_cell = sum(clusterResults_IkBo == 4) + sum(clusterResults_IkBo == 3) + sum(clusterResults_IkBo == 2);
                    
                    '2conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_IkBo_2conf,total_cell,[]);
                    
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_pairs',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '2confall, all'
                    plot_bar_conf(confusion_mat_IkBo_all2conf,1000,[]);
                    
                    ylim([0 1]);
                    set(gca,'YTick',0:0.5:1)%'XTick',0:1:5,
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_pairs_confall',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '3conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_IkBo_3conf,total_cell,[]);
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_triples',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '4conf, 2,3,4'
                    [ax,fig] = plot_bar_conf(confusion_mat_IkBo_4conf,total_cell,[]);
                    
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_4ligands',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    '2conf, 4'
                    [ax,fig] = plot_bar_conf(confusion_mat_IkBo_2conf,sum(clusterResults_IkBo == 4),clusterResults_IkBo == 4);
                    ylim([0 0.3]);
                    set(gca,'YTick',0:0.1:0.3)%'XTick',0:1:5,
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_pairs_SRS4',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                end
                
                if 0
                    
                    
                    conf_ligand_wt = sum(confusion_mat_wt_all2conf)/1000;
                    conf_ligand_IkBo = sum(confusion_mat_IkBo_all2conf)/1000;
                    
                    
                    
                    
                    figure(1)
                    paperpos = [0,0,100,50]*2;
                    papersize = [100,50]*2;
                    draw_pos = [10,10,90,30]*2;
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    bar([conf_ligand_wt;conf_ligand_IkBo]')
                    
                    ylim([0 0.25]);
                    set(gca,'YTick',0:0.05:0.25)%'XTick',0:1:5,
                    
                    % Set x-axis tick labels to none
                    xticklabels({});
                    
                    % Set y-axis tick labels to none
                    yticklabels({});
                    
                    
                    ylim([0 1]);
                    set(gca,'YTick',0:0.5:1)%'XTick',0:1:5,
                    saveas(gcf,strcat(fig_save_path,'SRS_confusion_IkBo_wt_pairs_confall',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                end
                
                
            end
            
            %%
            if 0
                
                epsilon = 2; % Define the epsilon value
                
                clusterResults_IkBo = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_IkBo, 'UniformOutput', false));
                clusterResults_wt = cell2mat(cellfun(@(x) applyDBSCAN(x,epsilon), codon_single_cell_wt, 'UniformOutput', false));
                
                % Preallocate the matrix for efficiency
                codon_single_cell_wt_Matrix = zeros(1000, 30);
                
                % Fill the matrix
                for i_row = 1:1000
                    % Reshape each 5x6 matrix into a 1x30 row vector and assign it to the result matrix
                    codon_single_cell_wt_Matrix(i_row, :) = reshape(codon_single_cell_wt{i_row}, 1, []);
                end
                
                
                codon_single_cell_IkBo_Matrix = zeros(1000, 30);
                
                % Fill the matrix
                for i_row = 1:1000
                    % Reshape each 5x6 matrix into a 1x30 row vector and assign it to the result matrix
                    codon_single_cell_IkBo_Matrix(i_row, :) = reshape(codon_single_cell_IkBo{i_row}, 1, []);
                end
                
                if 0
                    writematrix(codon_single_cell_wt_Matrix, 'codon_sngle_cell_wt.csv');
                    writematrix(codon_single_cell_IkBo_Matrix, 'codon_sngle_cell_SS.csv');
                    
                    
                    writematrix(clusterResults_wt, 'clusterResults_wt.csv');
                    writematrix(clusterResults_IkBo, 'clusterResults_SS.csv');
                    
                    writematrix(para_mat_wt(1:1000,:), 'para_mat_wt.csv');
                    writematrix(para_mat_IkBo(1:1000,1:size(para_mat_wt,2)), 'para_mat_IkBo.csv');
                    
                end
                
                
                if 0 % change to violin plots
                    
                    para_names = sim_data_tbl.parameter_reac_num(1,:);
                    for i_para = 1:size(para_mat_wt,2)
                        
                        % para_mat_wt_rescale = log((para_mat_wt(:,i_para)-0.01)./(1-para_mat_wt(:,i_para)));
                        for i_srs = 1:5
                            idx = (clusterResults_wt==i_srs);
                            para_vals{i_srs} = para_mat_wt(idx,i_para);%; % para_mat_wt_rescale(idx)
                        end
                        
                        para_val_plot{1} = [ para_vals{1}; para_vals{2}; para_vals{3}];
                        para_val_plot{2} = [ para_vals{4}; para_vals{5}];
                        
                        % permutation test:
                        
                        [p_val, ~, ~] = permutationTest(para_val_plot{1}, para_val_plot{2}, 10000, 'sidedness', 'both');
                        
                        
                        %
                        figure(1)
                        paperpos=[0,0,130,100]*0.8;
                        papersize=[130 100]*0.8;
                        draw_pos=[10,10,120,90]*0.8;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        % subplot(1,length(vis_data_field),i_data_field)
                        al_goodplot_cellinput(para_val_plot,[],0.5,ones(length(para_val_plot),1)*[ 0 0 0] ,'left',[],std(para_mat_wt(:,i_para))/2500);
                        al_goodplot_cellinput(para_val_plot,[],0.5,ones(length(para_val_plot),1)*[0 0 0],'right',[],std(para_mat_wt(:,i_para))/2500);
                        % Unilateral plots for 2 timepoints (left: before, right: after), 3 groups.
                        % One can produce multiple plots at once using a P cell, Nx1 mat for each cell input, P plots (1 per column).
                        % One can use different options for each column.
                        % If options are given only for 1 column, it is replicated for the others.
                        grid off
                        % Assuming you want to display the p-value between the first and second violin plots
                        % Calculate position for the p-value text
                        % Format the p-value as a string in scientific notation
                        pval_str = sprintf('p = %.2e', p_val);
                        pos = 1:2;
                        pValueX = mean(pos(1:2)); % Midpoint between the first two plots
                        maxYValue = max([max(para_val_plot{1}), max(para_val_plot{2})]); % Find the maximum Y value among the first two datasets
                        YL = ylim();
                        pValueY = (YL(2)-YL(1))*0.9 + YL(1);% + diff(ylim) * 0.05; % Position the text slightly above the max value
                        
                        % Display the p-value
                        text(pValueX, pValueY, pval_str, 'HorizontalAlignment', 'center', 'FontSize', 8);
                        
                        % Set x-axis tick labels to none
                        xticklabels({});
                        
                        % Set y-axis tick labels to none
                        yticklabels({});
                        %
                        xlim([0,3]);
                        % xticklabels({'\leq 3 ','>3'})
                        
                        
                        % set(gca,'fontsize',14,'fontname','Arial');
                        saveas(gcf,strcat(fig_save_path,'Para_violin_distrib_by_SRS_p',num2str(para_names(i_para))),'epsc');
                        close
                        
                    end
                    
                    
                end
                
                
                if 0
                    for i_para = 1:size(para_mat_wt,2)
                        for i_srs = 1:5
                            idx = (clusterResults_wt==i_srs);
                            mean_vals(i_para,i_srs) = mean(para_mat_wt(idx,i_para));
                            std_vals(i_para,i_srs) = std(para_mat_wt(idx,i_para));
                            
                            mean_vals_IkBo(i_para,i_srs) = mean(para_mat_IkBo(idx,i_para));
                            std_vals_IkBo(i_para,i_srs) = std(para_mat_IkBo(idx,i_para));
                        end
                        paperpos = [0,0,70,50]*1.8;
                        papersize = [70,50]*1.8;
                        draw_pos = [10,10,50,30]*1.8;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        % Plotting with error bars
                        errorbar(1:5, mean_vals(i_para, :), std_vals(i_para, :), 'LineWidth', 1.5); hold on;
                        % errorbar(1:5, mean_vals_IkBo(i_para, :), std_vals_IkBo(i_para, :), 'LineWidth', 1.5); hold on;
                        
                        %xlabel('SRS Label');
                        %ylabel('Parameter Value');
                        %title('Parameter Values with Error Bars');
                        % legend('WT', 'IkBo');
                        set(gca,'XTick',[1,2,3,4,5])
                        
                        saveas(gcf,strcat(fig_save_path,'Single_cell_wt_SRS_paraid_',num2str(i_para),'_',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        
                        % saveas(gcf,strcat(fig_save_path,'Single_cell_SRS_paraid_',num2str(i_para),'_',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                    end
                    
                    %plot(1:5, mean_vals(1,:)); hold on
                    % plot(1:5, mean_vals_IkBo(1,:)); hold on
                    
                    
                end
                
                a = 1;
                
                
                %                 dataMatrix = codon_single_cell_wt;
                %
                %                 % Assuming you have a function umap to apply UMAP reduction
                %                 [reducedData, umap] = run_umap(dataMatrix);
                %
                %                 % Assuming clusterLabels contains cluster IDs for each row in dataMatrix
                %                 specificCluster = reducedData(clusterLabels == 1, :);
                %
                %                 % Visualization
                %                 figure;
                %                 scatter(specificCluster(:,1), specificCluster(:,2), 10, 'filled');
                %                 title('Visualization of Cluster 1 using UMAP');
                %                 xlabel('UMAP 1');
                %                 ylabel('UMAP 2');
                %
                %                 a = 1;
            end
            
        end
    end
    
    %% draw heat map IkBo-/-
    if 1
        clear data_NFkB
        cells_inteval = 50;
        
        data_NFkB.info_ligand = data_IkBo.info_ligand;
        data_NFkB.info_dose_str = data_IkBo.info_dose_str;
        data_NFkB.info_dose_index = data_IkBo.info_dose_index;
        for i_ligand = 1:length(data_IkBo.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data_IkBo.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_plot = [];
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_plot = [data_plot,data_NFkB.model_sim{i_ligand}(1:cells_inteval:end,:)];
            %data_NFkB.exp = data_NFkB.model_sim;
        end
        
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        [~,order_to_plot] = sort(metrics{1}.(metric_order_name_{1})(1:cells_inteval:end,metric_order_column),'descend');
        
        figure(2)
        
        clear traj_mat;
        traj_mat = data_plot;
        
        % Step 1: Perform hierarchical clustering on the rows
        Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
        Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
        % Step 2: Determine the order of rows based on hierarchical clustering
        [H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
        close(gcf); % Close dendrogram figure
        
        % Step 3: Reorder the matrix based on the clustering result
        reordered_traj_mat = traj_mat(Outperm, :);
        
        % If you need to display the dendrogram alongside, you can plot it separately
        figure(2)
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
        saveas(gcf,strcat(savepath,'SRS_heatmap_IkBm_hierarchical_cluster_20cells'),'epsc');%_20cells
        close
        
        
        figure(1)
        
        paperpos=[0,0,246,64]*3;
        papersize=[246,64]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(reordered_traj_mat(order_to_plot,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        %
        XLabels = 0:5:((size(reordered_traj_mat,2)-1)*5);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/60);
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:size(data_plot,1);
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
        saveas(gcf,strcat(savepath,'SRS_heatmap_IkBm_20cells'),'epsc');%_20cells
        close
        
    end
    
    if 0
        data_NFkB.info_ligand = data_IkBo.info_ligand;
        data_NFkB.info_dose_str = data_IkBo.info_dose_str;
        data_NFkB.info_dose_index = data_IkBo.info_dose_index;
        for i_ligand = 1:length(data_IkBo.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data_IkBo.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        
        order_name_vec = {'matched_sc_srs_TNF','matched_sc_srs_LPS','matched_sc_srs_CpG','matched_sc_srs_PolyIC','matched_sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers_savefig);
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 1 % 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{i_order_i_ligand}.(metric_order_name_{1})(:,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
    %% draw heat map wt
    if 1
        clear data_NFkB
        cells_inteval = 50;
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        data_plot = [];
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_plot = [data_plot,data_NFkB.model_sim{i_ligand}(1:cells_inteval:end,:)];
            %data_NFkB.exp = data_NFkB.model_sim;
        end
        
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        [~,order_to_plot] = sort(metrics{length_metrics+1}.(metric_order_name_{1})(1:cells_inteval:end,metric_order_column),'descend');
        
        figure(2)
        clear traj_mat;
        traj_mat = data_plot;
        
        % Step 1: Perform hierarchical clustering on the rows
        Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
        Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
        % Step 2: Determine the order of rows based on hierarchical clustering
        [H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
        close(gcf); % Close dendrogram figure
        
        % Step 3: Reorder the matrix based on the clustering result
        reordered_traj_mat = traj_mat(Outperm, :);
        
        % If you need to display the dendrogram alongside, you can plot it separately
        figure(2)
        paperpos=[0,0,55,70];
        papersize=[55 70];
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)

        [dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
        set(gca, 'XDir', 'reverse','YDir','reverse');
        set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
        
        xticklabels('')
        yticklabels('')
        axis off
        saveas(gcf,strcat(savepath,'SRS_heatmap_wt_hierarchical_cluster_20cells'),'epsc');%_20cells
        close
     
        
        figure(1)
        
        paperpos=[0,0,246,64]*3;
        papersize=[246,64]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
                
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(data_plot(order_to_plot,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        %
        XLabels = 0:5:((size(data_plot,2)-1)*5);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/60);
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:size(data_plot,1);
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
        saveas(gcf,strcat(savepath,'SRS_heatmap_wt_20cells'),'epsc');%_20cells
        close
        
    end
    
    if 0
        clear data_NFkB
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        
        vers_savefig = '_matching_wt';
        order_name_vec = {'matched_sc_srs_TNF','matched_sc_srs_LPS','matched_sc_srs_CpG','matched_sc_srs_PolyIC','matched_sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers_savefig);
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 1% 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{length_metrics+i_order_i_ligand}.(metric_order_name_{1})(:,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
end

%%
function cluster_num = applyDBSCAN(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end

%%
function cluster_results = applyDBSCAN_clusters(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
i_cluster = 1;
idx_rest = idx';
for i_cluster = 1:max(idx)
    cluster_results{i_cluster} = find(idx_rest == i_cluster);
end
end

function confusion_mat = confusion_matrix_cal(clusters_genetype,all_combine)
confusion_mat = zeros(length(clusters_genetype),length(all_combine));
for i_cell = 1:length(clusters_genetype)
    for i_2combine = 1:length(all_combine)
        for i_clusters = 1:length(clusters_genetype{i_cell})
            if isequal(sort(all_combine{i_2combine}), sort(clusters_genetype{i_cell}{i_clusters}))
                confusion_mat(i_cell,i_2combine) = 1;
            end
            
        end
        
    end
    
end
end

function confusion_mat = confusion_all_matrix_cal(clusters_genetype,all_combine)
confusion_mat = zeros(length(clusters_genetype),length(all_combine));
for i_cell = 1:length(clusters_genetype)
    for i_2combine = 1:length(all_combine)
        for i_clusters = 1:length(clusters_genetype{i_cell})
            if all(ismember(all_combine{i_2combine}, clusters_genetype{i_cell}{i_clusters}))
                confusion_mat(i_cell,i_2combine) = 1;
            end
        end
    end
end

end

function [ax,fig] = plot_bar_conf(confusion_mat,total_cell,idx_cell)

if isempty(total_cell)
    total_cell = size(confusion_mat,1);
end

if isempty(idx_cell)
    conf_ligand = sum(confusion_mat)/total_cell;
    
else
    conf_ligand = sum(confusion_mat(idx_cell,:))/total_cell;
end

figure(1)
paperpos = [0,0,100,50]*2;
papersize = [100,50]*2;
draw_pos = [10,10,90,30]*2;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
conf_ligand
bar(conf_ligand)

ylim([0 0.25]);
set(gca,'YTick',0:0.05:0.25)%'XTick',0:1:5,

% Set x-axis tick labels to none
xticklabels({});

% Set y-axis tick labels to none
yticklabels({});
ax = gca;
fig = gcf;
end