clear all

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

%% draw 5 ligand single ligand stim

vers = '_r1';
vers_savefig = strcat(vers,'_SS_0x');%_IkBao _IkBao_onlyPeak
if 1 %draw 5 ligand single ligand stim
    % Sim16_IkBao_5_signle_ligand_codon_metric _r1?
    % Sim15_5_signle_ligand_codon_metric _r2
    % Sim16_IkBao_5_signle_ligand_codon_metric_0x r1
    % Sim16_IkBao_5_signle_ligand_codon_metric_p01x
    load(strcat('../raw_data2023/Sim16_IkBao_5_signle_ligand_codon_metric_0x',vers,'.mat'))
    
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
        
        load(strcat('../raw_data2023/Sim15_5_signle_ligand_codon_metric_r2','.mat'))
        
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
            codon_single_cell_SS=cell(2,1);

            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:length_data
                    for i_codon = 1:length(codon_list)
                        codon_single_cell_SS{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
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
            
            clusterResults_wt = cell2mat(cellfun(@applyDBSCAN, codon_single_cell_wt, 'UniformOutput', false));
            clusterResults_SS = cell2mat(cellfun(@applyDBSCAN, codon_single_cell_SS, 'UniformOutput', false));
            
            
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
            saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_hist_wt_cmp_0x_eps2'),'epsc')
            close
            
            figure(1)
            paperpos = [0,0,70,50]*1.8;
            papersize = [70,50]*1.8;
            draw_pos = [10,10,50,30]*1.8;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            % Calculating the mean
            meanValue = mean(clusterResults_SS)
            histogram(clusterResults_SS,0.5:1:5.5)
            set(gca,'XTick',[1,2,3,4,5])
            
            % Adding a line at the mean value
            hold on; % Retain the current plot when adding new plots
            xline(meanValue, 'r', 'LineWidth', 2); % Red line for the mean value
            xlabel('cluster nums')
            ylabel('counts')
            saveas(gcf,strcat(fig_save_path,'Single_cell_Cluster_hist_SS_0x_eps2'),'epsc')
            close

        end
    end
    
    %% draw heat map
    if 0
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
        
        order_name_vec = {'sc_srs_TNF','sc_srs_LPS','sc_srs_CpG','sc_srs_PolyIC','sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers_savefig);
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 1% 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{i_order_i_ligand}.(metric_order_name_{1})(:,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
end

%%
function cluster_num = applyDBSCAN(matrix)

epsilon = 2; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end