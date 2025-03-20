%%
% compare sim vs exp, and scramble sim vs exp
clear all
addpath('./exercise/')
vers_vec = {'_r3'};%,'_r4','_r5','_r6'
for i_vers = 1 :length(vers_vec)
    vers = vers_vec{i_vers};
    if 1
        data_save_file_path = '../raw_data2023/';%_fay_parameter/';
        
        %% load data
        % sim, metrics can be used
        load(strcat('../raw_data2023/Sim15_5_signle_ligand_codon_metric',vers,'.mat'))
        metrics_sim_weighted_match = metrics;
        data_sim_weighted_match = data;
        
        load(strcat('../raw_data2023/Sim8_5_signle_ligand_codon_metric',vers,'.mat'))
        metrics_sim = metrics;
        data_sim = data;
        %
        % scrambled sim, metrics can be used
        load(strcat('../raw_data2023/Sim14_5_signle_ligand_codon_metric',vers,'.mat'))
        metrics_scrambled_sim = metrics;
        data_scrambled_sim = data;
        
        % totally sampled sim, sample from the distribution
        load(strcat('../raw_data2023/Sim2',vers,'_srs_codon_metric_2024.mat'))
        metrics_sample_dist_sim = metrics;
        data_sample_dist_sim = data;
        
        % exp,metrics can be used
        % load('Supriya_data_metrics.mat')
        load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat')) % Ade's dataset
        metrics_exp = metrics;
        data_exp = data;
        
        index_exp_ligand_vects = {[3];%6
            [6];%7
            [11];%8
            [16];%9
            [18]};%31
        
        %% cal codon and sim vs exp diff
        codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
        i_index = 1;
        clear metrics_cal data_info
        field_list = fieldnames(metrics_exp{1});
        for i_cond = 1:length(index_exp_ligand_vects)
            index_exp_ligand = index_exp_ligand_vects{i_cond};
            sti_vec{i_cond} = data_sim_weighted_match.info_ligand(i_cond);
            
            for i_rpc = 1:length(index_exp_ligand)
                
                for i_field = 1:length(field_list)
                    metrics_cal{i_index}.(field_list{i_field}) = metrics_exp{index_exp_ligand(i_rpc)*2-1}.(field_list{i_field});
                    metrics_cal{i_index+1}.(field_list{i_field}) = metrics_sim_weighted_match{i_cond}.(field_list{i_field})(1:9:end,:);
                    metrics_cal{i_index+2}.(field_list{i_field}) = metrics_scrambled_sim{i_cond}.(field_list{i_field})(1:9:end,:);
                    metrics_cal{i_index+3}.(field_list{i_field}) = metrics_sim{i_cond}.(field_list{i_field})(1:9:end,:);
                    metrics_cal{i_index+4}.(field_list{i_field}) = metrics_sample_dist_sim{i_cond}.(field_list{i_field})(1:9:end,:);
                end
                
                data_info.info_ligand(i_index:i_index+4) =[data_exp.info_ligand(index_exp_ligand(i_rpc)),...
                    data_sim_weighted_match.info_ligand(i_cond),...
                    data_scrambled_sim.info_ligand(i_cond),...
                    data_sim.info_ligand(i_cond),...
                    data_sample_dist_sim.info_ligand(i_cond)];
                
                data_info.info_dose_str(i_index:i_index+4) =[data_exp.info_dose_str(index_exp_ligand(i_rpc)),...
                    data_sim_weighted_match.info_dose_str(i_cond),...
                    data_scrambled_sim.info_dose_str(i_cond),...
                    data_sim.info_dose_str(i_cond),...
                    data_sample_dist_sim.info_dose_str(i_cond)];
                
                data_info.data_label(i_index:i_index+4) = {'exp','sim_weighted_match','scramble_sim','sim','sample_dist_sim'};
                
                i_index = i_index+5;
                
            end
            % good_fit_pop_wdis(i_ligand,i_dose) = mean(W_diff{i_cond}(:,i_sti));
            
        end
        
        [collect_feature_vects_all_minmax_scaled,metrics_all_minmax_scaled] = calculate_codon_from_metric2023(data_info,metrics_cal); %,  parameter
        [collect_feature_vects_all,metrics_all] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %,  parameter
        save('codons_5_single_ligand_sampling_exp.mat','collect_feature_vects_all_minmax_scaled','metrics_all_minmax_scaled','collect_feature_vects_all','metrics_all')
    end
    i_index = 1;
    for i_cond = 1:length(index_exp_ligand_vects)
        
        for i_codon = 1:length(codon_list)
            
            %% exp vs weighted match sim
            W_diff_weighted_match_sim(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all.(codon_list{i_codon}){i_index+1});
            
            %% exp vs scramble sim
            W_diff_scrambled_sim(i_codon,i_cond) = ...
                cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all.(codon_list{i_codon}){i_index+2});
            
            %% exp vs sim
            W_diff_sim(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all.(codon_list{i_codon}){i_index+3});
            
            
            %% exp vs sample dist sim
            W_diff_sample_dist_sim(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all.(codon_list{i_codon}){i_index+4});
            
            %% minmax rescaled exp vs weighted match sim
            W_diff_weighted_match_sim_minmax(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+1});
            
            %% exp vs scramble sim
            W_diff_scrambled_sim_minmax(i_codon,i_cond) = ...
                cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+2});
            
            %% exp vs sim
            W_diff_sim_minmax(i_codon,i_cond) = ...
                cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+3});
            
            %% exp vs sample dist sim
            W_diff_sample_dist_sim_minmax(i_codon,i_cond) = ...
                cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
                collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+4});
            
        end
        i_index = i_index+5;
    end
    
    % scrambled
    W_diff_scrambled_sim_mat = zeros(6,1);
    W_diff_scrambled_sim_mat= mean(W_diff_scrambled_sim,2);
    W_diff_scrambled_sim_minmax_mat= mean(W_diff_scrambled_sim_minmax);
    W_diff_scrambled_sim_minmax_mat_codon= mean(W_diff_scrambled_sim_minmax,2);
    
    % sim
    W_diff_sim_mat = zeros(6,1);
    W_diff_sim_mat = mean(W_diff_sim,2);
    W_diff_sim_minmax_mat = mean(W_diff_sim_minmax);%% plot hist?
    W_diff_sim_minmax_mat_codon = mean(W_diff_sim_minmax,2);%% plot hist?
    
    % weighted match
    W_diff_weighted_match_sim_mat = zeros(6,1);
    W_diff_weighted_match_sim_mat = mean(W_diff_weighted_match_sim,2);
    W_diff_weighted_match_sim_minmax_mat = mean(W_diff_weighted_match_sim_minmax);%% plot hist?
    W_diff_weighted_match_sim_minmax_mat_codon = mean(W_diff_weighted_match_sim_minmax,2);%% plot hist?
    
    % sample dist sim
    W_diff_sample_dist_sim_mat = zeros(6,1);
    W_diff_sample_dist_sim_mat = mean(W_diff_sample_dist_sim,2);
    W_diff_sample_dist_sim_minmax_mat = mean(W_diff_sample_dist_sim_minmax);%% plot hist?
    W_diff_sample_dist_sim_minmax_mat_codon = mean(W_diff_sample_dist_sim_minmax,2);%% plot hist?
    
    W_diff_comp = [W_diff_scrambled_sim_minmax_mat;W_diff_sim_minmax_mat;W_diff_weighted_match_sim_minmax_mat;W_diff_sample_dist_sim_minmax_mat]';
    W_diff_comp_codon = [W_diff_scrambled_sim_minmax_mat_codon,W_diff_sim_minmax_mat_codon,W_diff_weighted_match_sim_minmax_mat_codon,W_diff_sample_dist_sim_minmax_mat_codon];
    
    Xlabel = categorical({'TNF','LPS','CpG','PolyIC','Pam3CSK'});
    Xlabel = reordercats(Xlabel,{'TNF','LPS','CpG','PolyIC','Pam3CSK'});
    
    Xlabel_codon = categorical({'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'});
    Xlabel_codon = reordercats(Xlabel_codon,{'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'});
    
    %
    % index = [1:8,10,14,16];
    % compare_sim = W_diff_scrambled_sim_mat - W_diff_sim_mat;
    % sum(sign(compare_sim))
    % sum(sign(compare_sim(:,index)),2)
    %
    % scramble_sim_good_fit = W_diff_scrambled_sim_mat<=1;
    % sim_good_fit = W_diff_sim_mat<=1;
    %
    % index = [1:8,10,14,16];
    % compare_sim_minmax = W_diff_scrambled_sim_minmax_mat - W_diff_sim_minmax_mat;
    % sum(sign(compare_sim_minmax))
    % sum(sign(compare_sim_minmax(:,index)),2)
    %
    % scramble_sim_good_fit_minmax = W_diff_scrambled_sim_minmax_mat<=0.5;
    % sim_good_fit_minmax = W_diff_sim_minmax_mat<=0.5;
    
    %% plot
    if 1 % W-distance, run this for figure
        clear sti_vec
        fig_save_path = '../SubFigures2023/';
        for i_cond = 1:5
            sti_vec{i_cond} = data_sim_weighted_match.info_ligand(i_cond);
        end
        mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
            ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
        
        paper_pos = [0,0,250,150];
        paper_size = [250,150];
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,W_diff_sample_dist_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
        caxis([0,0.4])
        
        saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand_wdis_sample_dist_sim_vs_exp',vers),'epsc')
        close()
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,W_diff_weighted_match_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
        caxis([0,0.4])
        
        saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand_wdis_weighted_match_sim_vs_exp',vers),'epsc')
        close()
        
        figure(1)
        set(gcf, 'PaperUnits','points')
        
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,W_diff_scrambled_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
        caxis([0,0.4])
        
        saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand_wdis_scramble_sim_vs_exp',vers),'epsc')
        close()
        
        figure(3)
        set(gcf, 'PaperUnits','points')
        
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,W_diff_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
        caxis([0,0.4])
        
        saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand_wdis_sim_vs_exp',vers),'epsc')
        close()
        
        
                paper_pos = [0,0,350,150];
        paper_size = [350,150];
        
        figure(4)
        set(gcf, 'PaperUnits','points')
        
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        bar(Xlabel,W_diff_comp)
        lgd = legend({'Scramble','Matched','Weighted Matched','Sampling dist'},'Location','NorthEast');
        set(lgd,'Box','off')
        ylabel('W-distance')
        ylim([0,0.3])
        saveas(gcf,strcat(fig_save_path,'Avg_Wdis_scramble_sim_weighted_vs_exp_by_ligand',vers),'epsc')
        close()
        
        figure(5)
        set(gcf, 'PaperUnits','points')
        
        
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        bar(Xlabel_codon,W_diff_comp_codon)
        lgd = legend({'Scramble','Matched','Weighted Matched','Sampling dist'},'Location','NorthEast');
        set(lgd,'Box','off')
        ylabel('W-distance')
        ylim([0,0.3])
        saveas(gcf,strcat(fig_save_path,'Avg_Wdis_scramble_sim_weighted_vs_exp_by_codon',vers),'epsc')
        close()
        
    end
    
    if 0
        fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
        fig_opt.paper_opt.papersize=[220 180]*3;
        fig_opt.save_file = strcat(fig_save_path,'codon_5_ligand_2023_12_distrib',vers);
        violin_plot_codon_2023_12(collect_feature_vects_all,fig_opt) % draw all codon and sti in one
        
    end
end

%% sub functions
function [w_dis] = cal_w_dist(exp_data,sim_data)

pts = linspace(min(min(exp_data),min(sim_data)),max(max(exp_data),max(sim_data)),50);

[~,~,bw_exp] = ksdensity(exp_data,pts,...
    'Function','pdf');
[~,~,bw_sim] = ksdensity(sim_data,pts,...
    'Function','pdf');%
bw = min(bw_exp,bw_sim);

[f_exp,xi_exp] = ksdensity(exp_data,pts,...
    'Function','pdf','Bandwidth',bw);%,'Bandwidth',bw
[f_sim,xi_sim] = ksdensity(sim_data,pts,...
    'Function','pdf','Bandwidth',bw);%

w_dis = w_distance(exp_data, sim_data, 2);
end
