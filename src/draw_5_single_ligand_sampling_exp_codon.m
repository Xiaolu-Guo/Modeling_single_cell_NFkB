%%
% compare sim vs exp, and scramble sim vs exp
clear all
addpath('./exercise/')
if 1
    data_save_file_path = '../raw_data2023/';%_fay_parameter/';
    
    
    %% load data
    % sim, metrics can be used
    load('../raw_data2023/Sim15_5_signle_ligand_codon_metric.mat')
    metrics_sim = metrics;
    data_sim = data;
    
%     load('../raw_data2023/Sim8_5_signle_ligand_codon_metric_r2.mat')
%     metrics_sim = metrics;
%     data_sim = data;
%     
    % scrambled sim, metrics can be used
    load('../raw_data2023/Sim14_5_signle_ligand_codon_metric_r2.mat')
    metrics_scrambled_sim = metrics;
    data_scrambled_sim = data;
    
    % exp,metrics can be used
    load('Supriya_data_metrics.mat')
    load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))
    metrics_exp = metrics;
    data_exp = data;
    
    index_exp_ligand_vects = {[3];%6
        [6];%7
        [12];%8
        [16];%9
        [19]};%31
    
    %% cal codon and sim vs exp diff
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    i_index = 1;
    clear metrics_cal data_info
    field_list = fieldnames(metrics_exp{1});
    for i_cond = 1:length(index_exp_ligand_vects)
        index_exp_ligand = index_exp_ligand_vects{i_cond};
        sti_vec{i_cond} = data_sim.info_ligand(i_cond);
        
        for i_rpc = 1:length(index_exp_ligand)
            for i_field = 1:length(field_list)
                metrics_cal{i_index}.(field_list{i_field}) = metrics_exp{index_exp_ligand(i_rpc)}.(field_list{i_field});
                metrics_cal{i_index+1}.(field_list{i_field}) = metrics_sim{i_cond}.(field_list{i_field})(1:9:end,:);
                metrics_cal{i_index+2}.(field_list{i_field}) = metrics_scrambled_sim{i_cond}.(field_list{i_field})(1:9:end,:);
                
            end
            data_info.info_ligand(i_index:i_index+2) =[data_exp.info_ligand(index_exp_ligand(i_rpc)),...
                data_sim.info_ligand(i_cond),...
                data_scrambled_sim.info_ligand(i_cond)];
            
            data_info.info_dose_str(i_index:i_index+2) =[data_exp.info_dose_str(index_exp_ligand(i_rpc)),...
                data_sim.info_dose_str(i_cond),...
                data_scrambled_sim.info_dose_str(i_cond)];
            
            data_info.data_label(i_index:i_index+2) = {'exp','sim','scramble_sim'};
            
            i_index = i_index+3;
            
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
        
        %% exp vs sim
        W_diff_sim(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
            collect_feature_vects_all.(codon_list{i_codon}){i_index+1});
        
        
        %% exp vs scramble sim
        W_diff_scrambled_sim(i_codon,i_cond) = ...
            cal_w_dist(collect_feature_vects_all.(codon_list{i_codon}){i_index},...
            collect_feature_vects_all.(codon_list{i_codon}){i_index+2});
        
        %% minmax rescaled exp vs sim
        W_diff_sim_minmax(i_codon,i_cond) = cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
            collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+1});
        
        %% exp vs scramble sim
        W_diff_scrambled_sim_minmax(i_codon,i_cond) = ...
            cal_w_dist(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index},...
            collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){i_index+2});
    end
    i_index = i_index+3;
end

W_diff_scrambled_sim_mat = zeros(6,1);
W_diff_sim_mat = zeros(6,1);
W_diff_scrambled_sim_mat= mean(W_diff_scrambled_sim,2);
W_diff_sim_mat = mean(W_diff_sim,2);
W_diff_sim_minmax_mat = mean(W_diff_sim_minmax,2);
W_diff_scrambled_sim_minmax_mat= mean(W_diff_scrambled_sim_minmax,2);

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
        sti_vec{i_cond} = data_sim.info_ligand(i_cond);
    end
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,W_diff_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand_weighted_match_wdis_sim_vs_exp'),'epsc')
    close()
    
    figure(1)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,W_diff_scrambled_sim_minmax,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_5_single_ligand__wdis_scramble_sim_vs_exp'),'epsc')
    close()
end

if 0
    fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
    fig_opt.paper_opt.papersize=[220 180]*3;
    fig_opt.save_file = strcat(fig_save_path,'codon_5_ligand_2023_12_distrib');
    violin_plot_codon_2023_12(collect_feature_vects_all,fig_opt) % draw all codon and sti in one
    
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
