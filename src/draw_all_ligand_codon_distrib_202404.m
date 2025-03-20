function [] = draw_all_ligand_codon_distrib_202404(data_save_file_path,fig_save_path)
% For Guo et al. Figure 2B-C, sim vs exp codon Wasserstein dist
% 
% tested 05/10/2024, Matlab 2020a

fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
fig_opt.paper_opt.papersize=[220 180]*3;

% data_proj_num_vec = data_proj_nums{i_ligand};
%     load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))


%% codon distribution
% violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
% codon and sti seperately, and save
fig_opt.save_file = strcat(fig_save_path,'All_ligand_codon_2023_06_distrib');

codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
sti_vec = cell(1,15);
ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
dose_label = {'L','M','H'};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

good_fit_pop = NaN(5,3);

for i_codon = 1:length(codon_list)
    i_sti_exp = 1;
    
    for i_ligand_exp = 1:length(ligand_vec)%TNF,LPS [27];%
        
        for i_dose_exp = 1:length(dose_vec{i_ligand_exp})
            
            i_data_exp = find(categorical(data.info_ligand)==ligand_vec{i_ligand_exp} & categorical(data_dose_str)==dose_vec{i_ligand_exp}{i_dose_exp});
            sti_vec{i_sti_exp} = strcat(data.info_ligand{i_data_exp},'-',dose_label{i_dose_exp});
            exp_data{i_codon,i_sti_exp} = collect_feature_vects.(codon_list{i_codon}){i_data_exp*2-1};
            sim_data{i_codon,i_sti_exp} = collect_feature_vects.(codon_list{i_codon}){i_data_exp*2};
            i_sti_exp = i_sti_exp +1;
        end
    end
end


w_dis = [];

for i_codon = 1:length(codon_list)
    for i_sti_exp = 1:size(exp_data,2)
        for i_sti_sim = 1:size(sim_data,2)
            
            if i_sti_exp == i_sti_sim
                W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                    w_distance(exp_data{i_codon,i_sti_exp}, sim_data{i_codon,i_sti_sim}, 2);
            elseif i_sti_exp <i_sti_sim
                W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                    w_distance(sim_data{i_codon,i_sti_exp}, sim_data{i_codon,i_sti_sim}, 2);
            else
                W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                    w_distance(exp_data{i_codon,i_sti_exp}, exp_data{i_codon,i_sti_sim}, 2);
            end
            
        end
    end
end

W_diff_mat = mean(W_diff,3);

index_reorder = [1,2,3,13,14,15,7,8,9,4,5,6,10,11,12];

W_diff_mat = W_diff_mat(index_reorder,index_reorder);

%% Figure 2C: tested 05/10/2024, Matlab 2020a
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

if 1
    figure (1)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,450]/2;
    paper_size = [450,450]/2;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    
    heatmap(W_diff_mat,'colormap',mymap,'ColorbarVisible','off')
    caxis([0,0.2])
    
    % Getting handle of the heatmap
    h = gca;
    
    % Removing tick labels from x and y axes
    h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
    h.YDisplayLabels = repmat({''}, size(h.YDisplayData));
    
    saveas(gcf,strcat(fig_save_path,'Wdis_exp_sim_202401'),'epsc')
    close()
    
end

%% Figure 2B: tested 05/10/2024, Matlab 2020a

if 1
    i_exp_sim = 1;
    i_sim = 1;
    i_exp = 1;
    for i_sti_exp = 1:size(exp_data,2)
        for i_sti_sim = 1:size(sim_data,2)
            
            if i_sti_exp == i_sti_sim
                W_diff_exp_sim(i_exp_sim) = W_diff_mat(i_sti_exp,i_sti_sim);
                i_exp_sim = i_exp_sim+1;
            elseif i_sti_exp <i_sti_sim
                W_diff_sim(i_sim) = W_diff_mat(i_sti_exp,i_sti_sim);
                i_sim = i_sim+1;
            else
                W_diff_exp(i_exp) = W_diff_mat(i_sti_exp,i_sti_sim);
                i_exp = i_exp+1;
            end
            
        end
    end
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,250,150];
    paper_size = [250,150];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    histogram(W_diff_exp_sim,'Normalization','probability'); hold on
    histogram(W_diff_sim,'Normalization','probability'); hold on
    histogram(W_diff_exp,'Normalization','probability'); hold on
    
    % Getting handle of the current axes
    h = gca;
    
    % Removing tick labels from x and y axes
    h.XTickLabel = [];
    h.YTickLabel = [];
    
    saveas(gcf,strcat(fig_save_path,'Wdis_exp_sim_hist_202401'),'epsc')
    
    close()
    
end




