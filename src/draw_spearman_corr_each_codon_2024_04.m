function [] = draw_spearman_corr_each_codon_2024_04(data_save_file_path,fig_save_path)
% For Guo et al. Figure 3B, spearman corr scatter plots
% 
% tested 05/12/2024, Matlab 2020a

parper_size_codon_para = [100,100];
parper_size_para_codon = [100,80];
std_thresh = 0.1;
% load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'};
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'10ug/mL';'33ug/mL';'100ug/mL'}};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

% fig_save_path = strcat(fig_save_path,'2023_03_');
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];


collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};

para_names = {'rcp syn','C deg','end','TAK ac', 'time d1','time d2','Tot NFkB'};
color_vec =  [[119 180 202]/255;% TNF_color =
    [229 129 56]/255;% P3CSK_color =
    [137 180 66]/255; % CpG_color =
    [222 78 66]/255; % LPS_color =
    [101 77 123]/255]; % PolyIC_cclor =


if 1
    corr_mat_para = NaN(6,5,7);
    feature_corr_mat_para =  NaN(6,5,7);
    data_color = double(0);
    
    i_data =1;
    
    for i_feature = 1:length(features)
        corr_mat_codon = NaN(7,5);
        pval_mat_codon = NaN(7,5);
        i_cond = 1;
        for i_sti = 1:length(ligand_vec)
            ligand = ligand_vec{i_sti};
            % dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
            for i_dose = length(dose_vec{i_sti})
                dose  = dose_vec{i_sti}{i_dose};
                para_vals = data.parameters_mode_nan{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose};
                para_vals = para_vals(:,1:end-1);
                codon_corr = collect_feature_vects.(features{i_feature}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'};
                sti_vec{i_cond} = strcat(ligand,'-High');
                sti_vec2{i_cond} = ligand;
                
                %% core para
                para_corr = para_vals(:,[1:3,size(para_vals,2)]);
                [corr_mat,core_pval] =  corr(para_corr,codon_corr,'Type','Spearman');
                
                para_names_str = data.para_name{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose}(:,1:end-1);
                para_names_str{end} = 'NFKBInit';
                
                corr_mat_codon(4:end,i_cond) = corr_mat;
                pval_mat_codon(4:end,i_cond) = core_pval;
                
                
                %% receptor para
                paras_corr= para_vals(:,setdiff(1:size(para_vals,2),[1:3,size(para_vals,2)],'stable'));
                codon_corr = collect_feature_vects.(features{i_feature}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'};
                [rcp_mat,rcp_pval] =  corr(paras_corr,codon_corr,'Type','Spearman');
                
                corr_mat_codon(1:length(rcp_mat),i_cond) = rcp_mat;
                pval_mat_codon(1:length(rcp_mat),i_cond) = rcp_pval;
                
                i_cond = i_cond+1;
            end
        end
        
        corr_rel = (abs(corr_mat_codon) >=0.5) .* sign(corr_mat_codon) +...
            (abs(corr_mat_codon) <0.5 & abs(corr_mat_codon) >= 0.3 )  .* sign(corr_mat_codon) * 0.5;
        
        for i_para = 1:size(corr_mat_para,3)
            corr_mat_para(i_feature,:,i_para) = corr_mat_codon(i_para,:);
            feature_corr_mat_para(i_feature,:,i_para) = i_feature;
            for i_ligand = 1:length(ligand_vec)
                data_color(i_feature,i_ligand,i_para) = i_ligand;
            end
        end
        
 
        
    end
    

    
    %% Figure scatter plots for spearman corr
    % Figure 3B: 
% 
% tested 05/12/2024, Matlab 2020a

    if 1
        
        Y_scatter = corr_mat_para(:);
        X_scatter = feature_corr_mat_para(:);
        y_num_bins = 5;
        x_max_width = 0.7;
        
        scatter_args = { 9,color_vec(data_color(:),:),'filled'};
        % scatter_args = { 9,[0,0,0],'filled'};
        
        [fig_gcf,fig_gca] = simple_violin_scatter( X_scatter, Y_scatter, y_num_bins, x_max_width, scatter_args );
        
        set(fig_gcf, 'PaperUnits','points')
        paper_pos = [0,0,[150,100]]*2;
        paper_size = [150,80]*2;
        set(fig_gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        %set(fig_gca,'XTick',1:6,'XTickLabel',{},'YTickLabel',{})
        set(fig_gca,'XTick',1:6,'XTickLabel',{'Speed','Peak','Duration','Total','EvL','Osc'})
        ylabel('spearman correlation')

        xlim([0,7]);
        saveas(fig_gcf,strcat(fig_save_path,'Fig3A_para_codon_spearman_corr_scatter'),'epsc')
        %saveas(fig_gcf,strcat(fig_save_path,'Fig3A_para_codon_spearman_corr_scatter_bk'),'epsc')
        close()
    end
    

end

