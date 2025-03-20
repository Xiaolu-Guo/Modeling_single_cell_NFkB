function [] = draw_spearman_corr_each_codon_2023_06(data_save_file_path,fig_save_path)
% For Guo et al. Figure S3B, spearman correlation between signaling codons and
% parameters
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


if 1
    corr_mat_para = NaN(6,5,7);
    
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
            corr_mat_para(i_feature,:,i_para) = corr_rel(i_para,:);
        end
        
 
        
        %% Figure S3B codon-para spearman corr
        if 1
            figure(2)
            set(gcf, 'PaperUnits','points')
            
            paper_pos = [0,0,parper_size_codon_para];
            paper_size = parper_size_codon_para;
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
            h = heatmap( sti_vec2,para_names,corr_rel,'Colormap',mymap,'CellLabelColor','none');%'k' % corr_mat_codon
            % xlabel(features{i_feature})
            caxis([-1,1])
            %h.CellLabelColor = [0,0,0];
            h.CellLabelFormat = '%0.2g';
            
            colorbar('off')
            set(gca,'fontsize',7,'fontname','Arial');
            
            switch features{i_feature}
                case 'OscVsNonOsc'
                    for i_index = 1:length(h.XDisplayLabels)
                        h.XDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                    end
                otherwise
                    for i_index = 1:length(h.XDisplayLabels)
                        %h.XDisplayLabels{i_index} = ['\color[rgb]{1,1,1}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                        h.XDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                    end
            end
            
            for i_index = 1:length(h.YDisplayLabels)
                h.YDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.YDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
            end
            
            saveas(gcf,strcat(fig_save_path,features{i_feature},'_para_corr_high_dose'),'epsc')
            close()
        end
        
    end
 
end

