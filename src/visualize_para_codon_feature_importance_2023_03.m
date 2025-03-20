function [] = visualize_para_codon_feature_importance_2023_03(data_save_file_path,data_save_path,fig_save_path)
% For Guo et al. Figure S3C, feature importance of parameter calculated
% from XGBoost regression
% 
% tested 05/12/2024, Matlab 2020a

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'};
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'10ug/mL';'33ug/mL';'100ug/mL'}};
dose_symbol = {'Low','Med','High'};

dose_vec = {{'10ng/mL'};
    {'1ug/mL'};
    {'333nM'};%'10nM';;'1uM'
    {'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'100ug/mL'}};
dose_symbol = {'High'};

data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

mymap = [ ones(21,1),(1:-0.05:0)',(1:-0.05:0)'];


features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};

stat_methods_vec = {'LinearRegression','RamdomForestRegression','XGBRegression'};

ColorScaling_ = {'log','scaled','scaled'};
c_axis_limit = {[1e-3,1e5],[0,0.8],[0.1,0.8]};

para_names = {'rcp syn','C deg','end','TAK ac', 'time d1','time d2','Tot NFkB'};

for i_stat_methods = length(stat_methods_vec)
    stat_methods = stat_methods_vec{i_stat_methods};
    score_ = NaN(length(features),5*length(dose_symbol));
    rand_feature_imp = NaN(length(features),5*length(dose_symbol));
    for i_feature = 1:length(features)
        corr_mat_codon = NaN(7,5*length(dose_symbol));
        i_cond = 1;
        for i_sti = 1:length(ligand_vec)
            ligand = ligand_vec{i_sti};
            % dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
            for i_dose = 1:length(dose_vec{i_sti})
                dose  = dose_vec{i_sti}{i_dose};
                dose_sym = dose_symbol{i_dose};
                
                A = readmatrix(strcat(data_save_path,stat_methods,'_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
                A_rand = readmatrix(strcat(data_save_path,stat_methods,'_rand_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
                
                score_(i_feature,i_cond) = A(1);
                len_rcp = length(A)-5;
                rand_feature_imp(i_feature,i_cond) = mean(A_rand(length(A)+1:end));
                corr_mat_codon(4:end,i_cond) = A(2:5);%core
                corr_mat_codon(1:len_rcp,i_cond) = A(6:5+len_rcp);%receptor
                sti_vec{i_cond} = strcat(ligand,'-',dose_sym);
                i_cond = i_cond+1;
            end
        end        
        
        % Figure S3C
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        %         paper_pos = [0,0,450,200];
        %         paper_size = [450,200];
        paper_pos = [0,0,200,200];
        paper_size = [200,200];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,para_names,corr_mat_codon,'Colormap',mymap,...
            'ColorScaling',ColorScaling_{i_stat_methods},...
            'CellLabelColor','none');%'k' % corr_mat_codon
       
        colorbar('off')

        for i_index = 1:length(h.XDisplayLabels)
            h.XDisplayLabels{i_index} = ['\color[rgb]{1,1,1}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
        end
        % xlabel(features{i_feature})
        if i_stat_methods>1
            caxis(c_axis_limit{i_stat_methods})
        end
        %h.CellLabelColor = [0,0,0];
        h.CellLabelFormat = '%0.2g';
        %set(gca,'fontsize',6,'fontweigt','b');

        saveas(gcf,strcat(fig_save_path,stat_methods,'_',features{i_feature},'_para_feature_importance'),'epsc')
        close()
         
    end
    
end

