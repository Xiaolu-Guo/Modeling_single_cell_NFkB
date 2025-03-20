
codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
dose_label = {'L','M','H'};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);


if 0
    % clear all
    load('mutual_info_cal_data_example.mat')
    
    data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))
    
    i_sti_exp = 1;
    data_save_path = './python/';
    
    for i_ligand_exp = 1:length(ligand_vec)%TNF,LPS [27];%
        
        nfkb_codon_exp_all = [];
        nfkb_exp_id_all = [];
        
        nfkb_codon_sim_all = [];
        nfkb_sim_id_all = [];
        
        for i_dose_exp = 1:length(dose_vec{i_ligand_exp})
            
            nfkb_codon_exp = [];
            nfkb_exp_id = [];
            
            nfkb_codon_sim = [];
            nfkb_sim_id = [];
            
            for i_codon = 1:length(codon_list)
                i_data_exp = find(categorical(data.info_ligand)==ligand_vec{i_ligand_exp} & categorical(data_dose_str)==dose_vec{i_ligand_exp}{i_dose_exp});
                nfkb_codon_exp = [nfkb_codon_exp,collect_feature_vects.(codon_list{i_codon}){i_data_exp*2-1}];
                nfkb_codon_sim = [nfkb_codon_sim,collect_feature_vects.(codon_list{i_codon}){i_data_exp*2}];
                
            end
            
            nfkb_exp_id = (i_dose_exp-1)*ones(size(nfkb_codon_exp,1),1);
            nfkb_codon_exp_all = [nfkb_codon_exp_all;nfkb_codon_exp];
            nfkb_exp_id_all = [nfkb_exp_id_all;nfkb_exp_id];
            
            nfkb_sim_id = (i_dose_exp-1)*ones(size(nfkb_codon_sim,1),1);
            nfkb_codon_sim_all = [nfkb_codon_sim_all;nfkb_codon_sim];
            nfkb_sim_id_all = [nfkb_sim_id_all;nfkb_sim_id];
            
            
        end
        
        writematrix(nfkb_codon_exp_all,strcat(data_save_path,'X_exp_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
        writematrix(nfkb_exp_id_all,strcat(data_save_path,'y_exp_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
        
        writematrix(nfkb_codon_sim_all,strcat(data_save_path,'X_sim_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
        writematrix(nfkb_sim_id_all,strcat(data_save_path,'y_sim_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
        
    end
    
end

if 1
    
    ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
    data_vec = {'exp','sim'};
    fig_save_path = '../SubFigures2023/';

    for i_ligand_exp = 1:length(ligand_vec)
  
        for i_data_vec = 1:length(data_vec)
%             readmatrix(strcat(data_save_path,'XGBClassification_',data_vec{i_data_vec},'_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
            confusion_mat = readmatrix(strcat(data_save_path,'XGBClassification_',data_vec{i_data_vec},'_dose_codon_stim_',ligand_vec{i_ligand_exp},'.csv'));
            confusion_mat = confusion_mat';% check ?? 
            confusion_score = confusion_mat./(sum(confusion_mat,2)*ones(1,size(confusion_mat,2)));
            
            figure(2)
            set(gcf, 'PaperUnits','points')
            mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
                ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
            
            paper_pos = [0,0,160,120]*0.75;
            paper_size = [160,120]*0.75;
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
            h = heatmap( dose_label,dose_label,confusion_score,'Colormap',mymap,'CellLabelColor','k');%'none'
            caxis([0,1])
            
            h.CellLabelColor = [0,0,0];
            h.CellLabelFormat = '%0.2g';
            h.FontColor = [0,0,0];
            
            set(gca,'FontSize',7,'FontName','Arial')
            saveas(gcf,strcat(fig_save_path,'percision_',data_vec{i_data_vec},'_',ligand_vec{i_ligand_exp}),'epsc')
            close()
        end
        
    end
end

