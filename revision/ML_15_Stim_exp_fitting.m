data_name_vec = {'Fitting', 'Ade'};
method_all = {'randomforest','randomforest'};
fig_save_path = '../SubFigures2023/';
% for i_data = 
% 
% end

for i_data_set = 1:length(data_name_vec)
    SaveFileName = strcat('./revision/data_signaling_codons_Machine_learning_format/','All_dose_RandomForest_conf_mat_',data_name_vec{i_data_set}, '.csv');
    conf_mat = readmatrix(SaveFileName);

    confusion_mat = conf_mat';
    data_set = data_name_vec{i_data_set};
    method = method_all{i_data_set};
    confusion_score = confusion_mat./(sum(confusion_mat,2)*ones(1,size(confusion_mat,2)));
    
    ligand_all = {'TNF-L','TNF-M','TNF-H','Pam-L','Pam-M','Pam-H','CpG-L','CpG-M','CpG-H','LPS-L','LPS-M','LPS-H','pIC-L','pIC-M','pIC-H'};
    %index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
    figure(2)
    set(gcf, 'PaperUnits','points')
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    paper_pos = [0,0,150,120]*1.3;
    paper_size = [150,120]*1.3;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( ligand_all,ligand_all,confusion_score,'Colormap',mymap,'CellLabelColor','k');%'none'
    caxis([0,1])
    % h.FontColor = [0,0,0];
    
    %         for i_index = 1:length(index)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.8,0.8,0.8}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8} {red}
    %         end
    %
    %         for i_index = 1:length(index_non_wide)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.4,0.4,0.4}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
    %         end
    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    h.FontColor = [0,0,0];

    set(gca,'FontSize',7,'FontName','Arial')
    saveas(gcf,strcat(fig_save_path,'All_dose_percision_',data_set,'_',method),'epsc')
    close()
    
    
    confusion_mat_2 = conf_mat;
    sensitivity_score = confusion_mat_2./(sum(confusion_mat_2,2)*ones(1,size(confusion_mat_2,2)));
    
    ligand_all = {'TNF-L','TNF-M','TNF-H','Pam-L','Pam-M','Pam-H','CpG-L','CpG-M','CpG-H','LPS-L','LPS-M','LPS-H','pIC-L','pIC-M','pIC-H'};
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    paper_pos = [0,0,150,120]*1.3;
    paper_size = [150,120]*1.3;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( ligand_all,ligand_all,sensitivity_score,'Colormap',mymap,'CellLabelColor','k');%'none'
    caxis([0,1])

    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    h.FontColor = [0,0,0];

    set(gca,'FontSize',7,'FontName','Arial')
    saveas(gcf,strcat(fig_save_path,'All_dose_sensitivity_',data_set,'_',method),'epsc')
    close()
end