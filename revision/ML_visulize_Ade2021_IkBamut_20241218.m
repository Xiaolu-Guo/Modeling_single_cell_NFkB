data_name_vec = {'WT','WT_Train_IkBamut'};
method_all = {'randomforest','randomforest'};
fig_save_path = '../SubFigures2023/revision_2025/';
% for i_data = 
% 
% end

for i_data_set = 1:length(data_name_vec)
    SaveFileName = strcat('./data_signaling_codons_Machine_learning_format/','Ade2021exp_RandomForest_conf_mat_',data_name_vec{i_data_set}, '.csv');
    conf_mat = readmatrix(SaveFileName);

    confusion_mat = conf_mat';
    data_set = data_name_vec{i_data_set};
    method = method_all{i_data_set};
    confusion_score = confusion_mat./(sum(confusion_mat,2)*ones(1,size(confusion_mat,2)));
    
    ligand_all = {'LPS','pIC','TNF'};
    ligand_order = [3,1,2];
    %index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
    figure(2)
    set(gcf, 'PaperUnits','points')
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    paper_pos = [0,0,150,120]*0.9;
    paper_size = [150,120]*0.9;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( ligand_all(ligand_order),ligand_all(ligand_order),confusion_score(ligand_order,ligand_order),'Colormap',mymap,'CellLabelColor','k');%'none'
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
    saveas(gcf,strcat(fig_save_path,'Ade2021exp_percision_',data_set,'_',method),'epsc')
    close()
    
    
    confusion_mat_2 = conf_mat;
    sensitivity_score = confusion_mat_2./(sum(confusion_mat_2,2)*ones(1,size(confusion_mat_2,2)));
    
      ligand_all = {'LPS','pIC','TNF'};

    figure(2)
    set(gcf, 'PaperUnits','points')
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    paper_pos = [0,0,150,120]*0.9;
    paper_size = [150,120]*0.9;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( ligand_all(ligand_order),ligand_all(ligand_order),sensitivity_score(ligand_order,ligand_order),'Colormap',mymap,'CellLabelColor','k');%'none'
    caxis([0,1])

    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    h.FontColor = [0,0,0];

    set(gca,'FontSize',7,'FontName','Arial')
    saveas(gcf,strcat(fig_save_path,'Ade2021exp_sensitivity_',data_set,'_',method),'epsc')
    close()
end