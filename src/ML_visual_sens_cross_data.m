% For Guo et al. Figure 2D & Figure S4D, corresponding sensitivity, and cross validation on different data sets
% Machine learning results
% resutls from python codes, for random forest classifier using signaling
% codon as input and ligand infomation as output
% 
% tested 05/30/2024, Matlab 2020a

i_data_set = 1;
data_name{i_data_set} = 'Ade';
data_set_all = {'Fitting_train_exp_test','exp_train_fitting_test'};
method_all = {'random_forest','random_forest'};
clear confusion_mat_all
confusion_mat_all{1} = [195   5   1   0   0
   8 127  12  33  13
   2  22  71  23  22
   1  48  13  63  10
   1  13  14   4 139];% random forest Fitting train, exp test

confusion_mat_all{2} = [172   6   2   4  17
  13 103  16  46  15
   7  26  51  29  27
   2  55  13  53  12
  15  30   5   9 112];% random forest exp train, fitting test


for i_data_set = 1:length(confusion_mat_all)
    confusion_mat = (confusion_mat_all{i_data_set})';
    data_set = data_set_all{i_data_set};
    method = method_all{i_data_set};
    confusion_score = confusion_mat./(sum(confusion_mat,2)*ones(1,size(confusion_mat,2)));
    
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    
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
    saveas(gcf,strcat(fig_save_path,'percision_',data_set,'_',method),'epsc')
    close()
    
    confusion_mat_2 = (confusion_mat_all{i_data_set});
    sensitivity_score = confusion_mat_2./(sum(confusion_mat_2,2)*ones(1,size(confusion_mat_2,2)));
    
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    
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
    saveas(gcf,strcat(fig_save_path,'sensitivity_',data_set,'_',method),'epsc')
    close()
end

% fscore(i_ligand) = confusion_mat(

