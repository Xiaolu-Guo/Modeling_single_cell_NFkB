% For Guo et al. Figure 2D & Figure S4D, Machine learning results
% resutls from python codes, for random forest classifier using signaling
% codon as input and ligand infomation as output
% 
% tested 05/30/2024, Matlab 2020a

i_data_set = 1;
data_name{i_data_set} = 'Ade';
data_set_all = {'exp','fitting','sampling'};
method_all = {'random_forest','random_forest','random_forest'};

confusion_mat_all{1} = [530  29  12   8  53
  38 327  42 110  70
  22 101 113  61  69
  16 168  45 152  31
  33  98  28  24 365];% random forest

confusion_mat_all{2} = [605  11   4   2  10
  22 391  36  99  39
   4  73 160  63  66
   4 148  40 185  35
   3  47  46  17 435];% random forest Fitting

confusion_mat_all{3} =  [581   5   1   7   6
  27 383  76  88  26
   6  75 332 100  87
  18  69  90 390  33
   5  26  81  15 473];% random forest Sampling

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

