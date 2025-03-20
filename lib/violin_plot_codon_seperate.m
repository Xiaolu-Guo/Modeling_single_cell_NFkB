function [] = violin_plot_codon_seperate(collect_feature_vects,save_path)
%repeated codes:
%% calculate features of replicates

% black represents the first type, experimental data
% red represents the second type, predicted (mode) data

FeatureListFile = 'FeatureList_v2.xlsx';
FeatureListTable = readtable(FeatureListFile, 'Sheet', 'Codewords');

CodewordList = table2cell(FeatureListTable(:,3));
codewords = unique(CodewordList);
ligand_all = unique(collect_feature_vects.info_ligand);
if length(ligand_all) == 5
    ligand_all = {'TNF';'CpG';'Pam3CSK';'LPS';'PolyIC'};
end
data_type_all = unique(collect_feature_vects.info_data_type);
data_type = 1:length(data_type_all);
%1 represents exp
%2 represents pred

%% subplot size

height_shift=0.7;
width_shift=0.25;
subplot_width=1/(length(ligand_all)+width_shift);
subplot_height=1/(length(codewords)+height_shift);

%%

figure(1)

paperpos=[0,0,200,150];
papersize=[200,150];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)

for i=1:length(codewords)
    vects=collect_feature_vects.(codewords{i});
    
    
    id_begin=1;
    
    for jj=1:length(ligand_all)
        
        %% set the value for posi left
        
        %% set the value for posi bottom
        
            
            
            
            %% set the value for width height
            %             if strcmp(codewords{i},'OscVsNonOsc')
            %                 position_subfig(4)=subplot_height*1.9;
            %             end
            
            %% plot at specific position
            figure(1)
            set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)

            ax = gca;
            %ax=gca;
            id_end=id_begin+sum(strcmp(collect_feature_vects.info_ligand,ligand_all{jj}))-1;
            ligand_dose_label=collect_feature_vects.info_dose_str...
                (strcmp(collect_feature_vects.info_data_type,data_type_all{1})...
                & strcmp(collect_feature_vects.info_ligand,ligand_all{jj}));% LPS
            
            buffer = 0.1*max(cell2mat(vects(:)));
            %ylimits = [min(cell2mat(vects(:)))-buffer max(cell2mat(vects(:)))+buffer];
            ylimits = [0-buffer 1+buffer];
            %violin(vects, [1:length(ids)]*2, 'Axes', ax, 'YLim', ylimits);

            violin(data_type,vects(id_begin:id_end), id_begin:id_end, 'Axes', ax, 'YLim', ylimits);
            
%             if heigh_pos == 0
                %title(codewords{i})
                set(ax,'XTick',(id_begin+.5):2:id_end,'XTickLabel',ligand_dose_label);
                ax.XTickLabelRotation = 45;
                xlabel(ligand_all{jj});
%             else
%                 set(ax,'XTick',[],'XTickLabel',{});
%             end
            
%             if jj==1

                    ylabel(codewords{i});
                    %'Min-Max-Scaled Zscore'

%             else
%                 set(ax,'YTick',[],'YTickLabel',{});
%             end
            
            set(gca,'FontSize',10,'FontWeight','b')
            id_begin=id_end+1;
                saveas(gcf,strcat(save_path,'codon_',ligand_all{jj},'_',codewords{i}),'epsc')
    close
        
    end
    
end
