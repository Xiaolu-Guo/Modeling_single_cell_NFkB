function [] = violin_plot_codon_2023_CRI(collect_feature_vects,figure_option)
%repeated codes:
%% calculate features of replicates

% black represents the first type, experimental data
% red represents the second type, predicted (mode) data

FeatureListFile = 'FeatureList_v2.xlsx';
FeatureListTable = readtable(FeatureListFile, 'Sheet', 'Codewords');

CodewordList = table2cell(FeatureListTable(:,3));
codewords = unique(CodewordList);
codewords = codewords([3,5],1);
ligand_all = unique(collect_feature_vects.info_ligand,'stable');
if length(ligand_all) == 5
    ligand_all = {'LPS'};
end
units_all = {'\mu g/mL'};
dose_val = { {'10'}};

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

if nargin ==2
    if isfield(figure_option,'paper_opt')
        paperpos = figure_option.paper_opt.paperpos;
        papersize = figure_option.paper_opt.papersize;
    end
else
    paperpos=[0,0,1200,700];
    papersize=[1200,700];
end
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)

for i=1:length(codewords)
    vects=collect_feature_vects.(codewords{i});
    position_subfig=[-1,-1,-1,subplot_height*0.85];
    vects = vects([21,22,29,30],1);
    
    id_begin=1;
    
    for jj=1:length(ligand_all)
        
        %% set the value for posi left
        position_subfig(1)=(jj-1+0.15+width_shift)*subplot_width;%*(jj~=1)
        
        %% set the value for posi bottom
        switch codewords{i}
            case 'Speed'
                heigh_pos=0;
                plot_ornot=1;
            case 'OscVsNonOsc'
                heigh_pos=1;
                plot_ornot=1;
        end
        
        if plot_ornot
            
            
            position_subfig(2)=(heigh_pos+0.1+height_shift)*subplot_height;
            
            %% set the value for width height
            position_subfig(3)=subplot_width * 0.85 ;%+ width_shift*(jj==1)
            %             if strcmp(codewords{i},'OscVsNonOsc')
            %                 position_subfig(4)=subplot_height*1.9;
            %             end
            
            %% plot at specific position
            ax = subplot('Position', position_subfig);
            %ax=gca;
            id_begin = find(strcmp(collect_feature_vects.info_ligand,ligand_all{jj}),1,'first');
            id_end = find(strcmp(collect_feature_vects.info_ligand,ligand_all{jj}),1,'last');
            id_begin = floor((id_begin+id_end)/2);
            id_end = id_begin+1;
%             id_end=id_begin+sum(strcmp(collect_feature_vects.info_ligand,ligand_all{jj}))-1;
%             ligand_dose_label=collect_feature_vects.info_dose_str...
%                 (strcmp(collect_feature_vects.info_data_type,data_type_all{1})...
%                 & strcmp(collect_feature_vects.info_ligand,ligand_all{jj}));% LPS
            ligand_dose_label= dose_val{jj};
            
            buffer = 0.1*max(cell2mat(vects(:)));
            %ylimits = [min(cell2mat(vects(:)))-buffer max(cell2mat(vects(:)))+buffer];
            ylimits = [0-buffer 1+buffer];
            %violin(vects, [1:length(ids)]*2, 'Axes', ax, 'YLim', ylimits);
id_begin =1;
id_end = 4;
            id_vec = id_begin:id_end;
% id_vec = [21,22,29,30,21,22,21,22 ]
            violin(data_type,vects(id_vec), id_vec, 'Axes', ax, 'YLim', ylimits);
            
            if heigh_pos == 0
                %title(codewords{i})
                set(ax,'XTick',(id_begin+.5):2:id_end,'XTickLabel',ligand_dose_label);
                ax.XTickLabelRotation = 45;
                xlabel({ligand_all{jj},strcat('(',units_all{jj},')')});
            else
                set(ax,'XTick',[],'XTickLabel',{});
            end
            
            if jj==1

                    ylabel(codewords{i});
                    %'Min-Max-Scaled Zscore'

            else
                set(ax,'YTick',[],'YTickLabel',{});
            end
            
            set(gca,'FontSize',10,'FontWeight','b')
%             id_begin=id_end+1;
            
        end
        
    end
    
end

if nargin == 2
    if isfield(figure_option,'save_file')
        saveas(gcf,figure_option.save_file,'epsc')
        close
    end
    
end