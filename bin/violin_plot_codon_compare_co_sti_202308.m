function [] = violin_plot_codon_compare_co_sti_202308(collect_feature_vects,figure_option)

% check whether figure_option has the codons to plot

%% calculate features of replicates

% black represents the first type, experimental data
% red represents the second type, predicted (mode) data

FeatureListFile = 'FeatureList_v2.xlsx';
FeatureListTable = readtable(FeatureListFile, 'Sheet', 'Codewords');

CodewordList = table2cell(FeatureListTable(:,3));
codewords = unique(CodewordList);
if isfield(figure_option,'codon')
    codewords = figure_option.codon;
end
fields_name = fieldnames(collect_feature_vects);
% index_valid = [1:12,19:24,27:38];
% for i_fields = 1:length(fields_name)
%     collect_feature_vects.(fields_name{i_fields}) =  collect_feature_vects.(fields_name{i_fields})(index_valid);
% end

ligand_all = unique(collect_feature_vects.info_ligand,'stable');
if length(ligand_all) == 5
    % ligand_all = {'TNF';'CpG';'Pam3CSK';'LPS';'PolyIC'};
    ligand_all = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
    
end

data_type_all = unique(collect_feature_vects.info_data_type,'stable');
data_type = 1:length(data_type_all);
data_type_each_plot = 1;
%1 represents exp
%2 represents pred

%% subplot size

height_shift=0.7;
width_shift=0.25;
subplot_width=1/(length(ligand_all)+width_shift);
subplot_height=1/(length(codewords)*length(data_type)+height_shift);

%%


if nargin ==2
    if isfield(figure_option,'paper_opt')
        paperpos = figure_option.paper_opt.paperpos;
        papersize = figure_option.paper_opt.papersize;
    end
    
    if isfield(figure_option,'distri_color')
        colors = figure_option.distri_color;
    else
        colors = {[0,1,0]};
        
    end
else
    paperpos=[0,0,2000,1000];
    papersize=[2000,1000];
end

for i=1:length(codewords)
    vects = collect_feature_vects.(codewords{i});
    
    ylimits = [-1.2,1.2];
    
    for i_data_type = 1:length(data_type)
        
        figure(1)
        
        ax = gca;
        id_begin = 1;
        id_end = length(collect_feature_vects.info_ligand);
        
        plot_pos = id_begin:id_end;
        plot_pos(1) = plot_pos(1)+0.5;
        plot_pos(3) = plot_pos(3)-0.5;
        
        violin_20230818(data_type_each_plot,vects(id_begin:id_end), plot_pos, 'Axes', ax, 'YLim', ylimits,'Color',colors(i_data_type),'Area',0.08,'MarkerSize',3);
        
        set(ax,'XTick',(id_begin+.5):2:id_end,'XTickLabel',{});
        XL = xlim();
        xlim([XL(1)-0.2,XL(2)+0.2])
        set(ax,'YTickLabel',{});%'YTick',[],
        
        set(gca,'FontSize',7,'FontName','Arial')
        % Set the box on so you can see all the edges
        set(gca, 'Box', 'on');
        % Change the line width of the left and bottom edges
        set(gca, 'LineWidth', 1);
        % Get the children of the current axes
        axChildren = get(gca, 'Children');
        
        % Change the line width of the top and right edges
        for child = 1:length(axChildren)
            if isa(axChildren(child), 'matlab.graphics.axis.decorator.Ruler')
                axChildren(child).LineWidth = 1;
            end
        end
        %             id_begin=id_end+1;
    end
    
    if nargin == 2
        if isfield(figure_option,'save_file')
            set(gcf, 'PaperUnits','points')
            figPosition = paperpos;
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position', figPosition)
            
            saveas(gcf,strcat(figure_option.save_file,'_',collect_feature_vects.info_ligand{end},'_',codewords{i}),'epsc')
            close
        end
        
    end
    
end

