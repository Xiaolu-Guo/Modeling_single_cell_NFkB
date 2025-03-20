load('../raw_data/link_cell_metrics_TNF_LPS.mat')

metric_fields = {'duration','oscpower','max_pos_integral','pos_pk1_time','max_value','time2HalfMaxPosIntegral'};% ,'time2HalfMaxValue'};

    
metric_spe = {[2],[],[],[],[],[],[]};
cell_diff = zeros(2000,length(metric_fields));
cell_diff_unmatch = zeros(2000,length(metric_fields));
pperm = randperm(2000);
threshold_vec = [0.1,0.25,0.5];
axis_lim = [0,8,0,8;
    0,3e-3,0,3e-3;
    0,2,0,3;
    0,2.5,0,5;
    0,1,0,0.8;
    0,4.5,0,4.5];
cell_specificity_counts = zeros(length(threshold_vec),length(metric_fields));

xlim_num = [0,20;0,2e4;0,10;0,10;0,10;0,5;];% 0,20];
for i_metric = 1:length(metric_fields)
    if isempty(metric_spe{i_metric})
        cell_diff(:,i_metric) =abs( metrics{2}.(metric_fields{i_metric})...
            - metrics{1}.(metric_fields{i_metric}))./metrics{1}.(metric_fields{i_metric});
%         for i_thresh = 1:length(threshold_vec)
%             cell_specificity_counts(i_thresh,i_metric) = sum(abs(metrics{2}.(metric_fields{i_metric})...
%                 - metrics{1}.(metric_fields{i_metric}))./metrics{1}.(metric_fields{i_metric}) > threshold_vec(i_thresh))/2000;
%         end
         cell_diff_unmatch(:,i_metric) = abs( metrics{2}.(metric_fields{i_metric})(pperm)...
            - metrics{1}.(metric_fields{i_metric}))./metrics{1}.(metric_fields{i_metric});
        
        index_TNF = metrics{1}.(metric_fields{i_metric})' < axis_lim(i_metric,2);
        index_LPS = metrics{2}.(metric_fields{i_metric})' < axis_lim(i_metric,4);
        
        % subplot(3,length(metric_fields),i_metric )
        figure(1)
        scatter(metrics{1}.(metric_fields{i_metric}),metrics{2}.(metric_fields{i_metric}),3,[0,0,1],'filled');hold on
        [r,m,b] = regression(metrics{1}.(metric_fields{i_metric})(index_TNF&index_LPS)',metrics{2}.(metric_fields{i_metric})(index_TNF&index_LPS)');
    else
        cell_diff(:,i_metric) = abs(metrics{2}.(metric_fields{i_metric})(:,metric_spe{i_metric})...
            - metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}))./metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric});
%         for i_thresh = 1:length(threshold_vec)
%             cell_specificity_counts(i_thresh,i_metric) = sum(abs(metrics{2}.(metric_fields{i_metric})(:,metric_spe{i_metric})...
%                 - metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}))./metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}) > threshold_vec(i_thresh))/2000;
%         end
        cell_diff_unmatch(:,i_metric) = abs(metrics{2}.(metric_fields{i_metric})(pperm,metric_spe{i_metric})...
            - metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}))./metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric});
 
        index_TNF = metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}) < axis_lim(i_metric,2);
        index_LPS = metrics{2}.(metric_fields{i_metric})(:,metric_spe{i_metric}) < axis_lim(i_metric,4);
       
        % subplot(3,length(metric_fields),i_metric)
        figure(1)
        scatter(metrics{1}.(metric_fields{i_metric})(:,metric_spe{i_metric}),metrics{2}.(metric_fields{i_metric})(:,metric_spe{i_metric}),3,[0,0,1],'filled');hold on
        [r,m,b] = regression(metrics{1}.(metric_fields{i_metric})(index_TNF&index_LPS,metric_spe{i_metric})',metrics{2}.(metric_fields{i_metric})(index_TNF&index_LPS,metric_spe{i_metric})');
        
    end
    
    xlim(axis_lim(i_metric,1:2));
        ylim(axis_lim(i_metric,3:4));
    XL = xlim();

    plot(XL,m*XL + b,'k','LineWidth',1.5)
    % title(strcat('R =',num2str(r)));
    r
    % title(replace(metric_fields{i_metric},'_',' '))
    xlabel({'TNF'});%, replace(metric_fields{i_metric},'_',' ')})
    ylabel({'LPS'});%, replace(metric_fields{i_metric},'_',' ')})
    
    
    % set(gca,'FontSize',8)
    % set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)
    set(gca,'fontsize',7,'fontweight','b')
    Set_figure_size_square
    saveas(gcf,strcat('../SubFigures/TNF_LPS_scatter_',metric_fields{i_metric}),'epsc')
    close
    edges = linspace(xlim_num(i_metric,1),xlim_num(i_metric,2),20);
    % subplot(3,length(metric_fields),i_metric + length(metric_fields));hold on
    figure(2)
    histogram(cell_diff(:,i_metric),edges );hold on
    xlabel( {'relative' 'metric differences'})
    ylabel('counts')
    xlim(xlim_num(i_metric,:));
    YL = ylim();
    plot([prctile(cell_diff(:,i_metric),50),prctile(cell_diff(:,i_metric),50)],YL,'Color','g');
%     X = categorical({'>10%','>25%','>50%'});
%     X = reordercats(X,{'>10%','>25%','>50%'});
%     
%     subplot(3,length(metric_fields),i_metric + 2* length(metric_fields))
%     bar(X, cell_specificity_counts(:,i_metric) )
%     xlabel({'differences larger than'})
%     ylabel('percentage of cell response specifically')

    % subplot(3,length(metric_fields),i_metric +  length(metric_fields))
    histogram(cell_diff_unmatch(:,i_metric),edges );hold on
    % xlabel({'differences of',replace(metric_fields{i_metric},'_',' ')});
    ylabel('counts')
        xlim(xlim_num(i_metric,:));
    plot([prctile(cell_diff_unmatch(:,i_metric),50),prctile(cell_diff_unmatch(:,i_metric),50)],YL,'Color','r');
    
    set(gca,'fontsize',12,'fontweight','b')
    Set_figure_size_square
    set(gcf, 'PaperPosition', paperpos*2,'PaperSize', papersize*2,'Position',draw_pos*2)

    if i_metric == 2
        legend('matched(m)','median m','unmatched(u)','median u')
    end
    saveas(gcf,strcat('../SubFigures/TNF_LPS_metric_diff_hist_',metric_fields{i_metric}),'epsc')
    close
    
% prctile(cell_diff_unmatch(:,i_metric),50) - prctile(cell_diff(:,i_metric),50)
end


