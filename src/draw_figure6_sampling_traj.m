% visulize codons for prediction and Supriya's data
load('./Supriya_data_metrics_figure6.mat')
supriya_data_metrics = metrics;
supriya_data_codon = collect_feature_vects;
supriya_data = data;

if 0 %debug
    load('Supriya_data_metrics_max_pos_integral.mat')
    supriya_data_metrics = metrics;
    supriya_data = data;
    
end
for i_cond = 1:length(supriya_data_metrics)
    supriya_data_metrics{i_cond}.duration = metrics{i_cond}.duration(:,2);
end

load('./sampling_supriya_data_figure6.mat')
sampling_supriya_data_metrics = metrics;
sampling_supriya_data_codon = collect_feature_vects;
sampling_supriya_data = data;
for i_cond = 1:length(sampling_supriya_data_metrics)
    sampling_supriya_data_metrics{i_cond}.duration = metrics{i_cond}.duration(:,3);
end
% OneDrivePath = getenv('OneDrive');
% load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\TNFKO paper\AllFigures\Modeling\Updated modeling figure\tnfop100o_oscpower.mat'])
metrics_to_vis = {'duration','oscpower','max_pos_integral','pos_pk1_time','max_value','time2HalfMaxPosIntegral','time2HalfMaxValue','peakfreq'};
% metrics_to_vis = {'max_pos_integral'};
% metrics_to_vis = {'duration'};

if draw_data
    for i_data = 1:length(sampling_supriya_data.model_sim)
        data_to_draw = sampling_supriya_data.model_sim{i_data};
        data_exp = supriya_data.model_sim{i_data};
        index_order_to_draw = 0;
        rest_index = 1:size(data_to_draw,1);
        for i_cell = 1:size(data_exp,1)
            index_time = find(isnan(data_exp(i_cell,1:size(data_to_draw,2))));
            if isempty(index_time)
                index_time = size(data_to_draw,2);
            else
                index_time = index_time -1;
            end
            diff = ones(length(rest_index),1)*data_exp(i_cell,1:index_time) - data_to_draw(rest_index,1:index_time);
            diff_traj = sum(diff.^2,2);
            index_order_to_draw(i_cell) = find(diff_traj== ones(size(diff_traj))*min(diff_traj));
            rest_index = setdiff(rest_index,index_order_to_draw(i_cell),'stable');
            if isempty(rest_index)
                break
            end
        end
        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        h=heatmap(data_to_draw(index_order_to_draw,:),'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.3]);
        cell_num = length(index_order_to_draw);
        % subplot(1,length(vis_data_field),i_data_field)
        % h=heatmap(data.(vis_data_field{i_data_field}){i_sti}(data.order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        
        XLabels = 0:5:((size(data_to_draw,2)-1)*5);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/60);
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:cell_num;
        % Convert each number in the array into a string
        YCustomXLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        YCustomXLabels(:) = " ";
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.YDisplayLabels = YCustomXLabels;
        
        % xlabel('Time (hours)');
        % ylabel(vis_data_field{i_data_field});
        % clb=colorbar;
        % clb.Label.String = 'NFkB(A.U.)';
        colorbar('off')
        
        set(gca,'fontsize',14,'fontname','Arial');
        %         saveas(gcf,strcat('./figures/',replace(filelist{i_file},{'./Supriya_data/2020-11-','/AllMeasurements.mat'},{'',''}),'_sortAmplitude'),'epsc');
        close
    end
end