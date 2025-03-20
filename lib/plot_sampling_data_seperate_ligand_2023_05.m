function [] = plot_sampling_data_seperate_ligand_2023_05(data,vis_data_field,savepath,filter_TNF)

paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;
fitting_dose = '';%'_alldose' '_lowdose'

switch nargin
    case 1
        vis_data_field=fieldnames(data);
end

% input
if nargin<4
    filter_TNF =1;
end
if filter_TNF % whether filter out TNF or not
    thresh_TNF=0.33;
    
    % thresh_field = {'model_sim_cv_order'};%pred_mode_cv_order osc
    for i_data = 1:3
        index = data.cv_order{i_data}(1: ceil(length(data.cv_order{i_data}) * thresh_TNF));
        data.(vis_data_field{1}){i_data} = data.(vis_data_field{1}){i_data}(index ,:);
        [~, data.order{i_data}]=sort(max(data.(vis_data_field{1}){i_data},[],2),'descend');
    end
    
    for i_data = 4:length( data.(vis_data_field{1}))
        [~, data.order{i_data}]=sort(max(data.(vis_data_field{1}){i_data},[],2),'descend');
    end
else
    for i_data = 1:length( data.(vis_data_field{1}))
        [~, data.order{i_data}]=sort(max(data.(vis_data_field{1}){i_data},[],2),'descend');
    end
end

for i_sti = 1:length(data.(vis_data_field{1}))
    for i_data_field = 1:length(vis_data_field)
        figure(i_sti)
        
        cell_num=length(data.order{i_sti});
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(data.(vis_data_field{i_data_field}){i_sti}(data.order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        %
        XLabels = 0:5:((size(data.(vis_data_field{i_data_field}){i_sti},2)-1)*5);
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
        
        %         xlabel('Time (hours)');
        %         ylabel(vis_data_field{i_data_field});
        % clb=colorbar;
        %
        % clb.Label.String = 'NFkB(A.U.)';
        colorbar('off')
        
        set(gca,'fontsize',14,'fontname','Arial');
        saveas(gcf,strcat(savepath,data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/',''),'_',vis_data_field{i_data_field},fitting_dose),'epsc');
        close
    end
end



end