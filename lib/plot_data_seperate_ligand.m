function [] = plot_data_seperate_ligand(data,vis_data_field,savepath)

paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

switch nargin
    case 1
        vis_data_field=fieldnames(data);
end


for i_sti = 1:length(data.(vis_data_field{1}))
    for i_data_field =1:length(vis_data_field)
        figure(i_sti)
        
        cell_num=length(data.order{i_sti});
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(data.(vis_data_field{i_data_field}){i_sti}(data.order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        
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
        saveas(gcf,strcat(savepath,data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/',''),'_',vis_data_field{i_data_field}),'epsc');
        close
    end
end



end