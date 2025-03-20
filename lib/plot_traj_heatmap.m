function [] = plot_traj_heatmap(data,vis_data_field,fig_opts)

% data.info_ligand, data.info_dose_str are required for naming the saved
% figures

if isfield(fig_opts,'fig_savepath')
    
    savepath = fig_opts.fig_savepath;
else
    savepath = './';
    
end

if isfield(fig_opts,'color_limits')
    color_limits = fig_opts.color_limits;
else
    color_limits = [-0.01,0.25]; %for NFkB
end

if isfield(fig_opts,'fig_name')
    save_fig_name = fig_opts.fig_name;
else
    save_fig_name = '';
end

if isfield(fig_opts,'fig_size')
    paperpos=fig_opts.fig_size.paperpos;
    papersize=fig_opts.fig_size.papersize;
    draw_pos=fig_opts.fig_size.draw_pos;
else
    paperpos=[0,0,100,130];
    papersize=[100 130];
    draw_pos=[10,10,90,120];
end


for i_sti = 1:length(data.(vis_data_field{1}))
    for i_data_field = 1:length(vis_data_field)
        figure(i_sti)
        
        cell_num=length(data.order{i_sti});
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(data.(vis_data_field{i_data_field}){i_sti}(data.order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limits);%[-0.001,0.2] for TNF
        
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
        saveas(gcf,strcat(savepath,save_fig_name,data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/',''),'_',vis_data_field{i_data_field}),'epsc');
        close
    end
end



end