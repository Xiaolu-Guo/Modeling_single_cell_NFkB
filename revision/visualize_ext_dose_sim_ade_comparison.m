


%% plot function

function fig_gcf = plot_heatmap(data_traj,data_order_index,color_limits)

if nargin <2
    data_first_half = data_traj(:,1:size(data_traj,2)/2);
    data_first_half(isnan(data_first_half)) = 0;
    [~,data_order_index] = sort(sum(data_first_half,2),"descend");
end

if nargin <3
    color_limits = [-0.001,0.3];
end

%
% p = inputParser;
% % Required: ID input
% %valid_data_traj = @(x) assert(isnumeric(input_data_traj));
% addRequired(p,'data_traj');
% % Optional parameters to be passed to metrics function
% % valid_color_limit = @(x) assert(isnumeric(x)&size(x) == [1,2]);
% addParameter(p,'color_limits',[-0.001,0.3]);%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
% addParameter(p,'data_order_index',data_order_index_default); %allows adjustment of convection shift (?)
%
% % Parse the inputs
% parse(p, data_traj, 'data_order_index', data_order_index, 'color_limits', color_limits);
% color_limits = p.Results.color_limits;
% data_order_index = p.Results.data_order_index;

figure(1)
paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

cell_num=size(data_traj,1);
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

% subplot(1,length(vis_data_field),i_data_field)
h=heatmap(data_traj(data_order_index,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limits);%[-0.001,0.2] for TNF

XLabels = 0:5:((size(data_traj,2)-1)*5);
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
fig_gcf = gcf;

end