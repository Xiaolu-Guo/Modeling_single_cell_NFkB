filename = 'SDE_simulation/data/concentration_noise_df.csv';
data_traj = readmatrix(filename);

data_ordered = sort(data_traj,'descend');
paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

savepath = '../SubFigures2023/';
i_sti = 1;
ligand = 'TNF';
dose = '10ngmL';
order_name = 'MaxValue'
figure(i_sti)

cell_num=size(data_traj,2);
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

% subplot(1,length(vis_data_field),i_data_field)
h=heatmap(data_ordered,'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
%
XLabels = 0:5:((size(data_ordered,2)-1)*5);
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
saveas(gcf,strcat(savepath,'IntrinsicNoiseSim_',ligand,'_',dose,'_',order_name),'epsc');
close
