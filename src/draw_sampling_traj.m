
function [] = draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path )
data_info.save_file_path = data_save_file_path;

data_info.save_file_name = save_filename; % only beginning, no .mat
gene_info.parameter_name_vec = {{'params66','params99','params100','params101','params54','params53'}};
fig_name = 'para_sample_';


for i_para_name = 1:length(gene_info.parameter_name_vec)
    
    %% load sim data
    sim_info.parameter_name = gene_info.parameter_name_vec{i_para_name};
    
    if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
        save_file_name = data_info.save_file_name;
        fig_save_name = strcat(fig_name,vers_fig);
        for i_para_name = 1:length(sim_info.parameter_name)
            save_file_name = strcat(save_file_name,'_',sim_info.parameter_name{i_para_name});
            fig_save_name = strcat(fig_save_name,'_',sim_info.parameter_name{i_para_name});
        end
        % save_file_name = strcat(save_file_name,'.mat');
        load(strcat(data_info.save_file_path,save_file_name,'.mat'),'sim_data_tbl','metrics','data');
        
    end
    
    %% get the data for drawing/plotting
    index_nuc_NFkB_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'wt';
    
    
    traj_wt = sim_data_tbl.trajectory(index_nuc_NFkB_wt,:);
    [~,traj_order] = sort(max(traj_wt,[],2),'descend');
    traj_wt = traj_wt(traj_order,:);
    traj_wt = traj_wt(:,1:5:end);
    
    
    
    %% draw the picture and save
    
    
    if 1 %
        draw_traj_heatmap(1,traj_wt)
        
        figure(1)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_wt'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_wt'),'svg')
        close
        
    end
    
    
    
end
end



%% draw heatmap
% to do: delete ylabel, xlabel only time 0,4,8
% proper figure size
% seperate the figures
% draw heatmap
function [] = draw_traj_heatmap(fig_num,traj_wt)

figure(fig_num(1))
h1=heatmap(traj_wt,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25]);

XLabels = 0:1:size(traj_wt,2)-1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels/12);
% Replace all but the fifth elements by spaces
CustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h1.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_wt,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h1.YDisplayLabels = YCustomXLabels;
colorbar('off')
Set_figure_size_square% clb=colorbar;
%
% clb.Label.String = 'NFkB(S.I.)';


end


