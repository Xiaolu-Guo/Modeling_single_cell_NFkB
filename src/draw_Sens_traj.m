% 'Module_Sens_TNF.mat'
% run Module_parameters_sim.m
function [] = draw_Sens_traj(data_save_file_path,fig_save_path)

% fig_save_path = '../../TNF_Knock_out/Figures/';
% fig_save_name = 'TNF_sens_para_nucNFkB_';

% save_file_name = strcat(save_file_name,'.mat');
module_all =  {'core','TNF','CpG','LPS','Pam3CSK','polyIC'};
for i_ligand = 1:length(module_all)
    module_name = module_all{i_ligand};
    clear sim_data_tbl
load(strcat(data_save_file_path,'Module_Sens_',module_name,'.mat'),'sim_data_tbl');

high_dose = max(sim_data_tbl.dose_val);
%% get the data for drawing/plotting

parameter_name_vec = unique(sim_data_tbl.parameter_name);
parameter_fold_vec = unique(sim_data_tbl.parameter_fold_change);

color_mapping = [ linspace(1,0.5,floor(length(parameter_fold_vec)/2))',zeros(floor(length(parameter_fold_vec)/2),2);
       0,0,0;
    zeros(floor(length(parameter_fold_vec)/2),2),linspace(0.5,1,floor(length(parameter_fold_vec)/2))'];
Line_wid = 0.75 * ones(length(parameter_fold_vec),1);
Line_wid(floor(length(parameter_fold_vec)/2)+1) = 1.5;

for i_parameter_name_vec = 1:length(parameter_name_vec)
    figure(i_parameter_name_vec)
    
    for i_parameter_fold_vec = 1:length(parameter_fold_vec)
        index_nuc_NFkB = sim_data_tbl.dose_val == high_dose ...
            & sim_data_tbl.species =='nucNFkB'...
            & sim_data_tbl.type == 'wt' ...
            & sim_data_tbl.parameter_name == parameter_name_vec(i_parameter_name_vec)...
            & sim_data_tbl.parameter_fold_change == parameter_fold_vec(i_parameter_fold_vec);
        
        nucNFkB_traj = sim_data_tbl.trajectory(index_nuc_NFkB,:);
        
        plot(1:length(nucNFkB_traj),nucNFkB_traj,'LineWidth',Line_wid(i_parameter_fold_vec),'Color',color_mapping(i_parameter_fold_vec,:));hold on
        
    end
    
    
    ylim([0,0.3])
    yt = yticks;
    ytlable = cell(1,length(yt));
    for i_yt = 1:length(yt)
        ytlable{i_yt} = '';
    end
    
    xlim([0,96*5])
    xt = 0:(24*5):(96*5);
    xtlable = cell(1,length(xt));
    for i_xt = 1:length(xt)
        xtlable{i_xt} = '';
    end
    
    set(gca,'FontSize',8)
    set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)
    set(gca,'fontsize',7,'fontweight','b')
    Set_figure_size_square_small

    saveas(gcf,strcat(fig_save_path,'Sens_traj_',module_name,'_',string(parameter_name_vec(i_parameter_name_vec))),'epsc')
    %saveas(gcf,strcat(fig_save_path,fig_save_name,string(parameter_name_vec(i_parameter_name_vec))),'svg')
    close
    
end

end
end


