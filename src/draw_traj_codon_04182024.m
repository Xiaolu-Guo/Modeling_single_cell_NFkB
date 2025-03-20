function [] = draw_traj_codon_04182024(data_save_file_path,fig_save_path,module_info)


vers = '_04172024';


features_metric_fields = {'duration','oscpower','max_pos_integral','pos_pk1_time','max_value','time2HalfMaxPosIntegral','time2HalfMaxValue'};%'max_pos_pk1_speed',,'pos_pk1_amp','time2HalfMaxPosIntegral'

data_info.save_file_path = data_save_file_path;%  '../../NFkB_para_estm_project/NFkB_figures/ParameterScan/data/';
data_info.species_outputname = {'nucNFkB'};%;'TNFR';'IKK'
data_info.species_composition = {{'NFkBn';'IkBaNFkBn'}}; % must ;{'TNFR'};{'IKK'}
% be r x 1, for each cell i must be ri x 1
data_info.flag = '';
data_info.type = 'wt';

if isfield(module_info,'module')
    data_info.save_file_name = strcat('Module_Sens_',module_info.module{1},'_',module_info.ligand{1},vers);
else
    data_info.save_file_name = strcat('Module_Sens_',module_info.ligand{1},'_',module_info.para_type{1},vers);
end

load(strcat(data_info.save_file_path,data_info.save_file_name));

%a = 1;

%% draw traj
if 1
    
    high_dose = max(sim_data_tbl.dose_val);
    %% get the data for drawing/plotting
    
    % parameter_name_vec = unique(sim_data_tbl.parameter_name);
    switch module_info.para_type{1}
        case 'NFkBtot'
            parameter_fold_vec = 10.^linspace(-1,1,21);
            
        otherwise
            parameter_fold_vec = unique(sim_data_tbl.parameter_fold_change);
    end
    color_mapping = [ linspace(1,0.5,floor(length(parameter_fold_vec)/2))',zeros(floor(length(parameter_fold_vec)/2),2);
        0,0,0;
        zeros(floor(length(parameter_fold_vec)/2),2),linspace(0.5,1,floor(length(parameter_fold_vec)/2))'];
    Line_wid = 0.75 * ones(length(parameter_fold_vec),1);
    Line_wid(floor(length(parameter_fold_vec)/2)+1) = 1.5;
    
    %for i_parameter_name_vec = 1:length(parameter_name_vec)
    figure(1)
    
    switch module_info.para_type{1}
        case 'NFkBtot'
            index_nuc_NFkB = find(sim_data_tbl.dose_val == high_dose ...
                & sim_data_tbl.species =='nucNFkB'...
                & sim_data_tbl.type == 'wt') ;
            
            for i_parameter_fold_vec = 1:length(parameter_fold_vec)
                
                
                nucNFkB_traj = sim_data_tbl.trajectory(index_nuc_NFkB(i_parameter_fold_vec),:);
                
                plot(1:length(nucNFkB_traj),nucNFkB_traj,'LineWidth',Line_wid(i_parameter_fold_vec),'Color',color_mapping(i_parameter_fold_vec,:));hold on
                
            end
                            ylim([0,3])

        otherwise
            for i_parameter_fold_vec = 1:length(parameter_fold_vec)
                index_nuc_NFkB = sim_data_tbl.dose_val == high_dose ...
                    & sim_data_tbl.species =='nucNFkB'...
                    & sim_data_tbl.type == 'wt' ... & sim_data_tbl.parameter_name == parameter_name_vec(i_parameter_name_vec)...
                    & sim_data_tbl.parameter_fold_change(:,1) == parameter_fold_vec(i_parameter_fold_vec);
                
                nucNFkB_traj = sim_data_tbl.trajectory(index_nuc_NFkB,:);
                
                plot(1:length(nucNFkB_traj),nucNFkB_traj,'LineWidth',Line_wid(i_parameter_fold_vec),'Color',color_mapping(i_parameter_fold_vec,:));hold on
                
            end
            
                ylim([0,0.3])

    end
    
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
    
    saveas(gcf,strcat(fig_save_path,'Sens_traj',vers,'_',module_info.ligand{1},'_',module_info.para_type{1}),'epsc')
    %saveas(gcf,strcat(fig_save_path,fig_save_name,string(parameter_name_vec(i_parameter_name_vec))),'svg')
    close
    %end
end

%% draw codon
switch module_info.para_type{1}
    case 'NFkBtot'
        metric_y_lim = {[0,8],[0,1e-3],[0,5],[0,1],[0,3],[0,4],[0,1]};
        
    otherwise
        metric_y_lim = {[0,8],[0,1e-3],[0,2],[0,4],[0,0.3],[0,6],[0,2]};
end
dose_names = unique(sim_data_tbl.dose_str);
marker_str = {'^','s','o'};

if 1
    for i_metric =1:length(features_metric_fields)
        figure
        plot([1,1],metric_y_lim{i_metric},'Color',[0.5,0.5,0.5],'LineWidth',1.5);hold on
        
        for i_dose = 1:length(dose_names)
            index = sim_data_tbl.dose_str == dose_names(i_dose);% &... sim_data_tbl.parameter_name == para_names(i_para_names) &...
            
            
            switch module_info.para_type{1}
                case 'NFkBtot'
                    fold_change = 10.^linspace(-1,1,21);
                    
                otherwise
                    fold_change = sim_data_tbl.parameter_fold_change(index);
            end
            
            
            % (metrics{i_metric})
            
            metric_data = sim_data_tbl.(features_metric_fields{i_metric})(index,:);
            if strcmp(features_metric_fields{i_metric},'duration')
                metric_data = metric_data(:,2);
            end
            
            features_metric_fields{i_metric}
            
            max(metric_data)
            min(metric_data)
            
            scatter(fold_change,metric_data,8,color_mapping,'filled',marker_str{i_dose});hold on
            % title({para_names{i_para_names},figure_info.metric_fields{i_metric}})
            Set_figure_size_square_small_2
            
        end
        
        XL = [0.1,10];%xlim();
        xlim(XL);
        xt = get(gca,'XTick');
        xtlable = cell(1,length(xt));
        for i_xt = 1:length(xt)
            xtlable{i_xt} = '';
        end
        
        ylim(metric_y_lim{i_metric})
        yt = yticks;
        ytlable = cell(1,length(yt));
        for i_yt = 1:length(yt)
            ytlable{i_yt} = '';
        end
        
        
        set(gca,'FontSize',8)
        set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)
        set(gca,'fontsize',7,'fontweight','b')
        
        set(gca,'FontSize',8,'XScale','log');%,'YLim',[-0.05,1.05]
        saveas(gcf,strcat(fig_save_path,'Codon',vers,'_',module_info.ligand{1},'_',module_info.para_type{1},'_',features_metric_fields{i_metric}),'epsc')
        
        close
    end
end