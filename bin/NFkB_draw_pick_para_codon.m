function [] = NFkB_draw_pick_para_codon(data_info,figure_info)
% draw the highest dose
% can only calculate with one stimuli
% doses str is only working for one stimuli

% addpath('../CommonUsedFunction/')

if nargin <1
    Module_sens = 'TNF';
end


if exist(figure_info.save_codon_figure_path,'dir')
else
    mkdir(figure_info.save_codon_figure_path)
end

load(strcat(data_info.save_file_path,data_info.save_file_name),'sim_data_tbl')


length_para_vec = 21;

color_mapping = [linspace(1,0.5,floor(length_para_vec/2))',zeros(floor(length_para_vec/2),2);
     0,0,0;
    zeros(floor(length_para_vec/2),2),linspace(0.5,1,floor(length_para_vec/2))'];
Line_wid =ones(length_para_vec,1);
Line_wid(floor(length_para_vec/2)+1)=2;


% para_names = unique(sim_table_data.parameter_name);
dose_names = unique(sim_data_tbl.dose_str);
% species_names = unique(sim_table_data.species);
% para_names = {'params55'};%'params54';'params65';'params63'};

para_names = unique(sim_data_tbl.parameter_name);

species_names = {'nucNFkB'};


vis_data_field = {'model_sim'};%,'sample'};
data_label = {'simulation'};%,'sample'};

marker_str = {'^','s','o'};
%
% prod.TAK1 = TAK1;
% prod.IKK = IKK;
% prod.NFkB = NFkB;
%
% prod_vec = {'TAK1','IKK','NFkB'};

% currentfolder = pwd;


for i_para_names = 1:length(para_names)
    
    for i_dose_names = 1:length(dose_names)
        
        %             index = find(strcmp(sim_table_data.parameter_name,para_names{i_para_names}) &...
        %                 strcmp(sim_table_data.dose_str,dose_names{i_dose_names}) &...
        %                 strcmp(sim_table_data.species,species_names{i_species_names}));
        
        index_nuc_NFkB = sim_data_tbl.dose_str == dose_names(i_dose_names) ...
            & sim_data_tbl.species =='nucNFkB'...
            & sim_data_tbl.type == 'wt' ...
            & sim_data_tbl.parameter_name == para_names(i_para_names);
        
        
        
        % & sim_data_tbl.parameter_fold_change == parameter_fold_vec(i_parameter_fold_vec)
        
        % NFkBn + IkBaNFkBn
        data.model_sim{i_dose_names} =sim_data_tbl.trajectory(index_nuc_NFkB,1:5:end);
        data.info_ligand{i_dose_names} = sim_data_tbl.ligand(index_nuc_NFkB(1));
        data.info_dose_index{i_dose_names} = i_dose_names;
        data.info_dose_str{i_dose_names} = sim_data_tbl.dose_str(index_nuc_NFkB(1));
        data.info_num_cells{i_dose_names} = size(data.model_sim{i_dose_names},1);
        
        data.order{i_dose_names} = (1:data.info_num_cells{i_dose_names})';
        % data.time_points_num = 97
        
        
    end
    
    
    data.exp = data.model_sim;
    
    %cd(codon_path)
    
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
    
    i_species_names =1;
    
    % cd(currentfolder)
    if length(dose_names)==5
        dose_names_draw = dose_names([1,3,5]);
    else
        dose_names_draw = dose_names;
    end
    if isfield(figure_info,'codon_fields')
        for i_codon =1:length(figure_info.codon_fields)
            figure
            for i_dose = 1:length(dose_names_draw)
                index = sim_data_tbl.parameter_name == para_names(i_para_names) &...
                    sim_data_tbl.dose_str == dose_names_draw(i_dose) &...
                    sim_data_tbl.species == species_names{i_species_names};
                
                %                         index_nuc_NFkB = sim_data_tbl.dose_str == dose_names(i_dose_names) ...
                %             & sim_data_tbl.species =='nucNFkB'...
                %             & sim_data_tbl.type == 'wt' ...
                %             & sim_data_tbl.parameter_name == para_names(i_para_names);
                
                fold_change = sim_data_tbl.parameter_fold_change(index);
                
                (figure_info.codon_fields{i_codon})
                scatter(fold_change,collect_feature_vects.(figure_info.codon_fields{i_codon}){i_dose},8,color_mapping,'filled',marker_str{i_dose});hold on
                % title({para_names{i_para_names},figure_info.codon_fields{i_codon}})
                Set_figure_size_square_small
                
            end
            set(gca,'FontSize',8,'XScale','log','YLim',[-0.05,1.05])
            saveas(gcf,strcat(figure_info.save_codon_figure_path,'Codon_',figure_info.codon_fields{i_codon},para_names(i_para_names)),'epsc')
            
            close
        end
    end
    
    
    if isfield(figure_info,'metric_fields')
        for i_metric =1:length(figure_info.metric_fields)
            figure
            plot([1,1],data_info.y_metric_lim{i_metric},'Color',[0.5,0.5,0.5],'LineWidth',1.5);hold on

            for i_dose = 1:length(dose_names_draw)
                index = sim_data_tbl.parameter_name == para_names(i_para_names) &...
                    sim_data_tbl.dose_str == dose_names_draw(i_dose) &...
                    sim_data_tbl.species == species_names{i_species_names};
                
                fold_change = sim_data_tbl.parameter_fold_change(index);
                
                % (metrics{i_metric})
                figure_info.metric_fields{i_metric}
                max(metrics{i_dose}.(figure_info.metric_fields{i_metric}))
                min(metrics{i_dose}.(figure_info.metric_fields{i_metric}))
                
                metric_data = metrics{i_dose}.(figure_info.metric_fields{i_metric});
                if strcmp(figure_info.metric_fields{i_metric},'duration')
                    metric_data = metric_data(:,2);
                end
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
            
            ylim(data_info.y_metric_lim{i_metric})
            yt = yticks;
            ytlable = cell(1,length(yt));
            for i_yt = 1:length(yt)
                ytlable{i_yt} = '';
            end
            
            
            set(gca,'FontSize',8)
            set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)
            set(gca,'fontsize',7,'fontweight','b')
            
            set(gca,'FontSize',8,'XScale','log');%,'YLim',[-0.05,1.05]
            saveas(gcf,strcat(figure_info.save_codon_figure_path,'Codon_',figure_info.metric_fields{i_metric},string(para_names(i_para_names))),'epsc')
            
            close
        end
    end
    
    
end





