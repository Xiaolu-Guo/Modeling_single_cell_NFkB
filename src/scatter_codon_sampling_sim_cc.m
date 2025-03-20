
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

color_vec = {[0.7,0.7,0.7],[1,0,0],[0,1,0],[0,0,1]};

dose_all = {'Unstim','Low','Medium','High'};

load(strcat(data_save_file_path_1,data_filename))

metrics_new = cell(1,2);
metric_names = fieldnames(metrics{1});
for i_metric_name = 1:length(metric_names)
    for i_metric = 1:length(metrics)
        
        metrics_new{i_metric}.(metric_names{i_metric_name}) = metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
        data_info.info_ligand{i_metric} = collect_feature_vects.info_ligand{i_metric};
        data_info.info_dose_str{i_metric} = collect_feature_vects.info_dose_str{i_metric};
        data_info.data_label{i_metric} = 'sampling';
    end
end

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_new); %,  parameter
collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);


%% save the signaling coodn data to mutual info format

codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
ligand_inf = collect_feature_vects.info_ligand{1};


switch ligand_inf
    case 'TNF'
        index_3_doses =[1,7,11,15];% d7,d11,d15
    case 'CpG' %: d9,d11,d13
        index_3_doses = [1,9,11,13];
    case 'LPS' %: d11,d13,d15
        index_3_doses = [1,11,13,15];
    case 'Pam3CSK' %: d7,d11,d15
        index_3_doses = [1,7,11,15];
    case 'PolyIC' %: d11,d13,d15
        index_3_doses = [1,11,13,15];
end


for i_index_data = 1:length(index_3_doses)%[2,3]%
    num_cells = length(collect_feature_vects.Duration{index_3_doses(i_index_data)});
    x_val = -ones(36*num_cells,1);
    y_val = x_val;
    for i_codon =1:length(codon_list)
        for j_codon = 1:length(codon_list)
            x_start = (i_codon-1)*6*num_cells + (j_codon-1) * num_cells + 1;
            y_start = (j_codon-1)*6*num_cells + (i_codon-1) * num_cells + 1;
            if i_codon == j_codon
                          
                x_val(x_start:x_start+num_cells - 1 ) = 0;
                y_val(y_start:y_start+num_cells - 1)= 0;
                
            else
                x_start = (i_codon-1)*6*num_cells + (j_codon-1) * num_cells + 1;
                y_start = (j_codon-1)*6*num_cells + (i_codon-1) * num_cells + 1;
                
                x_val(x_start:x_start+num_cells - 1 ) = collect_feature_vects.(codon_list{i_codon}){index_3_doses(i_index_data)}+(i_codon-1)*1.2;
                y_val(y_start:y_start+num_cells - 1)= collect_feature_vects.(codon_list{i_codon}){index_3_doses(i_index_data)}+(i_codon-1)*1.2;
            end
        end
        
        codon_data = collect_feature_vects.(codon_list{i_codon}){index_3_doses(i_index_data)}+(i_codon-1)*1.2;
        pts = linspace(min(min(codon_data),(i_codon-1)*1.2),max(max(codon_data),i_codon*1.2),100);
        
        
        [f_codon(i_index_data,:,i_codon),xi_codon(i_index_data,:,i_codon)] = ksdensity(codon_data,pts,...
            'Function','pdf','Bandwidth',0.05);%,'Bandwidth',bw
        
        
    end
    
    h(i_index_data) = scatter(x_val,y_val,1,color_vec{i_index_data},'filled');hold on
    
end

scale_factor = max(f_codon,[],[1 2]);

for i_index_data = 1:length(index_3_doses)
    
    for i_codon = 1:length(codon_list)
        plot(xi_codon(i_index_data,:,i_codon),f_codon(i_index_data,:,i_codon)/scale_factor(i_codon)+(i_codon-1)*1.2,'Color',color_vec{i_index_data},'LineWidth',1); hold on
    end
    
end


for i_x = 0:6
    plot([i_x*1.2-0.1,i_x*1.2-0.1],[-0.2,6*1.2],'Color',[0,0,0]);hold on
    plot([-0.2,6*1.2],[i_x*1.2-0.1,i_x*1.2-0.1],'Color',[0,0,0]);hold on
    
end

for i_codon = 1:length(codon_list)
    
    text(i_codon*1.2-0.7,-0.5,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    
    %             if mod(i_codon,2)
    %            text(i_codon*1.2-0.7,-0.5,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    %             else
    %               text(i_codon*1.2-0.7,-0.9,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    %
    %             end
    text(-0.8,i_codon*1.2-0.6,codon_list{i_codon},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    
end
text(3,7.5,ligand_inf,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',9)
xlim([-2,8])
ylim([-2,8])
legend([h(1),h(2),h(3),h(4)],dose_all,'Location','northeast','box','off');



% medium size figure
paperpos = [0,0,100,100]*5;
papersize = [100,100]*5;
draw_pos = [20,20,60,60]*5;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
%to delete
fig_save_path = '../SubFigures2023/';

set(gca, 'fontsize',7,'XTick',{},'YTick',{})
saveas(gcf,strcat(fig_save_path,'codon_distrib_sample_',ligand_inf),'epsc')
close()

