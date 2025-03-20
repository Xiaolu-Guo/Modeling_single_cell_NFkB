function [] = draw_all_ligand_sampling_codon_distrib_202306(data_save_file_path,fig_save_path)
% For Guo et al. Figure S4B-C, quality control for sampling parameter
% simulations
% 
% tested 05/12/2024, Matlab 2020a

fitting_dose ='alldose' ;%'lowdose','highdose','alldose','middose'
fig_opt.paper_opt.paperpos=[0,0,220,150]*1.5;
fig_opt.paper_opt.papersize=[220 150]*1.5;

load(strcat(data_save_file_path,'Sim2_fitting_',fitting_dose,'_codon_metric.mat'))

% 1. read data and metric
% asign data
metric_names = fieldnames(metrics{1});

for i_metric = 1:length(metrics)
    for i_metric_name = 1:length(metric_names)
        metrics_new{i_metric}.(metric_names{i_metric_name}) = metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
    end
end

i_ids = 1;
data_label = {'experiment','sampling'}; %,'sample'};

for i_data = 1:length(data.model_sim)
    for i_data_type = 1:2
        data_info.info_ligand{i_ids} = data.info_ligand{i_data};
        data_info.info_dose_str{i_ids} = data.info_dose_str{i_data};
        data_info.data_label{i_ids} = data_label{i_data_type};
        i_ids = i_ids+1;
        
    end
end

data_sim = data;
metrics_sim = metrics_new;

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

metrics_all = cell(1,length(metrics_sim)*2);

i_metric_index = 1;

metric_index_exp = [1:6,10:12,14:19]*2-1;
for i_metric = 1:length(metrics_sim)
    for i_metric_name = 1:length(metric_names)
        metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics{metric_index_exp(i_metric)}.(metric_names{i_metric_name})(:,:);
    end
    i_metric_index = i_metric_index +1;
    
    for i_metric_name = 1:length(metric_names)
        metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics_sim{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
    end
    i_metric_index = i_metric_index +1;
    
end

vis_data_field = {'experiment','sampling'}; %,'sample'};
[collect_feature_vects,~] = calculate_codon_from_metric2023(data_info,metrics_all); %,  parameter
collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

fieldlists = fieldnames(collect_feature_vects);
%TNF,Pam3CSK,CpG,LPS,PolyIC
index_high_dose = [6,5,30,29,18,17,12,11,24,23];

for i_field = 1:length(fieldlists)
    collect_feature_vects.(fieldlists{i_field}) = collect_feature_vects.(fieldlists{i_field})(index_high_dose);
end

fig_opt.save_file = strcat(fig_save_path,'Ligand_H_codon_exp_sampling_from_',fitting_dose,'_2023_06_distrib');

% violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
if 1 %run me for publish_figure
    %Figure S4B: tested 05/12/2024
    violin_plot_codon_sampling_2023_06(collect_feature_vects,fig_opt) % draw all codon and sti in one
end

if 1
    % fig_opt.save_file = strcat(fig_save_path,'All_ligand_fit_sampling_codon_2023_05_distrib');
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    sti_vec = cell(1,5);
    ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
    
    good_fit_pop_wdis = NaN(5,1);
    
    i_sti = 1;
    for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
        i_data_samp = i_ligand*2-1;
        i_data_exp = i_ligand*2;
        
        for i_codon = 1:length(codon_list)
            
            %  sti_vec{i_sti} = strcat(data.info_ligand{i_data},'-',data.info_dose_str{i_data});
            sti_vec{i_sti} = ligand_vec{i_ligand};
            exp_data = collect_feature_vects.(codon_list{i_codon}){i_data_exp};
            samp_data = collect_feature_vects.(codon_list{i_codon}){i_data_samp};

            
            W_diff(i_codon,i_sti) = w_distance(exp_data, samp_data, 2);
        end
        
        good_fit_pop_wdis(i_ligand,1) = mean(W_diff(:,i_sti));

        i_sti = i_sti+1;       
        
    end
    
    
    if 1 % W-distance
        %Figure S4C: tested 05/12/2024
        figure(1)
        paperpos = [0,0,200,120]/1.6;
        papersize = [200,120]/1.6;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(ligand_vec);
        c = reordercats(c,ligand_vec);
        bar_pct=good_fit_pop_wdis;
        b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
        ax2= gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        % ylabel({'discrepancy between' 'exp. and sampling'})
        XL = xlim();
        ylim([0,0.5]);
        %plot(XL,[0.1,0.1],'--r','LineWidth',1.5)
        % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)

        b(1).FaceColor = 'flat';
        
        TNF_color = [119 180 202]/255;
        LPS_color = [222 78 66]/255;
        CpG_color = [137 180 66]/255;
        PolyIC_cclor = [101 77 123]/255;
        P3CSK_color = [229 129 56]/255;
        b(1).CData = [TNF_color;LPS_color;CpG_color;PolyIC_cclor;P3CSK_color];        
        
        % set(gca,'fontsize',7,'FontName','Arial','YTick',[0,0.1,0.2,0.3,0.4,0.5],'YTickLabel',{},'XTickLabel',{})
        set(gca,'fontsize',7,'FontName','Arial','YTick',[0,0.1,0.2,0.3,0.4,0.5])

        saveas(gcf,strcat(fig_save_path,'Ligand_Hdose_fit_sample_',fitting_dose,'_bar_codon_good_fit_pop_wdis'),'epsc')
        close
    end
end



