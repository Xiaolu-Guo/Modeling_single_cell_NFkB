function [] =  draw_sampling_20doses_codons_distrib()

sample_data = 1; % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
vis_data_field = {'pred_mode_amp'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};
vis_cv_field = {'vis_cv'};
fig_save_path = '../SubFigures2023/';
vers = '_Sampling';

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

        % 'Sim3_codon_r5_metric.mat'; Pam3CSK
        % 'Sim3_codon_r4_metric.mat'; PolyIC
        % 'Sim3_codon_r3_metric.mat'; CpG
        % 'Sim3_codon_r2_metric.mat'; LPS
        % 'Sim3_codon_r1_metric.mat'; TNF
        
data_filename_list = {'Sim3_codon_r5_metric.mat',... Pam3CSK
    'Sim3_codon_r4_metric.mat',... PolyIC
    'Sim3_codon_r3_metric.mat',... CpG
    'Sim3_codon_r2_metric.mat',... LPS
    'Sim3_codon_r1_metric.mat'}; % TNF

stim_list = {'Pam3CSK','PolyIC','CpG','LPS','TNF'};

for i_data_filename = 1:length(data_filename_list)
data_filename = data_filename_list{i_data_filename}

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
ylabel_list = {'Total','Duration','EvL','Speed','Peak','Osc'};


sti_all = cellfun(@strcat,collect_feature_vects.info_ligand,collect_feature_vects.info_dose_str,'UniformOutput', false) ;
ligand_inf = collect_feature_vects.info_ligand{1};

for i_codon =1:length(codon_list)
    
    figure(1)
    paperpos=[0,0,130*4,100]*1.5;
    papersize=[130*4 100]*1.5;
    draw_pos=[10,10,120*4+30,90]*1.5;
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    y = [];
    for i_sti = 1:length(sti_all)
        
        y = [y,collect_feature_vects.(codon_list{i_codon}){i_sti}];
    end
    z = y;
    % subplot(1,length(vis_data_field),i_data_field)
    al_goodplot(y,[],0.5,ones(size(y,2),1)*[ 0 0 0] ,'left',[],std(y(:))/2500);
    al_goodplot(z,[],0.5,ones(size(z,2),1)*[0 0 0],'right',[],std(y(:))/2500);
    % Unilateral plots for 2 timepoints (left: before, right: after), 3 groups.
    % One can produce multiple plots at once using a NxP input, P plots (1 per column).
    % One can use different options for each column.
    % If options are given only for 1 column, it is replicated for the others.
    %         xlim([0.4 8.6])
    %         xticks([1 2 3 4 5 6 7 8])
    %         xticklabels({'1hr', '2hr','3hr','4hr','5hr','6hr','7hr','8hr'})
    xlim([0.4 21.6])
    
    xticks(1:21)
    % xticklabels({'1','2','3','4','5'})
    %xtickangle(45)
    
    % ylabel(ylabel_list{i_codon})
    
    ylim([0,1]);
    set(gca,'fontsize',18,'fontname','Arial');
    saveas(gcf,strcat(fig_save_path,'codon_distrib_',stim_list{i_data_filename},'_',codon_list{i_codon},vers,'_0612'),'epsc');
    close
    
end

end

end