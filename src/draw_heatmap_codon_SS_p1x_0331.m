
load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/simulation_denoise/';%_fay_parameter/';
% vers = '_p01x_r1';
% load(strcat(data_save_file_path_1,'Sim5_SS_codon_metric',vers,'.mat'))

% vers = '_p25x_r1';
% load(strcat(data_save_file_path_1,'Sim5_SS_codon_metric',vers,'.mat'))


fig_save_path = '../SubFigures2023/';


    
vers_vec = {'_Sampling_WT','_p25x_r1'};%,'_p1x_r1','_p01x_r1'};
% vers_vec = {'_p01x_r1'};

%% draw heat map

for i_vers = 1:length(vers_vec)
    vers = vers_vec{i_vers};
if 1
    
    switch vers
         case  '_Sampling_WT'
            
            load(strcat(data_save_file_path_1,'Sim18_wt_IkBo_codon_metric1.mat'))
        case '_p25x_r1'
            load(strcat(data_save_file_path_1,'Sim5_SS_codon_metric',vers,'.mat'))
            
            
    end
    
    data_NFkB.info_ligand = data.info_ligand;
    data_NFkB.info_dose_str = data.info_dose_str;
    data_NFkB.info_dose_index = data.info_dose_index;
    for i_ligand = 1:length(data_NFkB.info_ligand)
        data_NFkB.info_num_cells{i_ligand} = 1000;
        data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
        
    end
    
    data_NFkB.exp = data_NFkB.model_sim;
    vis_data_field = {'model_sim'};
    % integrals 97 oscfreq 1 oscpower 1
    metric_order_name_{1} = 'integrals';
    metric_order_column = 97;
    
    order_name_vec = {'SS_TNF','SS_LPS','SS_CpG','SS_PolyIC','SS_Pam3CSK'};
    for i_order_name_vec = 1:length(order_name_vec)
        codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers,'_0331');
        
    end
    order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
    
    % 'sc_srs_' represents single-cell stumilus response specificity
    
    for i_order_i_ligand = 1% 1:length(data_NFkB.info_ligand)
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            [~,data_NFkB.order{i_ligand}] = sort(metrics{i_ligand}.(metric_order_name_{1})(1:9:end,metric_order_column),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
    end
    
end

%%  check codon distribution
if 1%
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    
    switch vers
        case  '_Sampling_WT'
            a = readmatrix('../raw_data2023/signaling_codon_WT_SS/X_codon_stim_WT_Sampling_r1_0331.csv');
        case '_p25x_r1'
            a = readmatrix('../raw_data2023/signaling_codon_WT_SS/X_codon_stim_SS_Sampling_p25x_r1_0331.csv');
             
    end
    
    for i_codon = 1:length(codon_list)
        figure(1)
        paperpos=[0,0,130,100]*1.5;
        papersize=[130 100]*1.5;
        draw_pos=[10,10,120,90]*1.5;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        
        y = [a(1:1000,i_codon),a(4001:5000,i_codon),a(2001:3000,i_codon),a(1001:2000,i_codon),a(3001:4000,i_codon)];
        z = y;
        % subplot(1,length(vis_data_field),i_data_field)
        al_goodplot(y,[],0.5,ones(size(y,2),1)*[ 0 0 0] ,'left',[],std(y(:))/2500);
        al_goodplot(z,[],0.5,ones(size(z,2),1)*[0 0 0],'right',[],std(y(:))/2500);

        xlim([0.4 5.6])
        
        xticks([1 2 3 4 5])
        xticklabels({'TNF','Pam3CSK','CpG','LPS','PolyIC'})
        xtickangle(45) 
        
        title(codon_list{i_codon})
        
        ylim([0,1]);
        set(gca,'fontsize',12,'fontname','Arial');
        saveas(gcf,strcat(fig_save_path,'codon_distrib_',codon_list{i_codon},vers,'_0331'),'epsc');
        close
    end
end

end