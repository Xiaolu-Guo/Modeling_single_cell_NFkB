function [] =  draw_traj_RMSD_2023_05(fig_save_path)
% For Guo et al. Figure S2A example trajectories
% 
% tested 05/09/2024, Matlab 2020a


color_pct = {'r','g','b'};

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
load(strcat(data_save_file_path_1,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

rmsd_cas = 'RMSD';


ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'10ng/mL';'33ng/mL';'100ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);


plot_cas = 'percentage';
plot_cas = 'thresh';

good_fit_ratio = NaN(length(ligand_vec),3);
index = -ones(19,3);

for i_ligand = 4 % 1:length(ligand_vec)%TNF,LPS [27];%
    
    %     dose_index_vec = find(categorical(data.info_ligand)==ligand_vec{i_ligand});
    
    for i_dose = 2 % 1:length(dose_vec{i_ligand})
        i_sti = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        RMSD_vec{i_sti} = sqrt(sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2));
        % RMSD_vec{i_dose} = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2);
        
        nRMSD_vec{i_sti} = RMSD_vec{i_sti}./mean(data.exp_mode_filter_nan{i_sti},2);
        nRMSD_max_min_vec{i_sti} = RMSD_vec{i_sti}./(max(data.exp_mode_filter_nan{i_sti},[],2)-min(data.exp_mode_filter_nan{i_sti},[],2));
        mean_diff_vec{i_sti} = mean(abs(data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}),2);
        
        switch rmsd_cas
            case 'RMSD'
                RMSD_plot = [RMSD_vec{i_sti}];
                thresh_polyIC_med_dose_90pct = 3.92e-2;
                
            case 'nRMSD'
                RMSD_plot = [nRMSD_vec{i_sti}];
                thresh_polyIC_med_dose_90pct = 6.19e-1;
                
            case 'nRMSDr'
                RMSD_plot = [nRMSD_max_min_vec{i_sti} ];
                thresh_polyIC_med_dose_90pct = 2.20e-1;
            case 'mean_diff'
                RMSD_plot = [mean_diff_vec{i_sti} ];
                thresh_polyIC_med_dose_90pct = 3.09e-2;
                
        end
        
        thresh_polyIC_med_dose_90pct = 0.02;
%         thresh_polyIC_med_dose_90pct = 0.025;
%         thresh_polyIC_med_dose_90pct = 0.03;
%         thresh_polyIC_med_dose_90pct = 0.035;
%         thresh_polyIC_med_dose_90pct = 0.04;

        good_fit_ratio(i_ligand,i_dose) = sum(RMSD_plot <= thresh_polyIC_med_dose_90pct)/length(RMSD_plot);
       
        %; %;%;RMSD_vec{2};RMSD_vec{3}
        %         histogram(RMSD_TNF,'Normalization','cdf'); hold on
        
        [~,Ind] = sort(RMSD_plot, 'ascend');
        %         cdf_x = cumsum(col1); % Compute cdf
        %         cdf_x = cdf_x/cdf_x(end); % Normalize
        %         cdf = cumsum(Ind); % Compute cdf
        %         cdf = Ind/length(Ind); % Normalize
        % Find index of top 20 %
        pct = [0.2,0.5,0.8];
        pct = [0.1,0.5,0.9];
        thresh = [2e-2,3e-2,4e-2];% 2.5e-2,3.5e-2,
        
        %% plot traj
        if 1
            paperpos = [0,0,200,50]*1.8;
            papersize = [200,50]*1.8;
            draw_pos = [10,10,190,30]*1.8;
            
            figure(2)
            for i_plot =1:3
                
                subplot(1,3,i_plot)
                
                switch plot_cas
                    case 'percentage'
                        index(i_sti,i_plot) = Ind(ceil(pct(i_plot)*length(Ind))) ;%
                        plot(((1:length(data.pred_mode_filter_nan{i_sti}(index(i_sti,i_plot),:)))-1)/12,data.pred_mode_filter_nan{i_sti}(index(i_sti,i_plot),:),'r');hold on
                        plot(((1:length(data.exp_mode_filter_nan{i_sti}(index(i_sti,i_plot),:)))-1)/12,data.exp_mode_filter_nan{i_sti}(index(i_sti,i_plot),:),'k');hold on
                        % title(strcat(num2str(pct(i_plot)*100),'% nRMSDr=',num2str(prctile(RMSD_TNF,pct(i_plot)*100))),'Color',color_pct{i_plot})
                        title(strcat(num2str(pct(i_plot)*100),'% ',replace(rmsd_cas,'_',' '),'=',sprintf('%0.2e',prctile(RMSD_plot,pct(i_plot)*100))),'Color',color_pct{i_plot})
                        
                    case 'thresh'
                        
                        index(i_sti,i_plot) = find(RMSD_plot(Ind) <= thresh(i_plot), 1, 'last');
                        index(i_sti,i_plot) = Ind(index(i_sti,i_plot));
                        plot(((1:length(data.pred_mode_filter_nan{i_sti}(index(i_sti,i_plot),:)))-1)/12,data.pred_mode_filter_nan{i_sti}(index(i_sti,i_plot),:),'r');hold on
                        plot(((1:length(data.exp_mode_filter_nan{i_sti}(index(i_sti,i_plot),:)))-1)/12,data.exp_mode_filter_nan{i_sti}(index(i_sti,i_plot),:),'k');hold on
                        % title(strcat('RMSD=',num2str(thresh(i_plot))),'Color',color_pct{i_plot})
                        title(strcat(replace(rmsd_cas,'_',' '),'=',sprintf('%0.2e',RMSD_plot(index(i_sti,i_plot)))),'Color',color_pct{i_plot})
                        
                end
                
                xlabel('time')
                ylabel('NFkB')
                
                ylim([-0.05,0.3])
                xlim([0,8])
                set(gca,'fontsize',6,'fontweight','b')
                % data.info_ligand{i_sti}
                % data.info_dose_str{i_sti}
            end
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            saveas(gcf,strcat(fig_save_path,'trajectory_plots_',data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/','')),'epsc')
            close
        end
 
        
    end
  
end


