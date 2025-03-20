function [] =  draw_traj_RMSD_with_time_202401(fig_save_path)
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

for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
    
    %     dose_index_vec = find(categorical(data.info_ligand)==ligand_vec{i_ligand});
    
    for i_dose = 1:length(dose_vec{i_ligand})
        i_sti = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        RMSD_vec{i_sti} = (data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2;
        % RMSD_vec{i_dose} = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2);
        
        nRMSD_vec{i_sti} = RMSD_vec{i_sti}./(mean(data.exp_mode_filter_nan{i_sti},2)*ones(1,size(data.exp_mode_filter_nan{i_sti},2)));
        nRMSD_max_min_vec{i_sti} = RMSD_vec{i_sti}./((max(data.exp_mode_filter_nan{i_sti},[],2)-min(data.exp_mode_filter_nan{i_sti},[],2))*ones(1,size(data.exp_mode_filter_nan{i_sti},2)));
        mean_diff_vec{i_sti} = mean(abs(data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}),2);
        
        clear RMSD_accumulated nRMSD_accumulated nRMSD_accumulated_bycell
        %         for i_hour = 1:8
        %             RMSD_accumulated(:,i_hour) = sum(RMSD_vec{i_sti}(:,(2:13)+i_hour*12-12),2);
        %             nRMSD_accumulated(:,i_hour) = RMSD_accumulated(:,i_hour)/mean(mean(data.exp_mode_filter_nan{i_sti}(:,(2:13)+i_hour*12-12),2));
        %             nRMSD_accumulated_bycell(:,i_hour) = RMSD_accumulated(:,i_hour)./mean(data.exp_mode_filter_nan{i_sti}(:,(2:13)+i_hour*12-12),2);
        %
        %         end
        
        for i_hour = 1:4
            RMSD_accumulated(:,i_hour) = sqrt(sum(RMSD_vec{i_sti}(:,(2:25)+i_hour*24-24),2)/24);
            nRMSD_accumulated(:,i_hour) = RMSD_accumulated(:,i_hour)/mean(mean(data.exp_mode_filter_nan{i_sti}(:,(2:25)+i_hour*24-24),2));
            nRMSD_accumulated_bycell(:,i_hour) = RMSD_accumulated(:,i_hour)./mean(data.exp_mode_filter_nan{i_sti}(:,(2:25)+i_hour*24-24),2);
            
        end
        figure(1)
        paperpos=[0,0,130,100]*3;
        papersize=[130 100]*3;
        draw_pos=[10,10,120,90]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        y = nRMSD_accumulated;
        z = nRMSD_accumulated;
        % subplot(1,length(vis_data_field),i_data_field)
        al_goodplot(y,[],0.5,ones(size(RMSD_accumulated,2),1)*[ 0 0 0] ,'left',[],std(y(:))/2500);
        al_goodplot(z,[],0.5,ones(size(RMSD_accumulated,2),1)*[0 0 0],'right',[],std(y(:))/2500);
        % Unilateral plots for 2 timepoints (left: before, right: after), 3 groups.
        % One can produce multiple plots at once using a NxP input, P plots (1 per column).
        % One can use different options for each column.
        % If options are given only for 1 column, it is replicated for the others.
%         xlim([0.4 8.6])
%         xticks([1 2 3 4 5 6 7 8])
%         xticklabels({'1hr', '2hr','3hr','4hr','5hr','6hr','7hr','8hr'})
        xlim([0.4 4.6])

        xticks([1 2 3 4])
        xticklabels({'2hr','4hr','6hr','8hr'})


        title(data.info_ligand{i_sti})
        
        ylim([0,1]);
        set(gca,'fontsize',14,'fontname','Arial');
        saveas(gcf,strcat(fig_save_path,'RMSD_timept_distrib_',ligand_vec{i_ligand},'_',num2str(i_dose)),'epsc');
        close
        
    end
    
end


