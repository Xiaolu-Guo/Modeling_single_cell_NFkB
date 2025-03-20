function [] =  draw_traj_RMSD_2023(data_save_file_path,fig_save_path)
color_pct = {'r','g','b'};

for data_proj_num_vec=[2:6]%TNF,LPS [27];%
    data = read_monolix_data_to_matrix_2023(data_proj_num_vec,data_save_file_path);
    

    
    for i_sti = 1:length(data.pred_mode)
        
        RMSD_vec{i_sti} = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/length(data.pred_mode{i_sti});
        
        nRMSD_vec{i_sti} = RMSD_vec{i_sti}./mean(data.exp_mode_filter_nan{i_sti},2);
        nRMSD_max_min_vec{i_sti} = RMSD_vec{i_sti}./(max(data.exp_mode_filter_nan{i_sti},[],2)-min(data.exp_mode_filter_nan{i_sti},[],2));
        
        RMSD_TNF = [nRMSD_max_min_vec{i_sti} ];%[nRMSD_vec{i_sti} ]; %[RMSD_vec{i_sti}];%;RMSD_vec{2};RMSD_vec{3}
        %         histogram(RMSD_TNF,'Normalization','cdf'); hold on
        
        [col1,Ind] = sort(RMSD_TNF, 'ascend');
        cdf = cumsum(col1); % Compute cdf
        cdf = cdf/cdf(end); % Normalize
        % Find index of top 20 %
        pct = [0.2,0.5,0.8];
        pct = [0.01,0.5,0.99];
        thresh = [1e-4,2e-4,3e-4];
        
        %% plot traj
            paperpos = [0,0,200,50]*1.8;
    papersize = [200,50]*1.8;
    draw_pos = [10,10,190,30]*1.8;
    
        figure(2)
        for i_plot =1:3
            
            subplot(1,3,i_plot)
            if 1
                index = find(cdf >= pct(i_plot), 1, 'first');
                plot(((1:length(data.pred_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.pred_mode_filter_nan{i_sti}(Ind(index),:),'r');hold on
                plot(((1:length(data.exp_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.exp_mode_filter_nan{i_sti}(Ind(index),:),'k');hold on
                % title(strcat(num2str(pct(i_plot)*100),'% nRMSDr=',num2str(prctile(RMSD_TNF,pct(i_plot)*100))),'Color',color_pct{i_plot})
                title(strcat(num2str(pct(i_plot)*100),'% nRMSDr=',sprintf('%0.2e',prctile(RMSD_TNF,pct(i_plot)*100))),'Color',color_pct{i_plot})

            else
                index = find(RMSD_TNF >= thresh(i_plot), 1, 'first');
                plot(((1:length(data.pred_mode_filter_nan{i_sti}(index,:)))-1)/12,data.pred_mode_filter_nan{i_sti}(index,:),'r');hold on
                plot(((1:length(data.exp_mode_filter_nan{i_sti}(index,:)))-1)/12,data.exp_mode_filter_nan{i_sti}(index,:),'k');hold on
%                 title(strcat('RMSD=',num2str(thresh(i_plot))),'Color',color_pct{i_plot})
                title(strcat('RMSD=',sprintf('%0.2e',thresh(i_plot))),'Color',color_pct{i_plot})


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

        
        %% plot cdf
        
                
        paperpos = [0,0,60,50]*1.8;
        papersize = [60,50]*1.8;
        draw_pos = [10,10,40,30]*1.8;
        
        figure(1)
        % histogram(RMSD_TNF,'Normalization','cdf'); hold on
        h = cdfplot(RMSD_TNF); hold on
        h.Color = 'k';
        h.LineWidth = 1;
        title('');
        xlabel('nRMSDr')
        ylabel('cdf')
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        set(gca,'fontsize',6,'fontweight','b')
        if 1
            for i_pct = 1:length(pct)
                plot([prctile(RMSD_TNF,pct(i_pct)*100),prctile(RMSD_TNF,pct(i_pct)*100)],[0,1],'Color',color_pct{i_pct});hold on
            end
        else
            for i_thresh = 1:length( thresh)
                plot([thresh(i_thresh),thresh(i_thresh)],[0,1],'Color',color_pct{i_thresh});hold on
            end
        end
        % legend(data.info_dose_str{end})
        saveas(gcf,strcat(fig_save_path,'traj_RMSD_',data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/','')),'epsc')
        close
        
    end
    
    
end

