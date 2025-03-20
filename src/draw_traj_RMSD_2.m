function [] =  draw_traj_RMSD_2(data_save_file_path,fig_save_path)
color_pct = {'r','g','b'};

for data_proj_num_vec=[3:7]%TNF,LPS [27];%
    data = read_monolix_data_to_matrix_2(data_proj_num_vec,data_save_file_path);
   
    paperpos = [0,0,200,50]*1.8;
    papersize = [200,50]*1.8;
    draw_pos = [10,10,190,30]*1.8;
    
    for i_sti = 1:length(data.pred_mode)
        
        RMSD_vec{i_sti} = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/length(data.pred_mode{i_sti});
        
        nRMSD_vec{i_sti} = RMSD_vec{i_sti}./mean(data.exp_mode_filter_nan{i_sti},2);

        figure(1)
        
        RMSD_TNF = [RMSD_vec{i_sti}];%;RMSD_vec{2};RMSD_vec{3}
%         histogram(RMSD_TNF,'Normalization','cdf'); hold on
      
        [col1,Ind] = sort(RMSD_TNF, 'ascend');
        cdf = cumsum(col1); % Compute cdf
        cdf = cdf/cdf(end); % Normalize
        % Find index of top 20 %
        pct = [0.2,0.5,0.8];
        thresh = [1e-4,2e-4,3e-4];
        
        figure(2)
        for i_plot =1:3
            
            subplot(1,3,i_plot)
            index = find(cdf >= pct(i_plot), 1, 'first');
            plot(((1:length(data.pred_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.pred_mode_filter_nan{i_sti}(Ind(index),:),'r');hold on
            plot(((1:length(data.exp_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.exp_mode_filter_nan{i_sti}(Ind(index),:),'k');hold on
            %         index = find(RMSD_TNF >= thresh(i_plot), 1, 'first');
            %         plot(((1:length(data.pred_mode_filter_nan{i_sti}(index,:)))-1)/12,data.pred_mode_filter_nan{i_sti}(index,:),'r');hold on
            %         plot(((1:length(data.exp_mode_filter_nan{i_sti}(index,:)))-1)/12,data.exp_mode_filter_nan{i_sti}(index,:),'k');hold on
            
            xlabel('time')
            ylabel('NFkB')
            title(strcat(num2str(pct(i_plot)*100),'%'),'Color',color_pct{i_plot})
            
            ylim([-0.05,0.25])
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
        
    paperpos = [0,0,60,50]*1.8;
    papersize = [60,50]*1.8;
    draw_pos = [10,10,40,30]*1.8;
    
    figure(1)
    % histogram(RMSD_TNF,'Normalization','cdf'); hold on
    h = cdfplot(RMSD_TNF); hold on
    h.Color = 'k';
    h.LineWidth = 1;
    title('');
    xlabel('RMSD')
    ylabel('cdf')
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    set(gca,'fontsize',6,'fontweight','b')
    for i_pct = 1:length(pct)
        plot([prctile(RMSD_TNF,pct(i_pct)*100),prctile(RMSD_TNF,pct(i_pct)*100)],[0,1],'Color',color_pct{i_pct});hold on
    end
        % legend(data.info_dose_str{end})

    saveas(gcf,strcat(fig_save_path,'traj_RMSD_',data.info_ligand{i_sti}),'epsc')
    close
end

