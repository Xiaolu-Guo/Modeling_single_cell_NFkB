function [] =  draw_traj_RMSD_2023_03(data_save_file_path,fig_save_path)
color_pct = {'r','g','b'};

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
load(strcat(data_save_file_path_1,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

rmsd_cas = 'nRMSD';


ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);


plot_cas = 'percentage';
good_fit_ratio = NaN(length(ligand_vec),3);

for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
    
    %     dose_index_vec = find(categorical(data.info_ligand)==ligand_vec{i_ligand});
    
    for i_dose = 1:length(dose_vec{i_ligand})
        i_sti = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        RMSD_vec{i_dose} = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2);
        nRMSD_vec{i_dose} = RMSD_vec{i_dose}./mean(data.exp_mode_filter_nan{i_sti},2);
        nRMSD_max_min_vec{i_dose} = RMSD_vec{i_dose}./(max(data.exp_mode_filter_nan{i_sti},[],2)-min(data.exp_mode_filter_nan{i_sti},[],2));
        
        
        switch rmsd_cas
            case 'RMSD'
                RMSD_plot = [nRMSD_vec{i_dose}];
                thresh_polyIC_med_dose_90pct = 1.98e-2;
                
            case 'nRMSD'
                RMSD_plot = [RMSD_vec{i_dose}];
                thresh_polyIC_med_dose_90pct = 1.54e-3;
                
            case 'nRMSDr'
                RMSD_plot = [nRMSD_max_min_vec{i_dose} ];
                thresh_polyIC_med_dose_90pct = 7.68e-3;
                
        end
        good_fit_ratio(i_ligand,i_dose) = sum(RMSD_plot <= thresh_polyIC_med_dose_90pct)/length(RMSD_plot);
        
        
        %; %;%;RMSD_vec{2};RMSD_vec{3}
        %         histogram(RMSD_TNF,'Normalization','cdf'); hold on
        
        [col1,Ind] = sort(RMSD_plot, 'ascend');
        cdf = cumsum(col1); % Compute cdf
        cdf = cdf/cdf(end); % Normalize
        % Find index of top 20 %
        pct = [0.2,0.5,0.8];
        pct = [0.1,0.5,0.9];
        % thresh = [1e-4,2e-4,3e-4];
        
        %% plot traj
        if 0
            paperpos = [0,0,200,50]*1.8;
            papersize = [200,50]*1.8;
            draw_pos = [10,10,190,30]*1.8;
            
            figure(2)
            for i_plot =1:3
                
                subplot(1,3,i_plot)
                
                switch plot_cas
                    case 'percentage'
                        index = find(cdf >= pct(i_plot), 1, 'first');
                        plot(((1:length(data.pred_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.pred_mode_filter_nan{i_sti}(Ind(index),:),'r');hold on
                        plot(((1:length(data.exp_mode_filter_nan{i_sti}(Ind(index),:)))-1)/12,data.exp_mode_filter_nan{i_sti}(Ind(index),:),'k');hold on
                        % title(strcat(num2str(pct(i_plot)*100),'% nRMSDr=',num2str(prctile(RMSD_TNF,pct(i_plot)*100))),'Color',color_pct{i_plot})
                        title(strcat(num2str(pct(i_plot)*100),'% ',rmsd_cas,'=',sprintf('%0.2e',prctile(RMSD_plot,pct(i_plot)*100))),'Color',color_pct{i_plot})
                        
                    case 'thresh'
                        index = find(RMSD_plot >= thresh(i_plot), 1, 'first');
                        plot(((1:length(data.pred_mode_filter_nan{i_sti}(index,:)))-1)/12,data.pred_mode_filter_nan{i_sti}(index,:),'r');hold on
                        plot(((1:length(data.exp_mode_filter_nan{i_sti}(index,:)))-1)/12,data.exp_mode_filter_nan{i_sti}(index,:),'k');hold on
                        %                 title(strcat('RMSD=',num2str(thresh(i_plot))),'Color',color_pct{i_plot})
                        title(strcat(rmsd_cas,'=',sprintf('%0.2e',thresh(i_plot))),'Color',color_pct{i_plot})
                        
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
        
        
        %% plot cdf
        
        if 0
            paperpos = [0,0,60,50]*1.8;
            papersize = [60,50]*1.8;
            draw_pos = [10,10,40,30]*1.8;
            
            figure(1)
            % histogram(RMSD_TNF,'Normalization','cdf'); hold on
            h = cdfplot(RMSD_plot); hold on
            h.Color = 'k';
            h.LineWidth = 1;
            title('');
            xlabel(rmsd_cas)
            ylabel('cdf')
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            set(gca,'fontsize',6,'fontweight','b')
            switch plot_cas
                case 'percentage'
                    for i_pct = 1:length(pct)
                        plot([prctile(RMSD_plot,pct(i_pct)*100),prctile(RMSD_plot,pct(i_pct)*100)],[0,1],'Color',color_pct{i_pct});hold on
                    end
                case 'thresh'
                    for i_thresh = 1:length( thresh)
                        plot([thresh(i_thresh),thresh(i_thresh)],[0,1],'Color',color_pct{i_thresh});hold on
                    end
            end
            % legend(data.info_dose_str{end})
            saveas(gcf,strcat(fig_save_path,'traj_RMSD_',data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/','')),'epsc')
            close
        end
        
    end
    
    
    
end
%% figure 2B
if 1
figure(1)
paperpos = [0,0,120,50]*1.8;
papersize = [120,50]*1.8;
draw_pos = [10,10,100,30]*1.8;

set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
c = categorical(ligand_vec);
c = reordercats(c,ligand_vec);
bar_pct=good_fit_ratio*100;
b= bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);
ax2= gca;
ytickformat(ax2, '%g%%');
ylim([75,100])
ylabel({'fraction of', 'good fitting'})
% title(rmsd_cas)

b(3).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).FaceColor = 'flat';

TNF_color = [119 180 202]/255;
LPS_color = [222 78 66]/255;
CpG_color = [137 180 66]/255;
PolyIC_cclor = [101 77 123]/255;
P3CSK_color = [229 129 56]/255;
b(3).CData = [TNF_color;LPS_color;CpG_color;PolyIC_cclor;P3CSK_color];
b(2).CData = 1-(1-b(3).CData)*0.7;
b(1).CData = 1-(1-b(3).CData)*0.3;

set(gca,'fontsize',6,'fontweight','b')
saveas(gcf,strcat(fig_save_path,'bar_good_fit_ratio_',rmsd_cas),'epsc')
close

end


