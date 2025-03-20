function [] =  draw_traj_RMSD_2024_05(fig_save_path)
% For Guo et al. Figure S2, good fits ratio
% 
% tested 05/09/2024, Matlab 2020a

% change the threshold to get different figure S2A-B
thresh_good_fit = 0.02; 
% thresh_polyIC_med_dose_90pct = 0.03;
% thresh_polyIC_med_dose_90pct = 0.04;

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

good_fit_ratio = NaN(length(ligand_vec),3);

for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
    
    %     dose_index_vec = find(categorical(data.info_ligand)==ligand_vec{i_ligand});
    
    for i_dose = 1:length(dose_vec{i_ligand})
        i_sti = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        RMSD_vec{i_sti} = sqrt(sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2));
        
        RMSD_plot = [RMSD_vec{i_sti}];
        good_fit_ratio(i_ligand,i_dose) = sum(RMSD_plot <= thresh_good_fit)/length(RMSD_plot);
         
    end
  
end
%% figure S2B
if 1
    figure(1)
    paperpos = [0,0,120,50]*1.8;
    papersize = [120,50]*1.8;
    draw_pos = [10,10,100,30]*1.8;
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    index_ligand_reorder = [1,5,3,2,4];
    ligand_vec = ligand_vec(index_ligand_reorder);
    c = categorical(ligand_vec);
    c = reordercats(c,ligand_vec);
    bar_pct=good_fit_ratio(index_ligand_reorder,:)*100;
    b= bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);
    ax2= gca;
    ytickformat(ax2, '%g%%');
    ylim([0,100])
    ylabel({replace(strcat('fraction of',replace(rmsd_cas,'_',' ')),'of','of '), 'less than threshold'})
    % title(rmsd_cas)
    
    b(3).FaceColor = 'flat';
    b(2).FaceColor = 'flat';
    b(1).FaceColor = 'flat';
    
    TNF_color = [119 180 202]/255;
    LPS_color = [222 78 66]/255;
    CpG_color = [137 180 66]/255;
    PolyIC_cclor = [101 77 123]/255;
    P3CSK_color = [229 129 56]/255;
    b(3).CData = [TNF_color;P3CSK_color;CpG_color;LPS_color;PolyIC_cclor];
    b(2).CData = 1-(1-b(3).CData)*0.7;
    b(1).CData = 1-(1-b(3).CData)*0.3;
    
    set(gca,'fontsize',6,'fontweight','b')
    saveas(gcf,strcat(fig_save_path,'bar_good_fit_ratio_',replace(num2str(thresh_good_fit),'.','p')),'epsc')
    close
    
end

%% figure S2A
if 1
    i_sti =16;
    paperpos = [0,0,70,50]*1.8;
    papersize = [70,50]*1.8;
    draw_pos = [10,10,60,30]*1.8;
    
    RMSD_plot = [nRMSD_vec{i_sti}];
    [~,Ind] = sort(RMSD_plot, 'ascend');
  
    index_good = Ind(ceil(0.1*length(Ind))) ;%
    data_pred = data.pred_mode_filter_nan{i_sti};
    data_exp = data.exp_mode_filter_nan{i_sti};
    
    data_good_pred = data_pred(index_good,:);
    data_good_exp = data_exp(index_good,:);
    
    RMSD_good= sqrt(sum((data_good_pred - data_good_exp).^2,2)/size(data_good_exp,2))
    
    %nRMSD_good = RMSD_good{i_sti}./mean(data_good_exp,2)
    % nRMSD_max_good = RMSD_good{i_sti}./(max(data_good_exp,[],2)-min(data_good_exp,[],2))
    % mean_diff_good = mean(abs(data_good_pred - data_good_exp),2)
 
    time_vec = ((1:length(data_good_pred))-1)/12;
    
    figure(1)
    plot(time_vec,data_good_pred,'r');hold on
    plot(time_vec,data_good_exp,'k');hold on
    xlabel('time')
    ylabel('NFkB')
    
    ylim([-0.05,0.3])
    xlim([0,8])
    % title(strcat(replace(rmsd_cas,'_',' '),'=',sprintf('%0.2e',RMSD_plot(index_good))))
    set(gca,'fontsize',6,'fontweight','b')
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    saveas(gcf,strcat(fig_save_path,'PolyIC_good_fit'),'epsc')
    close
    
end


