function [] =  draw_traj_RMSD_distrib_2024_06(fig_save_path)
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
% ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
% dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
%     {'10ng/mL';'33ng/mL';'100ng/mL'};%;'33ng/mL';'100ng/mL'
%     {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
%     {'10ug/mL';'33ug/mL';'100ug/mL'};
%     {'10ng/mL';'100ng/mL';'1ug/mL'}};

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
        {'10ng/mL';'100ng/mL';'1ug/mL'};
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
     {'10ng/mL';'33ng/mL';'100ng/mL'};%;'33ng/mL';'100ng/mL'
    {'10ug/mL';'33ug/mL';'100ug/mL'}};

data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

RMSD_plot = [];
i_sti_index = 1;
for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
    
    % dose_index_vec = find(categorical(data.info_ligand)==ligand_vec{i_ligand});
    
    for i_dose = 1:length(dose_vec{i_ligand})
        i_sti = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        RMSD_vec{i_sti_index} = sqrt(sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/size(data.pred_mode_filter_nan{i_sti},2));
        i_sti_index = i_sti_index +1;
        % RMSD_plot = [RMSD_plot;RMSD_vec{i_sti_index}];
        % good_fit_ratio(i_ligand,i_dose) = sum(RMSD_plot <= thresh_good_fit)/length(RMSD_plot);
        
    end
    
end

if 1 % Figure S2A RMSD distribution: tested 05/15/2024
    figure(1)
    paperpos=[0,0,370,100]*1.5;
    papersize=[370 100]*1.5;
    draw_pos=[10,10,350,90]*1.5;
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    y = RMSD_vec';
    z = RMSD_vec';
    % subplot(1,length(vis_data_field),i_data_field)
    
    al_goodplot_pair_RMSD_diff_size(y,[],0.5,ones(length(y),1)*[0 0 255]/255 ,'left',[],std(cell2mat(y))/2500);
    al_goodplot_pair_RMSD_diff_size(z,[],0.5,ones(length(z),1)*[0 0 255]/255,'right',[],std(cell2mat(z))/2500);
    
    xlim([0.4 15.6])
    
    xticks([1:15])
    xticklabels({})
    %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
    
    ylim([0,0.06]);
%     for i_x = 1:15
%         plot([i_x,i_x],[0,5],'--','Color','k');hold on
%     end
    set(gca,'fontsize',14,'fontname','Arial');
    %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');
    
    saveas(gcf,strcat(fig_save_path,'RMSD_distrib_exp_sim'),'epsc');
    
    close
    
end



