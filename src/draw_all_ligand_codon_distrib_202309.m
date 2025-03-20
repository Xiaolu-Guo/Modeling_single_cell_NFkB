function [] = draw_all_ligand_codon_distrib_202309(data_save_file_path,fig_save_path)
% For Guo et al. Figure S2CDE, signaling codon distribution violin plot,
% Wasserstein dist
% 
% tested 05/10/2024, Matlab 2020a


fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
fig_opt.paper_opt.papersize=[220 180]*3;

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

%% codon distribution
% codon and sti seperately, and save
fig_opt.save_file = strcat(fig_save_path,'All_ligand_codon_2023_06_distrib');

codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
sti_vec = cell(1,15);
ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
dose_label = {'L','M','H'};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

good_fit_pop = NaN(5,3);

i_sti = 1;
for i_ligand = 1:length(ligand_vec)%TNF,LPS [27];%
    
    for i_dose = 1:length(dose_vec{i_ligand})
        i_data = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
        
        for i_codon = 1:length(codon_list)
            
            %  sti_vec{i_sti} = strcat(data.info_ligand{i_data},'-',data.info_dose_str{i_data});
            sti_vec{i_sti} = strcat(data.info_ligand{i_data},'-',dose_label{i_dose});
            exp_data = collect_feature_vects.(codon_list{i_codon}){i_data*2-1};
            sim_data = collect_feature_vects.(codon_list{i_codon}){i_data*2};
            pts = linspace(min(min(exp_data),min(sim_data)),max(max(exp_data),max(sim_data)),50);
            
            [~,~,bw_exp] = ksdensity(exp_data,pts,...
                'Function','pdf');
            [~,~,bw_sim] = ksdensity(sim_data,pts,...
                'Function','pdf');%
            bw = min(bw_exp,bw_sim);
            
            [f_exp,xi_exp] = ksdensity(exp_data,pts,...
                'Function','pdf','Bandwidth',bw);%,'Bandwidth',bw
            [f_sim,xi_sim] = ksdensity(sim_data,pts,...
                'Function','pdf','Bandwidth',bw);%
            
            if 0
                plot(xi_exp,f_exp,'k','LineWidth',2); hold on
                plot(xi_sim,f_sim,'r','LineWidth',2)
                legend('kernel-exp','kernel-sim',...
                    'Location','northwest')
            end
            
            w_dis = w_distance(exp_data, sim_data, 2);
            W_diff(i_codon,i_sti) = w_dis;
            
        end

        good_fit_pop_wdis(i_ligand,i_dose) = mean(W_diff(:,i_sti));
        i_sti = i_sti+1;
        
    end
    
end

%change the order or not
change_order = 0;

if change_order
    i_ligand_order = 1:length(sti_vec);
    sti_vec = sti_vec(i_ligand_order);

    good_fit_pop_wdis = good_fit_pop_wdis(i_ligand_order,:);

end


if 1 % Figure S2D: tested on 05/10/2024    W-distance, run this for figure
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    index_reorder = [1,2,3,13,14,15,7,8,9,4,5,6,10,11,12];
    sti_vec = sti_vec(index_reorder);
    W_diff = W_diff(:,index_reorder);
    h = heatmap( sti_vec,codon_list,W_diff,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.2])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_wdis_202309'),'epsc')
    close()
end

if 1 % Figure S2E: tested on 05/10/2024
    figure(1)
    paperpos = [0,0,120,50]*1.8;
    papersize = [120,50]*1.8;
    draw_pos = [10,10,100,30]*1.8;
    
    index_ligand_reorder = [1,5,3,2,4];
    ligand_vec = ligand_vec(index_ligand_reorder);
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    c = categorical(ligand_vec);
    c = reordercats(c,ligand_vec);
    
    index_reorder = [1,2,3,13,14,15,7,8,9,4,5,6,10,11,12];
    
    bar_pct=good_fit_pop_wdis(index_ligand_reorder,:);
    b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
    ax2= gca;
    % ytickformat(ax2, '%g%%');
    % ylim([75,100])
    ylabel({'Wasserstein Distance' 'of exp. and sim.'})
    XL = xlim();
    ylim([0,0.5]);
    %plot(XL,[0.1,0.1],'--r','LineWidth',1.5)
    % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)
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
    saveas(gcf,strcat(fig_save_path,'bar_codon_good_fit_pop_wdis'),'epsc')
    close
end

% run this for figure
if 1 % Figure S2C: tested on 05/10/2024
    
        violin_plot_codon_2023_09(collect_feature_vects,fig_opt) % draw all codon and sti in one
        
end

