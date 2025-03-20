function [] = draw_all_ligand_codon_fit_sampling_distrib_2023(data_save_file_path,fig_save_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
fig_opt.paper_opt.papersize=[220 180]*3;

% data_proj_num_vec = data_proj_nums{i_ligand};
%     load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

%
% for i_sti = 1:length(data.pred_mode_filter_nan)
%
%     data.pred_mode_cv{i_sti} = std(data.pred_mode_filter_nan{i_sti},[],2)./mean(data.pred_mode_filter_nan{i_sti},2);
%     [~,data.pred_mode_cv_order{i_sti}] = sort(data.pred_mode_cv{i_sti},'descend');
%     [~,data.pred_mode_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2},'descend');
%     [~,data.exp_osc_order{i_sti}] = sort(collect_feature_vects.OscVsNonOsc{i_sti*2-1},'descend');
%
% end

%
% % input
% thresh_TNF=0.33;
% thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
% for i_data = 1:3
%     index = data.(thresh_field{1}){i_data} (1: ceil(length(data.(thresh_field{1}){i_data}) * thresh_TNF));
%     data.pred_mode{i_data} = data.pred_mode_filter_nan{i_data}(index ,:);
%     data.exp{i_data} = data.exp_mode_filter_nan{i_data}(index,:);
%     [~, data.order{i_data}]=sort(max(data.exp{i_data},[],2),'descend');
% end
%
% codon_list = {'Duration','EarlyVsLate','OscVsNonOsc','PeakAmplitude','Speed','TotalActivity'};
%
% for i_codon =1:6
%     for i_data = 2:3
%         index = data.(thresh_field{1}){i_data} (1: ceil(length(data.(thresh_field{1}){i_data}) * thresh_TNF));
%
%         collect_feature_vects.(codon_list{i_codon}){i_data*2-1} = collect_feature_vects.(codon_list{i_codon}){i_data*2-1}( index,:);
%         collect_feature_vects.(codon_list{i_codon}){i_data*2} = collect_feature_vects.(codon_list{i_codon}){i_data*2}( index,:);
%     end
% end


%% codon distribution
% violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
% codon and sti seperately, and save
fig_opt.save_file = strcat(fig_save_path,'All_ligand_codon_2023_03_distrib');

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
            %         pts_exp = linspace(min(exp_data),max(exp_data),50);
            %         pts_sim = linspace(min(sim_data),max(sim_data),50);
            
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
            
            l2_diff_kde(i_codon,i_sti) = norm(f_exp-f_sim)*(xi_exp(2)-xi_exp(1));
        end
        
        
        good_fit_pop(i_ligand,i_dose) = mean(l2_diff_kde(:,i_sti));
        i_sti = i_sti+1;
        
    end
    
end

if 0
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,l2_diff_kde,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    
    
    saveas(gcf,strcat(fig_save_path,'codon_diff'),'epsc')
    close()
end

if 1
    figure(1)
    paperpos = [0,0,120,50]*1.8;
    papersize = [120,50]*1.8;
    draw_pos = [10,10,100,30]*1.8;
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    c = categorical(ligand_vec);
    c = reordercats(c,ligand_vec);
    bar_pct=good_fit_pop;
    b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
    ax2= gca;
    % ytickformat(ax2, '%g%%');
    % ylim([75,100])
    ylabel({'discrepancy between' 'exp. and sim.'})
    XL = xlim();
    ylim([0,0.5]);
    plot(XL,[0.1,0.1],'--r','LineWidth',1.5)
    % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)
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
    saveas(gcf,strcat(fig_save_path,'bar_codon_good_fit_pop'),'epsc')
    close
end
%    violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
if 0
        violin_plot_codon_2023_04(collect_feature_vects,fig_opt) % draw all codon and sti in one

end


if 0
    violin_plot_codon_2023_03_figure1(collect_feature_vects)
end
% violin_plot_codon_2023_CRI(collect_feature_vects,fig_opt)




