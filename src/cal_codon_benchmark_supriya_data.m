data_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/Experiments/Benchmark_Supriya/';
fig_save_path = '../SubFigures2023/';

cal_save_data = 0;
if cal_save_data
    data_file_list_r1 = {'rescaled_Bench_Supriya_CpG1.xlsx'
        'rescaled_Bench_Supriya_LPS1.xlsx'
        'rescaled_Bench_Supriya_P3K1.xlsx'
        'rescaled_Bench_Supriya_PIC1.xlsx'
        'rescaled_Bench_Supriya_TNF1.xlsx'};
    
    data_file_list_r2 = {'rescaled_Bench_Supriya_CpG2.xlsx'
        'rescaled_Bench_Supriya_LPS2.xlsx'
        'rescaled_Bench_Supriya_P3K2.xlsx'
        'rescaled_Bench_Supriya_PIC2.xlsx'
        'rescaled_Bench_Supriya_TNF2.xlsx'};
    
    ligand_str = {'CpG','LPS','Pam3CSK','PolyIC','TNF'};
    dose_str = {'x1','x1','x1','x1','x1'};
    
    for i_data_file = 1:length(data_file_list_r1)
        % for debugging
        %histogram(max(data,[],2));hold on
        data_exp = readmatrix(strcat(data_path,data_file_list_r1{i_data_file}));
        [row,~] = find(isnan(data_exp));
        index_non_NaN = setdiff(1:size(data_exp,1),row);
        data.exp_r1{i_data_file} = data_exp(index_non_NaN,:);
        data.exp{i_data_file} = data.exp_r1{i_data_file};
        data_exp  = readmatrix(strcat(data_path,data_file_list_r2{i_data_file}));
        [row,~] = find(isnan(data_exp));
        index_non_NaN = setdiff(1:size(data_exp,1),row);
        data.exp_r2{i_data_file} =data_exp(index_non_NaN,:);
        
        data.info_ligand{i_data_file} = ligand_str{i_data_file};
        data.info_dose_str{i_data_file} = dose_str{i_data_file};
    end
    
    vis_data_field = {'exp_r1','exp_r2'}; %,'sample'};
    data_label = {'experiment_r1','experiment_r2'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon_2023(data,vis_data_field,data_label); %,  parameter
    
    save(strcat(data_path,'benchmark_supriya_data_codon_metric.mat'),'data','metrics','collect_feature_vects');
else
    load(strcat(data_path,'benchmark_supriya_data_codon_metric.mat'));
end

codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };


i_sti = 1;
for i_data = 1:length(data.exp_r1)%TNF,LPS [27];%
    
    
    for i_codon = 1:length(codon_list)
        
        %  sti_vec{i_sti} = strcat(data.info_ligand{i_data},'-',data.info_dose_str{i_data});
        sti_vec{i_sti} = data.info_ligand{i_data};
        exp_data = collect_feature_vects.(codon_list{i_codon}){i_data*2-1};
        exp_data = exp_data*1;
        sim_data = collect_feature_vects.(codon_list{i_codon}){i_data*2};
        sim_data = sim_data*1;
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
            l1_diff_kde(i_codon,i_sti) = norm(f_exp-f_sim,1)*(xi_exp(2)-xi_exp(1));
            jsd = dist_by_JSD({exp_data,sim_data});
            JSD_diff(i_codon,i_sti) =  jsd(3);
            w_dis = w_distance(exp_data, sim_data, 2)
            W_diff(i_codon,i_sti) = w_dis;
    end
    
    
    good_fit_pop(i_data,1) = mean(l2_diff_kde(:,i_sti));
            good_fit_pop_l1(i_data,1) = mean(l1_diff_kde(:,i_sti));
good_fit_pop_JSD(i_data,1) = mean(JSD_diff(:,i_sti));
good_fit_pop_Wdis(i_data,1) = mean(W_diff(:,i_sti));
    i_sti = i_sti+1;
    
end

i_ligand_order = [5,2,1,4,3];
sti_vec = sti_vec(i_ligand_order);
l2_diff_kde = l2_diff_kde(:,i_ligand_order);
l1_diff_kde = l1_diff_kde(:,i_ligand_order); 
JSD_diff = JSD_diff(:,i_ligand_order);
W_diff = W_diff(:,i_ligand_order);

good_fit_pop =  good_fit_pop(i_ligand_order,:);
good_fit_pop_l1 = good_fit_pop_l1(i_ligand_order,:);
good_fit_pop_JSD = good_fit_pop_JSD(i_ligand_order,:);
good_fit_pop_Wdis =  good_fit_pop_Wdis(i_ligand_order,:);


if 0 %l2 norm
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,250,200];
    paper_size = [250,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,l2_diff_kde,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.2])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_supriya_exp_benchmark'),'epsc')
    close()
end


if 0 % jsd
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,JSD_diff,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_jsd_supriya_exp_benchmark'),'epsc')
    close()
end


if 1 % W-distance
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,250,200];
    paper_size = [250,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,W_diff,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.2])
    max(max(W_diff))
    mean(W_diff)
    %saveas(gcf,strcat(fig_save_path,'codon_diff_supriya_exp_benchmark_Wdis'),'epsc')
    close()
end