function [good_fit_pop,diff_distr,sti_vec] = cal_codon_distri_diff(collect_feature_vects,index_data1,index_data2,diff_method)
% fig_opt.save_file = strcat(fig_save_path,'All_ligand_fit_sampling_codon_2023_05_distrib');

if length(index_data1) ~= length(index_data2)
    error('the lengths of comparing data sets have to be same')
end

if nargin < 4
    diff_method = 'L2';
end

codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
sti_vec = cell(1,length(index_data1));
good_fit_pop_l2 = NaN(length(index_data1),1); % l2 norm difference between two codon distribution
good_fit_pop_l1 = NaN(length(index_data1),1);% l1 norm difference between two codon distribution
good_fit_pop_JSD = NaN(length(index_data1),1);% jsd difference between two codon distribution

i_sti = 1;
for i_data = 1:length(index_data1)%TNF,LPS [27];%
    
    % i_data = find(categorical(data.info_ligand)==ligand_vec{i_ligand} & categorical(data_dose_str)==dose_vec{i_ligand}{i_dose});
    
    for i_codon = 1:length(codon_list)
        
        %  sti_vec{i_sti} = strcat(data.info_ligand{i_data},'-',data.info_dose_str{i_data});
        sti_vec{i_sti} = collect_feature_vects.info_ligand{index_data1(i_data)};
        data1 = collect_feature_vects.(codon_list{i_codon}){index_data1(i_data)};
        data2 = collect_feature_vects.(codon_list{i_codon}){index_data2(i_data)};
        pts = linspace(min(min(data1),min(data2)),max(max(data1),max(data2)),50);
        %         pts_exp = linspace(min(exp_data),max(exp_data),50);
        %         pts_sim = linspace(min(sim_data),max(sim_data),50);
        
        [~,~,bw_exp] = ksdensity(data1,pts,...
            'Function','pdf');
        [~,~,bw_sim] = ksdensity(data2,pts,...
            'Function','pdf');%
        bw = min(bw_exp,bw_sim);
        
        [f_exp,xi_exp] = ksdensity(data1,pts,...
            'Function','pdf','Bandwidth',bw);%,'Bandwidth',bw
        [f_sim,xi_sim] = ksdensity(data2,pts,...
            'Function','pdf','Bandwidth',bw);%
        
        if 0
            plot(xi_exp,f_exp,'k','LineWidth',2); hold on
            plot(xi_sim,f_sim,'r','LineWidth',2)
            legend('kernel-exp','kernel-sim',...
                'Location','northwest')
        end
        
        l2_diff_kde(i_codon,i_sti) = norm(f_exp-f_sim)*(xi_exp(2)-xi_exp(1));
        l1_diff_kde(i_codon,i_sti) = norm(f_exp-f_sim,1)*(xi_exp(2)-xi_exp(1));
        jsd = dist_by_JSD({data1,data2});
        JSD_diff(i_codon,i_sti) =  jsd(3);
    end
    
    
    good_fit_pop_l2(i_data) = mean(l2_diff_kde(:,i_sti));
    good_fit_pop_l1(i_data) = mean(l1_diff_kde(:,i_sti));
    good_fit_pop_JSD(i_data) = mean(JSD_diff(:,i_sti));
    i_sti = i_sti+1;
    
    
end

switch diff_method
    case 'L2'
        good_fit_pop = good_fit_pop_l2;
        diff_distr = l2_diff_kde;
    case 'L1'
        good_fit_pop = good_fit_pop_l1;
        diff_distr = l1_diff_kde;
    case 'JSD'
        good_fit_pop = good_fit_pop_JSD;
        diff_distr = JSD_diff;
    otherwise
        error(strcat('do not support the method',diff_method))
end

end

