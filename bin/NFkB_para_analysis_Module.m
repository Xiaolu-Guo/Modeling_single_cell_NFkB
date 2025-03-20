function [] = NFkB_para_analysis_Module(module_tocal,vers,fold_change_vec,para_setting_sheet,data_save_file_path)

% module_all = {'polyIC'};
sim_info.para_setting_sheet = para_setting_sheet;
% senstive_para.TNF = cell(1,1);
% sens_order.TNF = cell(1,1);
% sim_info.Module = 'TNF'; %'core','coreup','PolyIC'
%only
sim_info.Module = module_tocal;
% sim_info.
% sim_info.parameter = {'params54','params55','params56'};

if isfield(sim_info,'parameter')
    for i_parameter_num = 1:length(sim_info.parameter)
        sim_info.para_fold{i_parameter_num} = fold_change_vec;
    end
else
    sim_info.para_fold{1} = fold_change_vec;
    sim_info.default_ratio = [0.1,1,10];
end


sim_info.stimuli = {'TNF','CpG','Pam3CSK','LPS','polyIC'};
%% data_info initializing
data_info.species = {'IkBaNFkBn','NFkBn'};% {'TNFR','IKK','IkBaNFkBn','NFkBn'};{'TNF','TNFR','TNFR_TNF','TTR','C1_off','C1','IkBa','IkBan','IKKIkBa','IKKIkBaNFkB','IkBaNFkB','IkBaNFkBn','IKK','NFkBn'};
data_info.save_file_path = data_save_file_path;

if isfield(sim_info,'parameter')
    data_info.save_file_name =strcat('TNFo_',sim_info.stimuli);
    for i_para =1:length(sim_info.parameter )
        data_info.save_file_name = strcat(data_info.save_file_name,sim_info.parameter{i_para},'_');
    end
    data_info.save_file_name = strcat(data_info.save_file_name,vers,'.mat');
else
    data_info.save_file_name = strcat('TNFo_',sim_info.Module,vers,'.mat');
end

%% run sim
NFkB_para_sens_analysis_sim(sim_info,data_info)%(module,vers,names);
%
%  [sens_order.(module), senstive_para.(module)] = NFkB_sens_para_pick(module,vers);
% sens_para_common.(module) = NFkB_sens_para_common(senstive_para.(module));
%
%% figure_info initializing
figure_info.save_figure_path = strcat('../../NFkB_para_estm_project/NFkB_figures/',vers,'_Sensitive_parameter_analysis/');

%         NFkB_sens_para_draw_results(data_info,figure_info)
codon_path = '../NFkB_codon_local/';
figure_info.codon_fields = {'OscVsNonOsc','Duration','EarlyVsLate','PeakAmplitude','Speed','TotalActivity'};
figure_info.metric_fields = {'oscpower'};
figure_info.save_codon_figure_path = strcat('../../NFkB_para_estm_project/NFkB_figures/',vers,'_Sensitive_parameter_codon_analysis/');

%% run codon cal
% [collect_feature_vects,metrics,sim_table_data] = TNFo_p54_codon_draw_results_alldose(data_info,figure_info,codon_path);


%%
% rand_mat = randn(20,1);

%
% fold_change_log = 0 + rand_mat* (log(2)/3);
%
% fold_change_vec = exp(fold_change_log);
%
% fold_change_tnfo_log = log(1.26) + rand_mat* (log(2)/3);
%
% fold_change_tnfo_vec = exp(fold_change_tnfo_log);
%
% fold_change_tnfo_smVar_log = 0 + rand_mat* (log(1.3)/3);
%
% fold_change_tnfo_smVar_vec = exp(fold_change_tnfo_smVar_log);
%
% histogram(fold_change_vec);hold on
% histogram(fold_change_tnfo_vec);hold on
% histogram(fold_change_tnfo_smVar_vec);hold on
%
% set(gca,'XScale','log')
