function [] = draw_Sens_codon(data_save_file_path,fig_save_path)
sti_all = {'core','TNF','CpG','Pam3CSK','LPS','PolyIC'};%,'core','coreup 'TNF','CpG','Pam3CSK','LPS',
% module_all = {'polyIC'};

% senstive_para.TNF = cell(1,1);
% sens_order.TNF = cell(1,1);
% module_all =  {'core','TNF','CpG','LPS','Pam3CSK','polyIC'};%{'polyIC'};%

% sim_info.Module = 'TNF'; %'core','coreup','PolyIC'
% only

for i_sti_all = 1 :length(sti_all)
    sim_info.stimuli = sti_all{i_sti_all};
    switch sim_info.stimuli
        case 'TNF'
            sim_info.parameter = {'params54';'params61'};
            
        case 'CpG'
            sim_info.parameter = {'params85';'params93'};
            
        case 'LPS'
            sim_info.parameter = {'params35';'params40';'params44'};
            
        case 'polyIC'
            sim_info.parameter = {'params77';'params79';'params83'};
            
        case 'Pam3CSK'
            sim_info.parameter = {'params68';'params75'};
            
        case 'core'
            sim_info.parameter = {'params66';'params99';'params100';'params101'};
            
    end
    data_info.species = {'nucNFkB'};% {'TNF','TNFR','TNFR_TNF','TTR','C1_off','C1','IkBa','IkBan','IKKIkBa','IKKIkBaNFkB','IkBaNFkB','IkBaNFkBn','IKK','NFkBn'};
    % data_info.save_file_path = './data/';
    data_info.save_file_path = data_save_file_path;
    
    %     if isfield(sim_info,'parameter')
    %         data_info.save_file_name =strcat('data_sens_',sim_info.stimuli);
    %         for i_para =1:length(sim_info.parameter )
    %             data_info.save_file_name = strcat(data_info.save_file_name,sim_info.parameter{i_para},'_');
    %         end
    %         data_info.save_file_name = strcat(data_info.save_file_name,'.mat');
    %     else
    sim_info.Module = sim_info.stimuli;
    data_info.save_file_name = strcat('Module_Sens_',sim_info.Module,'.mat');
    %     end
    
    %         NFkB_sens_para_analysis_sim(sim_info,data_info)%(module,vers,names);
    %
    %  [sens_order.(module), senstive_para.(module)] = NFkB_sens_para_pick(module,vers);
    % sens_para_common.(module) = NFkB_sens_para_common(senstive_para.(module));
    %
    figure_info.save_figure_path = fig_save_path;
    
    %         NFkB_sens_para_draw_results(data_info,figure_info)
    % codon_path = '../NFkB_codon/';
    % figure_info.codon_fields = {'Duration'};
    figure_info.metric_fields = {'duration','oscpower','max_pos_integral','pos_pk1_time','max_value','time2HalfMaxPosIntegral','time2HalfMaxValue'};%'max_pos_pk1_speed',,'pos_pk1_amp','time2HalfMaxPosIntegral'
    figure_info.save_codon_figure_path =fig_save_path;
    
    data_info.y_metric_lim = {[0,8],[0,1e-3],[0,2],[0,4],[0,0.3],[0,6],[0,2]};
    
    NFkB_draw_pick_para_codon(data_info,figure_info)
end
end

