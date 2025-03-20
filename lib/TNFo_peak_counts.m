
%% get data
if 0
    vers = '20220504_2';
    addpath('../CommonUsedFunction')
    codon_path = '../NFkB_codon/';
    data_info.save_file_path = './data/';

    % data_info.save_file_name = strcat('tnfo_single_para_sample_',vers); % only beginning, no .mat
    % gene_info.parameter_name_vec = {{'params54'},{'params61'},{'params62'},{'params63'},{'params65'}};
    
    
    data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers); % only beginning, no .mat
    gene_info.parameter_name_vec = {{'params54','params61'},{'params54','params62'},{'params54','params63'},{'params54','params65'}};
    
    
    
    i_para_name = 1;
    
    sim_info.parameter_name = gene_info.parameter_name_vec{i_para_name};
    
    if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
        save_file_name = data_info.save_file_name;
        for i_para_name = 1:length(sim_info.parameter_name)
            save_file_name = strcat(save_file_name,'_',sim_info.parameter_name{i_para_name});
        end
        % save_file_name = strcat(save_file_name,'.mat');
        load(strcat(data_info.save_file_path,save_file_name,'.mat'),'sim_data_tbl');
    end
    
    
    index_nuc_NFkB_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'wt';
    
    index_nuc_NFkB_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'tnfo';
    
    index_TNFR_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='TNFR'...
        & sim_data_tbl.type == 'wt';
    
    index_TNFR_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='TNFR'...
        & sim_data_tbl.type == 'tnfo';
    
    index_IKK_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IKK'...
        & sim_data_tbl.type == 'wt';
    
    index_IKK_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IKK'...
        & sim_data_tbl.type == 'tnfo';
    IKK_wt = mean(sim_data_tbl.trajectory(index_IKK_wt,:));
    IKK_tnfo = mean(sim_data_tbl.trajectory(index_IKK_tnfo,:));
    
    TNFR_wt = mean(sim_data_tbl.trajectory(index_TNFR_wt,1));
    TNFR_tnfo = mean(sim_data_tbl.trajectory(index_TNFR_tnfo,1));
    
    traj_wt = sim_data_tbl.trajectory(index_nuc_NFkB_wt,:);
    [~,traj_order] = sort(max(traj_wt,[],2),'descend');
    traj_wt = traj_wt(traj_order,:);
    traj_wt = traj_wt(1:5:end,:);
    traj_tnfo = sim_data_tbl.trajectory(index_nuc_NFkB_tnfo,:);
    [~,traj_order] = sort(max(traj_tnfo,[],2),'descend');
    traj_tnfo = traj_tnfo(traj_order,:);
    traj_tnfo = traj_tnfo(1:5:end,:);
    
    
    parameters_wt = sim_data_tbl.parameter_value(index_nuc_NFkB_wt,:);
    parameters_tnfo = sim_data_tbl.parameter_value(index_nuc_NFkB_tnfo,:);
    
    genot_vec = {'wt','tnfo'};
    for i_genotype = 1:2
        genotype = genot_vec{i_genotype};
        current_fold = pwd;
        traj = ...
            sim_data_tbl.trajectory(sim_data_tbl.dose_val == 10 ...
            & sim_data_tbl.species =='nucNFkB'...
            & sim_data_tbl.type == genotype, :);
        data.info_ligand{i_genotype} = 'TNF';
        
        data.model_sim{i_genotype} = traj(:,1:5:end);
        
        data.info_dose_index{i_genotype} = 3;
        data.info_dose_str{i_genotype} = '10ng/mL';
        data.info_num_cells{i_genotype} = size(data.model_sim{i_genotype},1);
        data.order{i_genotype} = (1:data.info_num_cells{i_genotype})';
    end
    
    data.exp = data.model_sim;
    
    cd(codon_path)
    vis_data_field = {'model_sim'};%,'sample'};
    data_label = {'simulation'};%,'sample'};
    [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
    
    clear osc_ratio
    for i_genotype = 1:2
        
        osc_tag = osc_cats(metrics{1}.peakfreq,metrics{1}.off_times);
        osc_ratio(i_genotype) = sum( osc_tag == 'osc')/length(osc_tag);
    end
    
    cd(current_fold)
end

%% Ade peakfreq
if 1
    [~,traj_order] = sort(metrics{1}.peakfreq,'descend');
    traj_wt = data.model_sim{1}(traj_order,:);
    
    [~,traj_order] = sort(metrics{2}.peakfreq,'descend');
    traj_tnfo = data.model_sim{2}(traj_order,:);
    figure(1)
    subplot(1,2,1)
    h=heatmap(traj_wt,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25],...
        'Title',{'wt',sim_info.parameter_name{1},sim_info.parameter_name{2}});
    
    subplot(1,2,2)
    h=heatmap(traj_tnfo,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25],...
        'Title',{'tnfo',sim_info.parameter_name{1},sim_info.parameter_name{2}});
    figure(2)
    h1 = histogram(metrics{1, 1}.peakfreq  ); hold on
    h1.Normalization = 'countdensity';
    h2 = histogram(metrics{2}.peakfreq  ); hold on
    h2.Normalization = 'countdensity';
    
end

if 0
    traj_wt = sim_data_tbl.trajectory(index_nuc_NFkB_wt,:);
    [~,traj_order] = sort(max(traj_wt,[],2),'descend');
    traj_wt = traj_wt(traj_order,:);
    traj_wt = traj_wt(1:5:end,:);
    traj_tnfo = sim_data_tbl.trajectory(index_nuc_NFkB_tnfo,:);
    [~,traj_order] = sort(max(traj_tnfo,[],2),'descend');
    traj_tnfo = traj_tnfo(traj_order,:);
    traj_tnfo = traj_tnfo(1:5:end,:);
    data.model_sim{1} =traj_wt(:,1:5:end);
    data.model_sim{2}  = traj_tnfo(:,1:5:end);
    
    pk_num = NaN(size(data.model_sim{1},1),2);
    osc_ratio = NaN(2,3);
    for i_genotype = 1:2
        traj = data.model_sim{i_genotype};
        for i_traj = 1:size(traj,1)
            [pks,~] = findpeaks(traj(i_traj,:));
            % plot(1:size(traj,2),traj(i_traj,:)); hold on
            [trghs,~] = findpeaks(-traj(1,:));
            trghs = -trghs;
            
            clear p2t
            p2t(1:2:length(pks)*2) = pks;
            p2t(2:2:length(trghs)*2) = trghs;
            
            pk_num(i_traj,i_genotype) = sum(abs(p2t(1:end-1) - p2t(2:end))>0.025) /2 ;
        end
        
        for i_threshold =1:3
            osc_ratio(i_genotype,i_threshold) = sum(pk_num(:,i_genotype) > i_threshold )/size(pk_num,1);
        end
    end
    
    
    [~,traj_order] = sort(pk_num(:,1),'descend');
    traj_wt = data.model_sim{1}(traj_order,:);
    
    [~,traj_order] = sort(pk_num(:,2),'descend');
    traj_tnfo = data.model_sim{2}(traj_order,:);
    
    % figure(1)
    % subplot(1,2,1)
    % h1=heatmap(traj_wt,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25],...
    %     'Title',{'wt',sim_info.parameter_name{1},sim_info.parameter_name{2}});
    %
    % subplot(1,2,2)
    % h2=heatmap(traj_tnfo,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25],...
    %     'Title',{'tnfo',sim_info.parameter_name{1},sim_info.parameter_name{2}});
    
    figure(2)
    
    for i_osratio = 1:size(osc_ratio,2)
        
        subplot(1,3,i_osratio)
        
        b=bar(osc_ratio(:,i_osratio)*100);
        %xlabel()
        b.FaceColor = 'flat';
        b.CData(1,:) = [0 0 1];
        b.CData(2,:) = [1 0 0];
        
        ytickformat(gca, '%g%%');
        grid on
        ylabel('Osc ratio')
        set(gca,'XTickLabel',{'WT','Tnf^{-/-}'})
        set(gca,'fontsize',7,'fontweight','b')
    end
    %         Set_figure_size
    %         %         saveas(gcf,strcat('./',save_file_name,'_tnfr'),'epsc')
    %         close
    %
end

%%
if 1
    
