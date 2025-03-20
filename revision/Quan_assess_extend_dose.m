%% initialization 
addpath('../raw_data2023/2024/')
savepath = '../SubFigures2023/revision_2024/';
data_savepath = '../raw_data2023/2024/';
vers_savefig = '_1126';
data_file_list = {'1000nMCpG_test_trajectories';
    '1000ngP3C4_test_trajectories';
    '100ngLPS_test_trajectories';
    '100ugPIC_test_trajectories';
    '10ngTNF_test_trajectories';%5
    'CpG_100nM_trajectories';
    'CpG_1uM_trajectories';
    'CpG_333nM_trajectories';
    'CpG_33nM_trajectories';
    'CpG_TNF_trajectories';%10
    'Ctrl_trajectories';
    'IkBaKO_TNF_3ng_trajectories';
    'IkBaMut_LPS_100ng_r1_trajectories';
    'IkBaMut_PolyIC_50ug_r1_trajectories';
    'IkBaMut_TNF_10ng_r1_trajectories';%15
    'IkBaWt_LPS_100ng_r1_trajectories';
    'IkBaWt_PolyIC_50ug_r1_trajectories';
    'IkBaWt_TNF_10ng_r1_trajectories';
    'IkBbeKO_LPS_10ng_trajectories';
    'IkBbeKO_TNF_3ng_trajectories';%20
    'LPS_0_3ng_trajectories';
    'LPS_100ng_trajectories';
    'LPS_10ng_trajectories';
    'LPS_1ng_trajectories';
    'LPS_333ng_trajectories';%25
    'LPS_33ng_r1_trajectories';
    'LPS_33ng_r2_trajectories';
    'LPS_33ng_r3_trajectories';
    'LPS_33ng_trajectories';
    'LPS_3ng_trajectories';%30
    'Pam3CSK4_100ng_trajectories';
    'Pam3CSK4_10ng_trajectories';
    'Pam3CSK4_1ng_trajectories';
    'Pam3CSK4_33ng_trajectories';
    'Pam3CSK4_3ng_trajectories';%35
    'PolyIC_100ug_trajectories';
    'PolyIC_10ug_trajectories';
    'PolyIC_33ug_trajectories';
    'PolyIC_3ug_trajectories';
    'TNF_0_3ng_trajectories';%40
    'TNF_10ng_trajectories';
    'TNF_1ng_trajectories';
    'TNF_33ng_trajectories';
    'TNF_3ng_trajectories'};

%TNF 33ng; Pam 33ng; CpG 1000nM; LPS 333 ng; pIC 3ug;
index_compare = [43,34,7,25,39];
extract_info_from_file_list

i_data_file = 22;
data_traj = readmatrix(strcat(data_file_list{i_data_file},'.csv'));
rescale_factor = rescale_factor_cal(data_traj);

%% visualize heatmap extended dose simulation for comparing with ade2021 experiments

if 1 % before 11/29: to delete before submission
    %TNF 33ng; LPS 333 ng; pIC 50ug;
    index_compare = [43,25,17];
    vers_data = '_2023_05';%'_202402';
    data_filename = strcat('Sim_ext_dose_codon_metric',vers_data,'_revision_1130.mat');

    data_save_file_path_1 = '../raw_data2023/2024/';%_fay_parameter/';
    load(strcat(data_save_file_path_1,data_filename))

    i_codon_cal_file = 1;
    for i_ligand_index = 1:length(data.exp)

        if 1 % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
            data_plot.model_sim{i_ligand_index} = data.model_sim{i_ligand_index}(1:9:end,:);
        end

        color_lmt = [0,0.25];
        data_first_half = data_plot.model_sim{i_ligand_index}(:,1:size(data.model_sim{i_ligand_index},2));
        data_first_half(isnan(data_first_half)) = 0;
        [~,data_order_index] = sort(sum(data_first_half,2),"descend");

        fig_gcf = plot_heatmap(data_plot.model_sim{i_ligand_index},  data_order_index,  color_lmt);
        saveas(gcf,strcat(savepath,data_file_list{index_compare(i_ligand_index)},'_model_sim',vers_data),'pdf');
        close

        data_traj = readmatrix(strcat(data_file_list{index_compare(i_ligand_index)},'.csv'));
        data_traj = data_traj * rescale_factor;
        % fig_gcf = plot_heatmap(data_traj);
        % saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
        % close

        data_traj = data_traj(:,1:97);%97 for 8hours, 139 for 11.5hours
        [row,~] = find(isnan(data_traj));
        index_non_NaN = setdiff(1:size(data_traj,1),row);
        if isempty(index_non_NaN)
        else
            data_plot.exp{i_codon_cal_file} = data_traj(index_non_NaN,:);
        end

        %[~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
        data_first_half = data_plot.exp{i_codon_cal_file}(:,1:size(data_plot.exp{i_codon_cal_file},2));
        data_first_half(isnan(data_first_half)) = 0;
        [~,data_order_index] = sort(sum(data_first_half,2),"descend");

        fig_gcf = plot_heatmap(data_plot.exp{i_codon_cal_file},  data_order_index,  color_lmt);
        saveas(gcf,strcat(savepath,data_file_list{index_compare(i_ligand_index)},'_exp'),'pdf');
        close

        i_codon_cal_file = i_codon_cal_file + 1;
    end
end

if 0 % before 11/29: to delete before submission
    %TNF 33ng; Pam 33ng; CpG 1000nM; LPS 333 ng; pIC 3ug;
    dose_index_vec = [17,9,15,21,9];
    index_all_doses_vec =[0,21*4,21*2,21*1,21*3];
    file_vec = [1,5,3,2,4];
    i_codon_cal_file = 1;
    for i_file = 1:5
        i_dose_index = dose_index_vec(i_file);
        data_filename = strcat('Sim3_codon_r',num2str(file_vec(i_file)),'_metric.mat');

        vis_data_field = {'pred_mode_total_1117'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};

        data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

        load(strcat(data_save_file_path_1,data_filename))

        if 1 % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
            data_plot.model_sim{i_file} = data.model_sim{i_dose_index}(1:9:end,:);
        end

        color_lmt = [0,0.25];
        data_first_half = data_plot.model_sim{i_file}(:,1:size(data.model_sim{i_file},2)/2);
        data_first_half(isnan(data_first_half)) = 0;
        [~,data_order_index] = sort(sum(data_first_half,2),"descend");

        fig_gcf = plot_heatmap(data_plot.model_sim{i_file},  data_order_index,  color_lmt);
        % saveas(gcf,strcat(savepath,data_file_list{i_data_file},'_model_sim'),'epsc');
        close
    end
end

%% visualize heatmap 5 doses for model extended doses simulation

if 0 % for debugging Pam

    dose_index_vec = [2,31:35];
    i_codon_cal_file =  1;
    color_lmt = [0,0.25];
    for i_data = 1:length(dose_index_vec) %length(data_file_list)
        i_data_file = dose_index_vec(i_data);
        data_traj = readmatrix(strcat(data_file_list{i_data_file},'.csv'));
        data_traj = data_traj * rescale_factor;
        % fig_gcf = plot_heatmap(data_traj);
        % saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
        % close

        data_traj = data_traj(:,1:97);%97 for 8hours, 139 for 11.5hours
        [row,~] = find(isnan(data_traj));
        index_non_NaN = setdiff(1:size(data_traj,1),row);
        if isempty(index_non_NaN)
        else
            data.exp{i_codon_cal_file} = data_traj(index_non_NaN,:);
        end

        %[~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
        data_first_half = data.exp{i_codon_cal_file}(:,1:size(data.exp{i_codon_cal_file},2)/2);
        data_first_half(isnan(data_first_half)) = 0;
        [~,data_order_index] = sort(sum(data_first_half,2),"descend");

        fig_gcf = plot_heatmap(data.exp{i_codon_cal_file},  data_order_index,  color_lmt);
        saveas(gcf,strcat(savepath,'debug_pam_',data_file_list{i_data_file}),'pdf');
        close

        i_codon_cal_file = i_codon_cal_file +1;

    end
end

if 0
    cal_codon = 0;
    if cal_codon
        i_codon_cal_file = 1;
        for i_file = 1:length(index_compare)
            i_data_file = index_compare(i_file);
            data_traj = readmatrix(strcat(data_file_list{i_data_file},'.csv'));
            data_traj = data_traj * rescale_factor;
            % fig_gcf = plot_heatmap(data_traj);
            % saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
            % close

            data_traj = data_traj(:,1:97);%97 for 8hours, 139 for 11.5hours
            [row,~] = find(isnan(data_traj));
            index_non_NaN = setdiff(1:size(data_traj,1),row);
            if isempty(index_non_NaN)
            else
                data.exp{i_codon_cal_file} = data_traj(index_non_NaN,:);
                data.info_ligand{i_codon_cal_file} = ligand_vec{i_data_file};
                data.info_dose_str{i_codon_cal_file} = dose_vec{i_data_file};
                data.info_genotype{i_codon_cal_file} = genotype_vec{i_data_file};
                data.info_test_data{i_codon_cal_file} = test_or_not_vec{i_data_file};
                data.info_replicate{i_codon_cal_file} = replicate_vec{i_data_file};
                i_codon_cal_file = i_codon_cal_file +1;
            end
        end

        vis_data_field = {'exp'}; %,'sample'};
        data_label = {'experiments'}; %,'sample'};
        [collect_feature_vects,metrics] = calculate_codon_2023(data,vis_data_field,data_label); %,  parameter
        a = 1;
        save(strcat(data_savepath,'Ade2021_ext_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');
    else
        %visualize ade's dataset for mutant
        load(strcat(data_savepath,'Ade2021_ext_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');
    end

    i_codon_cal_file =  1;
    color_lmt = [0,0.25];
    for i_file = 1:length(index_compare) %length(data_file_list)
        i_data_file = index_compare(i_file);

        %[~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
        color_lmt = [0,0.25];
        data_first_half = data.exp{i_codon_cal_file}(:,1:size(data.exp{i_codon_cal_file},2)/2);
        data_first_half(isnan(data_first_half)) = 0;
        [~,data_order_index] = sort(sum(data_first_half,2),"descend");

        fig_gcf = plot_heatmap(data.exp{i_codon_cal_file},  data_order_index,  color_lmt);
        % saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
        close

        i_codon_cal_file = i_codon_cal_file +1;

    end

end


%% quantitative assessment of extend dose between sim and exp
if 0

    cal_codon = 0;
    if cal_codon

        %TNF 33ng; Pam 33ng; CpG 1000nM; LPS 333 ng; pIC 3ug;
        ligand_vec_model = {'TNF', 'Pam', 'CpG', 'LPS', 'pIC'};
        dose_vec_model = {'33ng', '33ng', '1000nM', '333ng', '3ug'};

        dose_index_vec = [17,9,15,21,9];
        index_all_doses_vec =[0,21*4,21*2,21*1,21*3];
        file_vec = [1,5,3,2,4];
        i_codon_cal_file = 1;
        for i_file = 1:5
            i_dose_index = dose_index_vec(i_file);
            data_filename = strcat('Sim3_codon_r',num2str(file_vec(i_file)),'_metric.mat');

            vis_data_field = {'pred_mode_total_1117'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};

            data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

            load(strcat(data_save_file_path_1,data_filename))

            if 1 % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
                data_all.exp{i_file} = data.model_sim{i_dose_index}(1:9:end,:);
                data_all.info_ligand{i_codon_cal_file} = ligand_vec_model{i_file};
                data_all.info_dose_str{i_codon_cal_file} = dose_vec_model{i_file};
                data_all.info_genotype{i_codon_cal_file} = 'sim';
                data_all.info_test_data{i_codon_cal_file} = 'na';
                data_all.info_replicate{i_codon_cal_file} = '1';

                i_codon_cal_file = i_codon_cal_file +1;
            end
        end

        for i_file = 1:length(index_compare)
            i_data_file = index_compare(i_file);
            data_traj = readmatrix(strcat(data_file_list{i_data_file},'.csv'));
            data_traj = data_traj * rescale_factor;
            % fig_gcf = plot_heatmap(data_traj);
            % saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
            % close

            data_traj = data_traj(:,1:97);%97 for 8hours, 139 for 11.5hours
            [row,~] = find(isnan(data_traj));
            index_non_NaN = setdiff(1:size(data_traj,1),row);
            if isempty(index_non_NaN)
            else
                data_all.exp{i_codon_cal_file} = data_traj(index_non_NaN,:);
                data_all.info_ligand{i_codon_cal_file} = ligand_vec{i_data_file};
                data_all.info_dose_str{i_codon_cal_file} = dose_vec{i_data_file};
                data_all.info_genotype{i_codon_cal_file} = genotype_vec{i_data_file};
                data_all.info_test_data{i_codon_cal_file} = test_or_not_vec{i_data_file};
                data_all.info_replicate{i_codon_cal_file} = replicate_vec{i_data_file};
                i_codon_cal_file = i_codon_cal_file +1;
            end
        end

        vis_data_field = {'exp'}; %,'sample'};
        data_label = {'exp_and_sim'}; %,'sample'};
        [collect_feature_vects,metrics] = calculate_codon_2023(data_all,vis_data_field,data_all.info_genotype); %,  parameter
        collect_feature_vects.info_data_type = data_all.info_genotype';
        save(strcat(data_savepath,'Extend_dose_Ade2021_ext_Sim_codon_metric.mat'),'data_all','metrics','collect_feature_vects');
    else
        %visualize ade's dataset for mutant
        load(strcat(data_savepath,'Extend_dose_Ade2021_ext_Sim_codon_metric.mat'),'data','metrics','collect_feature_vects');
    end
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};


    if 0 % plot codon distributions for 5 ligands
        for i_codon = 1:length(codon_list)
            figure(1)
            paperpos=[0,0,350,100]*1.5;
            papersize=[350 100]*1.5;
            draw_pos=[10,10,330,90]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            y = collect_feature_vects.(codon_list{i_codon})(1:5); % y{i} size should be cellnum by 1
            z = collect_feature_vects.(codon_list{i_codon})(6:10); % 1000 by 10
            std_cal_y = [];
            for i_y = 1:length(y)
                std_cal_y = [std_cal_y;y{i_y}]
            end
            std_y = std(std_cal_y(:))/3500;
            % subplot(1,length(vis_data_field),i_data_field)

            al_goodplot_pair_RMSD_diff_size(y,[1:2:10],0.5,ones(size(y,2),1)*[255 0 0]/255 ,'bilateral',[],std_y); %left
            al_goodplot_pair_RMSD_diff_size(z,[2:2:10],0.5,ones(size(z,2),1)*[0 0 0]/255,'bilateral',[],std_y);

            xlim([0.4 10.6])

            xticks([1:10])
            xticklabels({})
            %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

            ylim([-0.1,1.1]);
            for i_x = 1:10
                plot([i_x,i_x],[0,5],'--','Color','k');hold on
            end
            set(gca,'fontsize',14,'fontname','Arial');
            %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

            saveas(gcf,strcat(savepath,'Codon_',codon_list{i_codon},'_ext_dose_signaling_codon',vers_savefig),'pdf');

            close
        end
    end

    if 0 % plot codon distributions for 3 ligands
        for i_codon = 1:length(codon_list)
            figure(1)
            paperpos=[0,0,350,100]*1.5;
            papersize=[350 100]*1.5;
            draw_pos=[10,10,330,90]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            y = collect_feature_vects.(codon_list{i_codon})([1,3,4]); % y{i} size should be cellnum by 1
            z = collect_feature_vects.(codon_list{i_codon})([6,8,9]); % 1000 by 10
            std_cal_y = [];
            for i_y = 1:length(y)
                std_cal_y = [std_cal_y;y{i_y}]
            end
            std_y = std(std_cal_y(:))/3500;
            % subplot(1,length(vis_data_field),i_data_field)

            al_goodplot_pair_RMSD_diff_size(y,[1:3:7],0.5,ones(size(y,2),1)*[255 0 0]/255 ,'bilateral',[],std_y); %left
            al_goodplot_pair_RMSD_diff_size(z,[2:3:8],0.5,ones(size(z,2),1)*[0 0 0]/255,'bilateral',[],std_y);

            xlim([0 9])

            xticks([0:9])
            xticklabels({})
            %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

            ylim([-0.1,1.1]);
            % for i_x = 1:10
            %     plot([i_x,i_x],[0,5],'--','Color','k');hold on
            % end
            set(gca,'fontsize',14,'fontname','Arial');
            %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

            saveas(gcf,strcat(savepath,'Codon_',codon_list{i_codon},'_ext_dose_signaling_codon_3ligand',vers_savefig),'pdf');

            close
        end
    end

    if 1 % plot heatmap differences across different conditions

        w_dis = [];

        for i_codon = 1:length(codon_list)
            for i_sti_exp = 1:5
                for i_sti_sim = 1:5

                    if i_sti_exp == i_sti_sim
                        W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                            w_distance(collect_feature_vects.(codon_list{i_codon}){i_sti_exp+5},...
                            collect_feature_vects.(codon_list{i_codon}){i_sti_sim}, 2);
                    elseif i_sti_exp <i_sti_sim
                        W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                        w_distance(collect_feature_vects.(codon_list{i_codon}){i_sti_exp},...
                            collect_feature_vects.(codon_list{i_codon}){i_sti_sim}, 2);    
                        %w_distance(sim_data{i_codon,i_sti_exp}, sim_data{i_codon,i_sti_sim}, 2);
                    else
                        W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                        w_distance(collect_feature_vects.(codon_list{i_codon}){i_sti_exp+5},...
                            collect_feature_vects.(codon_list{i_codon}){i_sti_sim+5}, 2); 
                            %w_distance(exp_data{i_codon,i_sti_exp}, exp_data{i_codon,i_sti_sim}, 2);
                    end

                end
            end
        end

        W_diff_mat = mean(W_diff,3);

        index_reorder = [1,3,4];

        W_diff_mat = W_diff_mat(index_reorder,index_reorder);

        %% Figure xx: 
        mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
            ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

        if 1
            figure (1)
            set(gcf, 'PaperUnits','points')

            paper_pos = [0,0,450,450]/6;
            paper_size = [450,450]/6;
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]


            heatmap(W_diff_mat,'colormap',mymap,'ColorbarVisible','off')
            caxis([0,0.25])

            % Getting handle of the heatmap
            h = gca;

            % Removing tick labels from x and y axes
            h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
            h.YDisplayLabels = repmat({''}, size(h.YDisplayData));

            saveas(gcf,strcat(savepath,'ext_dose_Wdis_exp_sim_3ligand_202411'),'pdf')
            close()

        end
    end

end

%% subfunctions
function [rescale_factor] = rescale_factor_cal(data)
% the macrophage cell volume is 6 pl, measured by Stefanie Luecke
NFkB_max_range = [0.25, 1.75]/6;
rescale_factor = NFkB_max_range(2) / prctile(max(data,[],2),90) ;% 100ng/mL LPS
end


function fig_gcf = plot_heatmap(data_traj,data_order_index,color_limits)

if nargin <2
    data_first_half = data_traj(:,1:size(data_traj,2)/2);
    data_first_half(isnan(data_first_half)) = 0;
    [~,data_order_index] = sort(sum(data_first_half,2),"descend");
end

if nargin <3
    color_limits = [-0.001,0.3];
end

%
% p = inputParser;
% % Required: ID input
% %valid_data_traj = @(x) assert(isnumeric(input_data_traj));
% addRequired(p,'data_traj');
% % Optional parameters to be passed to metrics function
% % valid_color_limit = @(x) assert(isnumeric(x)&size(x) == [1,2]);
% addParameter(p,'color_limits',[-0.001,0.3]);%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
% addParameter(p,'data_order_index',data_order_index_default); %allows adjustment of convection shift (?)
%
% % Parse the inputs
% parse(p, data_traj, 'data_order_index', data_order_index, 'color_limits', color_limits);
% color_limits = p.Results.color_limits;
% data_order_index = p.Results.data_order_index;

figure(1)
paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

cell_num=size(data_traj,1);
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

% subplot(1,length(vis_data_field),i_data_field)
h=heatmap(data_traj(data_order_index,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limits);%[-0.001,0.2] for TNF

XLabels = 0:5:((size(data_traj,2)-1)*5);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/60);
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:cell_num;
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig_gcf = gcf;

end