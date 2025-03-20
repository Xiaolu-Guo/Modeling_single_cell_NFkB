plot_traj_heatmap = 0; % plot heatmaps or not
cal_codon = 0; %calculate codon or not

vers = 'downsampling_';
data_savepath = '../raw_data2023/2024/';
%vers = 'original_';
savepath = '../SubFigures2023/revision_2024/';
ligand = 'Pam';
file_path = '../raw_data2023/SAEM_proj_2024_downsampling/XGESN2024002/';
%file_path = '../raw_data2023/SAEM_proj_2024_downsampling/XGESN2023005/';
exp_data_filename = 'down_sampling_SAEM_data_Pam3CSK480min.txt';
%exp_data_filename = 'down_sampling_SAEM_data_polyIC480min.txt';
%exp_data_filename = 'SAEM_data_polyIC480min.txt';

data_tbl = readtable(strcat(file_path,'results/predictions','.txt'));
exp_tbl =  readtable(strcat(file_path,exp_data_filename));
dose_all = unique(exp_tbl.AMT);
dose_index_vec = [];
for i_dose_all = 1:length(dose_all)
    if dose_all{i_dose_all} == '.'
    else
        dose_index_vec = [dose_index_vec,i_dose_all];
    end
end

dose_all = dose_all(dose_index_vec);
data_traj_sim = {};
data_traj_exp = {};
i_data_all = 1;
for i_dose = 1:length(dose_all)
    cell_id_start(i_dose) = exp_tbl.ID(find(strcmp(exp_tbl.AMT, dose_all{i_dose}),1,"first"));
    cell_id_end(i_dose) = exp_tbl.ID(find(strcmp(exp_tbl.AMT, dose_all{i_dose}),1,"last"));

    data_traj_sim{i_dose} = [];
    data_traj_exp{i_dose} = [];
    time_pts_num = 97;%??
    for i_cell = cell_id_start(i_dose):cell_id_end(i_dose)
        cell_traj_sim = data_tbl.indivPred_mode(data_tbl.id ==i_cell);
        cell_traj_exp = data_tbl.Y(data_tbl.id ==i_cell);
        if length(cell_traj_sim)==time_pts_num
            data_traj_sim{i_dose} = [data_traj_sim{i_dose} ;cell_traj_sim'];
            data_traj_exp{i_dose} = [data_traj_exp{i_dose} ;cell_traj_exp'];
        end
    end

    if plot_traj_heatmap
        data_traj = data_traj_sim{i_dose};
        fig_gcf = plot_heatmap(data_traj);
        saveas(gcf,strcat(savepath,vers,ligand,'_dose',num2str(i_dose),'_','sim'),'epsc');
        close

        data_traj = data_traj_exp{i_dose};
        fig_gcf = plot_heatmap(data_traj);
        saveas(gcf,strcat(savepath,vers,ligand,'_dose',num2str(i_dose),'_','exp'),'epsc');
        close
    end

    if cal_codon
        %remove nan
        index1 = isnan(sum(data_traj_exp{i_dose},2));
        index2 = isnan(sum(data_traj_sim{i_dose},2));
        index = (~index1)&(~index2);

        data_all.exp{i_data_all} = data_traj_exp{i_dose}(index,:);
        data_all.info_ligand{i_data_all} = ligand;
        data_all.info_dose_str{i_data_all} = dose_all{i_dose};
        data_all.info_genotype{i_data_all} = 'downsampling_exp';
        i_data_all = i_data_all +1;

        data_all.exp{i_data_all} = data_traj_sim{i_dose}(index,:);
        data_all.info_ligand{i_data_all} = ligand;
        data_all.info_dose_str{i_data_all} = dose_all{i_dose};
        data_all.info_genotype{i_data_all} = 'downsampling_sim';
        i_data_all = i_data_all+1;

    end

end

if cal_codon
    vis_data_field = {'exp'}; %,'sample'};
    data_label = {'exp_and_sim'}; %,'sample'};
    [collect_feature_vects,metrics] = calculate_codon_2023(data_all,vis_data_field,data_label); %,  parameter
    collect_feature_vects.info_data_type = data_all.info_genotype';
    save(strcat(data_savepath,vers,ligand,'_codon_metric.mat'),'data_all','metrics','collect_feature_vects');

    metrics_all = metrics;

    data_info.info_ligand =data_all.info_ligand;
    data_info.info_dose_str =data_all.info_dose_str;
    data_info.data_label = data_all.info_genotype;
    i_index = length(metrics) +1;

    data_save_file_path = '../raw_data2023/';
    load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'));


    field_list = fieldnames(metrics{1});

    for i_rpc = 1:length(metrics)

        for i_field = 1:length(field_list)
            metrics_all{i_index}.(field_list{i_field}) = metrics{i_rpc}.(field_list{i_field});
        end

        data_info.info_ligand(i_index) =data.info_ligand(ceil(i_rpc/2));
        data_info.info_dose_str(i_index) =data.info_dose_str(ceil(i_rpc/2));

        if mod(i_rpc,2)
            data_info.data_label(i_index) = {'exp'};
        else
            data_info.data_label(i_index) = {'sim'};
        end

        i_index = i_index+1;

    end


    [collect_feature_vects_all_minmax_scaled,metrics_all_minmax_scaled] = calculate_codon_from_metric2023(data_info,metrics_all); %,  parameter
    save(strcat(savepath,vers,ligand,'_with_all_exp_sim_codons.mat'),'collect_feature_vects_all_minmax_scaled','metrics_all_minmax_scaled','data_info')
else
    load(strcat(savepath,vers,ligand,'_with_all_exp_sim_codons.mat'))
end
codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

if 1 % plot codon distributions for 3 doses
    for i_codon = 1:length(codon_list)
        figure(1)
        paperpos=[0,0,350,100];
        papersize=[350 100];
        draw_pos=[10,10,330,90];
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        % load down sampled data
        y = collect_feature_vects_all_minmax_scaled.(codon_list{i_codon})([1,3,5]); % y{i} size should be cellnum by 1
        z = collect_feature_vects_all_minmax_scaled.(codon_list{i_codon})([2,4,6]); % 1000 by 10

        % load original data
        % y = collect_feature_vects_all_minmax_scaled.(codon_list{i_codon})([39,41,43]); % y{i} size should be cellnum by 1
        % z = collect_feature_vects_all_minmax_scaled.(codon_list{i_codon})([40,42,44]); % 1000 by 10

        std_cal_y = [];
        for i_y = 1:length(y)
            std_cal_y = [std_cal_y;y{i_y}]
        end
        std_y = std(std_cal_y(:))/3500;
        % subplot(1,length(vis_data_field),i_data_field)

        al_goodplot_pair_RMSD_diff_size(y,[1:3:7],0.5,ones(size(y,2),1)*[0 0 0]/255 ,'bilateral',[],std_y); %left
        al_goodplot_pair_RMSD_diff_size(z,[2:3:8],0.5,ones(size(z,2),1)*[255 0 0]/255,'bilateral',[],std_y);

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

        % save down sampled data
        saveas(gcf,strcat(savepath,vers,ligand,'Codon_',codon_list{i_codon}),'pdf');%,'_original_alldata_'

        % save original data
        % saveas(gcf,strcat(savepath,vers,'_original_alldata_',ligand,'Codon_',codon_list{i_codon}),'pdf');%,'_original_alldata_'

        close
    end
end


a = 1;

if 0 % plot heatmap differences across different conditions

    w_dis = [];
    basal_index = 19;
    for i_codon = 1:length(codon_list)
        for i_sti_exp = 1:3
            for i_sti_sim = 1:3

                if i_sti_exp == i_sti_sim
                    W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                        w_distance(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_exp)*2-1},...
                        collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_sim)*2}, 2);
                elseif i_sti_exp <i_sti_sim
                    W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                        w_distance(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_exp)*2},...
                        collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_sim)*2}, 2);
                    %w_distance(sim_data{i_codon,i_sti_exp}, sim_data{i_codon,i_sti_sim}, 2);
                else
                    W_diff(i_sti_exp,i_sti_sim,i_codon) =...
                        w_distance(collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_exp)*2-1},...
                        collect_feature_vects_all_minmax_scaled.(codon_list{i_codon}){(basal_index + i_sti_sim)*2-1}, 2);
                    %w_distance(exp_data{i_codon,i_sti_exp}, exp_data{i_codon,i_sti_sim}, 2);
                end

            end
        end
    end

    W_diff_mat = mean(W_diff,3);

    index_reorder = [1,2,3];

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

        saveas(gcf,strcat(savepath,'downsmapling_original_pam_wdis_202411'),'pdf')
        close()

    end
end
%% subfuncitons
function fig_gcf = plot_heatmap(data_traj)
[~,data_order_index] = sort(sum(data_traj,2),"descend");


figure(1)

color_limits = [-0.001,0.3];
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