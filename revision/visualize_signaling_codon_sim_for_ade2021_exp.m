clear all
load('./example_data_format/mutual_info_cal_data_example.mat')
savepath = '../SubFigures2023/revision_2024/';
%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
data_savepath = '../raw_data2023/2024/';
load(strcat(data_savepath,'Sim18_wt_IkBo_codon_metric_wt_revision_1115.mat'),'data','metrics','collect_feature_vects');

collect_feature_vects_wt = collect_feature_vects;
data_wt = data;
metrics_wt = metrics;

load(strcat(data_savepath,'Sim18_wt_IkBo_codon_metric_ikbamut_revision_1115.mat'),'data','metrics','collect_feature_vects');

collect_feature_vects_mut = collect_feature_vects;
data_mut = data;
metrics_mut = metrics;
%index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
% need to be revised if change
index_IkBamut = [1,2,4];
index_WT = [1,2,4];

metrics_cal = cell(1,length(metrics));
metric_names = fieldnames(metrics{1});

i_data = 1;
for i_index_wt = 1:length(index_WT)
    data_info.info_ligand{i_data} = data_wt.info_ligand{index_WT(i_index_wt)};
    data_info.info_dose_str{i_data} = data_wt.info_dose_str{index_WT(i_index_wt)};
    data_info.data_label{i_data} = 'wt';
    % metrics_cal{i_data} = metrics_wt{index_WT(i_index_wt)};

    for i_metric_name = 1:length(metric_names)
        metrics_cal{i_data}.(metric_names{i_metric_name}) = metrics_wt{index_WT(i_index_wt)}.(metric_names{i_metric_name})(1:9:end,:);
    end

    i_data = i_data+1;
end

for i_index_mut = 1:length(index_IkBamut)
    data_info.info_ligand{i_data} = data_mut.info_ligand{index_IkBamut(i_index_mut)};
    data_info.info_dose_str{i_data} = data_mut.info_dose_str{index_IkBamut(i_index_mut)};
    data_info.data_label{i_data} = 'IkBamut';
    % metrics_cal{i_data} = metrics_mut{index_IkBamut(i_index_mut)};
    for i_metric_name = 1:length(metric_names)
        metrics_cal{i_data}.(metric_names{i_metric_name}) = metrics_mut{index_IkBamut(i_index_mut)}.(metric_names{i_metric_name})(1:9:end,:);
    end
    i_data = i_data+1;
end

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled


data_save_path = './data_signaling_codons_Machine_learning_format/';

% Check if the directory exists
if ~exist(data_save_path, 'dir')
    % If the directory does not exist, create it
    mkdir(data_save_path);
end


index_WT = [1,2,3];
index_IkBamut = [4,5,6];

data_name = {'IkBamut','WT'};

%% % plot codon distributions for 3 ligands

if 1 % plot codon distributions for 3 ligands

    codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};


    for i_codon = 1:length(codon_list)
        figure(1)
        paperpos=[0,0,200,100]/1.4;
        papersize=[200 100]/1.4;
        draw_pos=[10,10,180,90]/1.4;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        wt_codon = collect_feature_vects.(codon_list{i_codon})([1,2,3]); % y{i} size should be cellnum by 1
        IkBamut_codon = collect_feature_vects.(codon_list{i_codon})([4,5,6]); % 1000 by 10


        std_cal_y = [];
        for i_y = 1:length(wt_codon)
            std_cal_y = [std_cal_y;wt_codon{i_y}]
        end

        if 1 % for better visualization
        switch codon_list{i_codon}
            case 'OscVsNonOsc'
                std_y = std(std_cal_y(:))/8000;
            case 'EarlyVsLate'
                std_y = std(std_cal_y(:))/2000;
            case 'TotalActivity'
                std_y = std(std_cal_y(:))/2000;
            otherwise
                std_y = std(std_cal_y(:))/3000;
        end

        end

        % subplot(1,length(vis_data_field),i_data_field)

        al_goodplot_pair_RMSD_diff_size(wt_codon,[1:3:7],0.5,ones(size(wt_codon,2),1)*[0 0 255]/255 ,'bilateral',[],std_y); %left
        al_goodplot_pair_RMSD_diff_size(IkBamut_codon,[2:3:8],0.5,ones(size(IkBamut_codon,2),1)*[255 0 0]/255,'bilateral',[],std_y);

        xlim([0 9])

        xticks([0:9])
        xticklabels({})
        %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

        switch codon_list{i_codon}
            case 'OscVsNonOsc'
                yticks([0:0.2:0.4])
                ylim([-0.1,1.1]/2.5);
            case 'EarlyVsLate'
                yticks([0:0.4:0.8])
                ylim([-0.1,1.1]/1.25);
            case 'TotalActivity'
                yticks([0:0.6:1.2])
       %  yticklabels({'','0','','0.2','','0.4','','0.6','','0.8','','1.0',''})

                ylim([-0.1,1.1]/0.8);
            otherwise
                ylim([-0.1,1.1]);
        end

        % for i_x = 1:10
        %     plot([i_x,i_x],[0,5],'--','Color','k');hold on
        % end
        set(gca,'fontsize',11,'fontname','Arial');
        %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

        saveas(gcf,strcat(savepath,'Codon_',codon_list{i_codon},'_wt_IkBmut_signaling_codon_3ligand_sim'),'pdf');
        close
    end
end


%% plot heatmaps for trajectories
if 0
for i_data_file = 1:6 %length(data_file_list)

    %[~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
    color_lmt = [0,0.25];
    data_plot = metrics{i_data_file}.time_series;
    data_first_half = data_plot(:,1:size(data_plot,2)/2);
    data_first_half(isnan(data_first_half)) = 0;
    [~,data_order_index] = sort(sum(data_first_half,2),"descend");

    fig_gcf = plot_heatmap(data_plot,  data_order_index,  color_lmt);
    saveas(gcf,strcat(savepath,'Sim_for_Ade_wtvsIkBmut_',collect_feature_vects.info_data_type{i_data_file},...
        '_',collect_feature_vects.info_ligand{i_data_file}),'epsc');
    close

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