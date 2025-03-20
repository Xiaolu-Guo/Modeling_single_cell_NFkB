addpath('../raw_data2023/2024/')
savepath = '../SubFigures2023/revision_2024/';
data_savepath = '../raw_data2023/2024/';
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

i_data_file = 16; %22 for others; 16 for IkBa mut exp
data_traj = readmatrix(strcat(data_file_list{i_data_file},'.csv'));
rescale_factor = rescale_factor_cal(data_traj);


%% [need to be tested] visualize 3 doses for WT and IkBamut
if 1
cal_codon = 1;
if cal_codon
    i_codon_cal_file = 1;
    for i_data_file = 13:18 %length(data_file_list)

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
    save(strcat(data_savepath,'Ade2021_IkBamut_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');
else
    %visualize ade's dataset for mutant
    load(strcat(data_savepath,'Ade2021_IkBamut_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');
end

i_codon_cal_file =  1;
color_lmt = [0,0.2];
for i_data_file = 13:18 %length(data_file_list)


    if 0
        [~,order_index] = sort(metrics{i_codon_cal_file}.oscpower,"descend");
        fig_gcf = plot_heatmap(data.exp{i_codon_cal_file},order_index,color_lmt);
        saveas(gcf,strcat(savepath,data_file_list{i_data_file},'_orderby_oscpower'),'epsc');
        close
    end

    if 0
        [~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
        fig_gcf = plot_heatmap(data.exp{i_codon_cal_file}, order_index,color_lmt);
        saveas(gcf,strcat(savepath,data_file_list{i_data_file},'_orderby_oscfreq'),'epsc');
        close
    end

    %[~,order_index] = sort(metrics{i_codon_cal_file}.oscfreq,"descend");
    color_lmt = [0,0.25];
    data_first_half = data.exp{i_codon_cal_file}(:,1:size(data.exp{i_codon_cal_file},2)/2);
    data_first_half(isnan(data_first_half)) = 0;
    [~,data_order_index] = sort(sum(data_first_half,2),"descend");

    fig_gcf = plot_heatmap(data.exp{i_codon_cal_file},  data_order_index,  color_lmt);
    saveas(gcf,strcat(savepath,data_file_list{i_data_file}),'epsc');
    close

    i_codon_cal_file = i_codon_cal_file +1;

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