function [] = draw_single_cell_confusion_revision_20241027()
% For Guo et al. Figure 7 S7: single-cell stimulus response specificity
%
% tested 05/15/2024, Matlab 2020a

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/2024/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/revision_2024/';
savepath = fig_save_path;
addpath('./lib/')
addpath('./src/')
addpath('./bin/')
addpath('./revision/')

%% draw 5 ligand single ligand stim

vers = '_r1';
vers_savefig = strcat(vers,'_matching_SS_p25x');%_IkBao _IkBao_onlyPeak ,'_SS_0x' _p1x
if 1 %draw 5 ligand single ligand stim
    % Sim16_IkBao_5_signle_ligand_codon_metric _r1?
    % Sim15_5_signle_ligand_codon_metric _r2
    % Sim16_IkBao_5_signle_ligand_codon_metric_0x r1
    % Sim16_IkBao_5_signle_ligand_codon_metric_p01x
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_0x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p1x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p01x _r1
    load(strcat('../raw_data2023/Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x',vers,'.mat'))
    data_IkBo = data;
    para_mat_IkBo = sim_data_tbl.parameter_value(1:9:end,:);
    ligand_index = [1,5,3,2,4]; % reoder the stim

    %% draw codons

    % codon_mat = double;
    single_cell_specificity_plot = 1;

    for i_data = 1:length(data.model_sim)

        data_info.info_ligand{i_data} = data.info_ligand{i_data};
        data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
        data_info.data_label{i_data} = 'SamplingComb';

    end
    metrics_cal = cell(1,length(metrics));
    metric_names = fieldnames(metrics{1});

    for i_metric_index = 1:length(metrics)
        for i_metric_name = 1:length(metric_names)
            metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
        end
    end

    length_data = length(data.model_sim);
    length_metrics = length(metrics);

    load(strcat('../raw_data2023/Sim8_5_signle_ligand_codon_metric_r3','.mat'))
    para_mat_wt = sim_data_tbl.parameter_value(1:9:end,:);

    for i_data = 1:length(data.model_sim)

        data_info.info_ligand{i_data+length_data} = data.info_ligand{i_data};
        data_info.info_dose_str{i_data+length_data} = data.info_dose_str{i_data};
        data_info.data_label{i_data+length_data} = 'SamplingComb';

    end

    for i_metric_index = 1:length(metrics)
        for i_metric_name = 1:length(metric_names)
            metrics_cal{i_metric_index+length_metrics}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
        end
    end


    %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
    [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled

    codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
    weights_codon = [1,1,1,1,1,1];
    %codon_list = {'TotalActivity'};

    %%

    clear codon_single_cell_wt codon_single_cell_SS
    codon_single_cell_wt=cell(2,1);
    codon_single_cell_IkBo=cell(2,1);

    for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
        for i_ligand = 1:length_data
            for i_codon = 1:length(codon_list)
                codon_single_cell_IkBo{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
            end
        end
    end

    for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
        for i_ligand = 1:length_data
            for i_codon = 1:length(codon_list)
                codon_single_cell_wt{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand+length_data}(i_cell,:);
            end
        end
    end


    if 1
        dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
        dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));

        %%

        % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
        %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
        %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
        %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]

        %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
        %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
        %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]

        dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
        dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);


        if 1 % revision figure [tested] 11/16 matlab 2024

            clear ligand_pairs_vec
            i_ligand_pairs = 1;
            for i_sti_1 = 1:length(ligand_index)
                for i_sti_2 = (i_sti_1+1):length(ligand_index)
                    ligand_pairs_vec{i_ligand_pairs} = [i_sti_1,i_sti_2];
                    i_ligand_pairs = i_ligand_pairs+1;
                end
            end

            ligand_pairs_vec = ligand_pairs_vec([1,2,4,3,5,7,9,6,8,10]);

            for i_codon = 1:length(codon_list)
                clear ligand_pair_single_cell_codon_diff
                for i_ligand_pairs = 1:length(ligand_pairs_vec)
                    for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                        i_ligand = ligand_pairs_vec{i_ligand_pairs}(1);
                        j_ligand = ligand_pairs_vec{i_ligand_pairs}(2);
                        ligand_pair_single_cell_codon_diff(i_cell,i_ligand_pairs) = codon_single_cell_wt{i_cell}(i_ligand,i_codon) - codon_single_cell_wt{i_cell}(j_ligand,i_codon);
                        ligand_pair_single_cell_codon_diff_IkBamut(i_cell,i_ligand_pairs) = codon_single_cell_IkBo{i_cell}(i_ligand,i_codon) - codon_single_cell_IkBo{i_cell}(j_ligand,i_codon);
                    end
                end
                ligand_pair_codon_diff{i_codon} = ligand_pair_single_cell_codon_diff;
                ligand_pair_codon_diff_IkBamut{i_codon} = ligand_pair_single_cell_codon_diff_IkBamut;

            end

            for i_codon = 1:length(codon_list)
                figure(1)
                paperpos=[0,0,350,100]*1.5;
                papersize=[350 100]*1.5;
                draw_pos=[10,10,330,90]*1.5;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

                y = abs(ligand_pair_codon_diff{i_codon});
                z = abs(ligand_pair_codon_diff{i_codon});
                % subplot(1,length(vis_data_field),i_data_field)

                al_goodplot_pair_RMSD(y,[],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'left',[],std(y(:))/3500);
                al_goodplot_pair_RMSD(z,[],0.5,ones(size(z,2),1)*[0 0 255]/255,'right',[],std(z(:))/3500);

                xlim([0.4 10.6])

                xticks([1:10])
                xticklabels({})
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

                ylim([0,5]);
                for i_x = 1:10
                    plot([i_x,i_x],[0,5],'--','Color','k');hold on
                end
                set(gca,'fontsize',14,'fontname','Arial');
                %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

                saveas(gcf,strcat(fig_save_path,'Codon_',codon_list{i_codon},'_codon_diff_distrib_wt',vers_savefig),'epsc');

                close
            end
        end

        if 0  % revision figure, test confusion pairs integral smaller than specific cells?
            % ligand_index = [1,5,3,2,4]

            dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));


            % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
            %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
            %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
            %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]

            %TNF-P3C, TNF-CpG, TNF-LPS, TNF-PIC [4,2,1,3]
            %CpG-P3C, LPS-P3C, P3C-PIC,   [9,7,10]
            %LPS-CpG, CpG-PIC, LPS-PIC   [5,8,6]

            dist_mat_wt = dist_mat_wt(:,[4,2,1,3,9,7,10,5,8,6]);

            conf_cells = dist_mat_wt < 1;

            clear ligand_pairs_vec
            i_ligand_pairs = 1;
            for i_sti_1 = 1:length(ligand_index)
                for i_sti_2 = (i_sti_1+1):length(ligand_index)
                    ligand_pairs_vec{i_ligand_pairs} = [i_sti_1,i_sti_2];
                    i_ligand_pairs = i_ligand_pairs+1;
                end
            end

            ligand_pairs_vec = ligand_pairs_vec([1,2,4,3,5,7,9,6,8,10]);


            i_codon_integral = 4;
            i_codon_peak = 2;
            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:5
                    cell_integral(i_cell,i_ligand) = codon_single_cell_wt{i_cell}(i_ligand,i_codon_integral);
                    cell_peak(i_cell,i_ligand) = codon_single_cell_wt{i_cell}(i_ligand,i_codon_peak);
                end
            end


            for i_ligand_pairs = 1:length(ligand_pairs_vec)
                figure(1)% integral
                paperpos=[0,0,80,100]*1.5;
                papersize=[80 100]*1.5;
                draw_pos=[10,10,60,90]*1.5;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                %conf_cells(i_ligand_pairs)

                y = abs(max(cell_integral(~conf_cells(:,i_ligand_pairs),ligand_pairs_vec{i_ligand_pairs}),[],2));
                
                z = abs(max(cell_integral(conf_cells(:,i_ligand_pairs),ligand_pairs_vec{i_ligand_pairs}),[],2));
                % subplot(1,length(vis_data_field),i_data_field)

                al_goodplot_pair_RMSD(y,[],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'left',[],std(y(:))/2500);
                al_goodplot_pair_RMSD(z,[],0.5,ones(size(z,2),1)*[255 0 0]/255,'right',[],std(z(:))/2500);

                xlim([0.4 1.6])

                xticks([1])
                xticklabels({})
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

                ylim([0,5]);
                for i_x = 1:10
                    plot([i_x,i_x],[0,2],'--','Color','k');hold on
                end
                set(gca,'fontsize',14,'fontname','Arial');
                %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

                saveas(gcf,strcat(fig_save_path,'Pair',num2str(i_ligand_pairs),'_conf_vs_spec_integral_wt',vers_savefig),'epsc');

                close

            end

            for i_ligand_pairs = 1:length(ligand_pairs_vec)
                figure(1)% integral
                paperpos=[0,0,80,100]*1.5;
                papersize=[80 100]*1.5;
                draw_pos=[10,10,60,90]*1.5;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                %conf_cells(i_ligand_pairs)

                y = abs(max(cell_peak(~conf_cells(:,i_ligand_pairs),ligand_pairs_vec{i_ligand_pairs}),[],2));
                
                z = abs(max(cell_peak(conf_cells(:,i_ligand_pairs),ligand_pairs_vec{i_ligand_pairs}),[],2));
                % subplot(1,length(vis_data_field),i_data_field)

                al_goodplot_pair_RMSD(y,[],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'left',[],std(y(:))/2500);
                al_goodplot_pair_RMSD(z,[],0.5,ones(size(z,2),1)*[255 0 0]/255,'right',[],std(z(:))/2500);

                xlim([0.4 1.6])

                xticks([1])
                xticklabels({})
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

                ylim([0,5]);
                for i_x = 1:10
                    plot([i_x,i_x],[0,2],'--','Color','k');hold on
                end
                set(gca,'fontsize',14,'fontname','Arial');
                %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

                saveas(gcf,strcat(fig_save_path,'Pair',num2str(i_ligand_pairs),'_conf_vs_spec_peak_wt',vers_savefig),'epsc');

                close

            end

        end
    end



end

end


%%
function conf_mat = confusion_mat(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3
conf_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        conf_mat(i_col) = (sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2))<epsilon);
        i_col = i_col + 1;
    end
end

end

%%
function dist_mat = dist_mat(matrix)

% epsilon = 2; % Set the epsilon value 1,2,3
dist_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        dist_mat(i_col) = sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2));
        i_col = i_col + 1;
    end
end

end



%%
function [fig1,fig2,fig3,Outperm_col] = hierarchical_row_column_plot_kw_0411(traj_mat,color_limit,color_map,conf_pairs,Outperm_col)

traj_mat = traj_mat>1;

figure(1)

% Step 1: Perform hierarchical clustering on the rows
Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(1)
paperpos=[0,0,55,70];
papersize=[55 70];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig1 = gcf;

figure(2)
% Step 1: Perform hierarchical clustering on the rows
Y_col = pdist(traj_mat', 'euclidean'); % Compute the pairwise distances between rows
Z_col = linkage(Y_col, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm_col0] = dendrogram(Z_col, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(2)
paperpos=[0,0,70,55];
papersize=[70 55];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h_col,~,~] = dendrogram(Z_col,0, 'Orientation', 'top'); % Plotting dendrogram separately
% set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h_col,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig2 = gcf;

if nargin <5
    Outperm_col = Outperm_col0;
end

% Step 3: Reorder the matrix based on the clustering result
reordered_traj_mat = traj_mat(Outperm, Outperm_col);



figure(3)

paperpos=[0,0,100,64]*3;
papersize=[100,64]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

% subplot(1,length(vis_data_field),i_data_field)
if nargin <2
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,1]);%[-0.001,0.2] for TNF
elseif nargin <3
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
else
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',color_map,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF

end
%
XLabels = 1:size(reordered_traj_mat,2);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/1);%conf_pairs
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

conf_pairs = conf_pairs(Outperm_col);


CustomXLabels = string(conf_pairs);

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_mat,1);
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
fig3 = gcf;


end


%%
function [fig1,fig2,fig3,Outperm_col] = hierarchical_row_column_plot_diff_thresh_kw_0618(traj_mat,threshold,color_limit,color_map,conf_pairs,Outperm_col)

traj_mat = traj_mat>threshold;

figure(1)

% Step 1: Perform hierarchical clustering on the rows
Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(1)
paperpos=[0,0,55,70];
papersize=[55 70];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig1 = gcf;

figure(2)
% Step 1: Perform hierarchical clustering on the rows
Y_col = pdist(traj_mat', 'euclidean'); % Compute the pairwise distances between rows
Z_col = linkage(Y_col, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm_col0] = dendrogram(Z_col, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(2)
paperpos=[0,0,70,55];
papersize=[70 55];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h_col,~,~] = dendrogram(Z_col,0, 'Orientation', 'top'); % Plotting dendrogram separately
% set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h_col,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig2 = gcf;

if nargin <6
    Outperm_col = Outperm_col0;
end

% Step 3: Reorder the matrix based on the clustering result
reordered_traj_mat = traj_mat(Outperm, Outperm_col);



figure(3)

paperpos=[0,0,100,64]*3;
papersize=[100,64]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

% subplot(1,length(vis_data_field),i_data_field)
if nargin <3
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,1]);%[-0.001,0.2] for TNF
elseif nargin <4
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
else
    h=heatmap(double(reordered_traj_mat(:,:)),'ColorMap',color_map,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF

end
%
XLabels = 1:size(reordered_traj_mat,2);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/1);%conf_pairs
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

conf_pairs = conf_pairs(Outperm_col);


CustomXLabels = string(conf_pairs);

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_mat,1);
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
fig3 = gcf;


end