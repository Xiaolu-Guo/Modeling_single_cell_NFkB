% For Guo et al. Figure S3A, Heatmap vs parameters
% 
% tested 05/12/2024, Matlab 2020a

% read the predicted data for all. and then order the heatmap based on fano
% facotr / CV or oscillation power. and then draw the plots of parameter
% distribution.

data_save_file_path = '../raw_data2023/';%_fay_parameter/';

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))% data,collect_feature_vects,metrics));
%All_ligand_codon_2023_t33_cv_filtered_TNF.mat
% data.parameters_mode_nan
% data.exp_mode_filter_nan
% data.pred_mode_filter_nan
% are the TNF cv filtered trajectories

fig_save_path = '../SubFigures2023/';

paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

paperpos_para = [0,0,50,130]*3;
papersize_para = [50 130]*3;
draw_pos_para = [10,10,40,120]*3;

features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

data_order=struct();
for i_feature = 1:length(features)
    for i_sti = 1:length(data.pred_mode_filter_nan)
        
        [~,data_order.(features{i_feature}){i_sti}] = sort(collect_feature_vects.(features{i_feature}){i_sti*2},'descend');
        
    end
end
vis_data_field = {'pred_mode_filter_nan'};


ligand_high_index = [3,6,12,16,19];

for i_feature = 1:length(features)
    if ~isfolder(strcat(fig_save_path,features{i_feature},'/'))
        mkdir(strcat(fig_save_path,features{i_feature},'/'));
    end
    for i_sti = ligand_high_index %1:length(data.(vis_data_field{1}))
        for i_data_field =1:length(vis_data_field)
            % Figure S3A heatmaps
            figure(i_sti)
            
            cell_num = length(data_order.(features{i_feature}){i_sti});
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            % subplot(1,length(vis_data_field),i_data_field)
            h=heatmap(data.(vis_data_field{i_data_field}){i_sti}(data_order.(features{i_feature}){i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
            
            XLabels = 0:5:((size(data.(vis_data_field{i_data_field}){i_sti},2)-1)*5);
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
            
            %         xlabel('Time (hours)');
            %         ylabel(vis_data_field{i_data_field});
            % clb=colorbar;
            %
            % clb.Label.String = 'NFkB(A.U.)';
            colorbar('off')
            
            
            set(gca,'fontsize',14,'fontname','Arial');
            saveas(gcf,strcat(fig_save_path,features{i_feature},'/','heatmap_sim_orderedby_',data.info_ligand{i_sti},'_',replace(data.info_dose_str{i_sti},'/','')),'epsc');
            close
            if 1
                for i_para = 1:(size(data.parameters_mode_nan{i_sti},2)-1)
                    % Figure S3A scatter plots
                    figure(i_sti*2)
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos_para,'PaperSize', papersize_para,'Position',draw_pos_para)
                    
                    scatter(data.parameters_mode_nan{i_sti}(data_order.(features{i_feature}){i_sti},i_para),data_order.(features{i_feature}){i_sti},10,'filled');
                    ylim([min(data_order.(features{i_feature}){i_sti}),max(data_order.(features{i_feature}){i_sti})]);
                    set(gca,'XTickLabel',{},'YTickLabel',{},'XScale','log','YDir','reverse');
                    
                    
                    saveas(gcf,strcat(fig_save_path,features{i_feature},'/',data.info_ligand{i_sti},'para_sim_orderedby_',replace(data.info_dose_str{i_sti},'/',''),...
                        '_',data.para_name{i_sti}{i_para}),'epsc');
                    close
                end
            end
        end
    end
    
end




