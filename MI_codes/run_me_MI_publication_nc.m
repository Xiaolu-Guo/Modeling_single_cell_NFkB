% run under the folder '/Users/admin/Documents/my document/Postdoc
% projects/MatlabCodes/NFkB_Mutual_Information/'
% please change the path to your own local file folder path

addpath(genpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information/')); % please change the path to your own local file folder path


%% [MI calculation] doses study
% channel capacity for codon 20 doses vs exp 3 doses + zero dose

if 0
    
    clear all
    addpath(genpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information/')); % please change the path to your own local file folder path
    
    
    local_flag = 0;
    ligand_all = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
    for i_ligand = 1:length(ligand_all)
        ligand = ligand_all{i_ligand};
        data_file = {strcat('mutual_info_format_codon_',ligand,'_3doses_zerodose_20230707.mat')};
        %                 data_file = {strcat('mutual_info_format_codon_',ligand,'_3doses_zerodose_20230707.mat');
        %             strcat('mutual_info_format_codon_',ligand,'_20doses_20230614.mat')};
        
        data_label = {'exp-doses'};
        %                 data_label = {'distinct-3-doses','21-doses'};
        
        if 1 % calculate mutual info or load mutual info mat
            
            for i_data = 1:length(data_file)
                [~,sc_info{i_data}] = run_script_fxn_signaling_codons_AS_XG(local_flag, data_file{i_data});
            end
        end
        
        for i_data = 1:length(data_file)
            info_mat(i_data) = sc_info{i_data}.I;
        end
        
        if 0
            fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
            figure(1)
            paperpos = [0,0,120,80];
            papersize = [120,80];
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
            c = categorical(data_label);
            c = reordercats(c,data_label);
            bar_pct=info_mat;
            b= bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);
            ax2= gca;
            % ytickformat(ax2, '%g%%');
            % ylim([75,100])
            ylabel({'Chanel Capacity'})
            
            % title(rmsd_cas)
            set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
            saveas(gcf,strcat(fig_save_path,ligand,'_chanel_capacity_codon_3_21doses'),'epsc')
            close
        end
    end
end

if 0 % experimental data with unstim
    
    clear all
    addpath(genpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information'));
    
    local_flag = 0;
    ligand_all = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
    for i_ligand = 1:length(ligand_all)
        ligand = ligand_all{i_ligand}
        data_file = {strcat('mutual_info_format_ade_exp_',ligand,'.mat')};
        %                 data_file = {strcat('mutual_info_format_codon_',ligand,'_3doses_zerodose_20230707.mat');
        %             strcat('mutual_info_format_codon_',ligand,'_20doses_20230614.mat')};
        
        data_label = {'exp-data-3doses'};
        %                 data_label = {'distinct-3-doses','21-doses'};
        
        if 1 % calculate mutual info or load mutual info mat
            
            for i_data = 1:length(data_file)
                [~,sc_info{i_data}] = run_script_fxn_signaling_codons_AS_XG(local_flag, data_file{i_data});
            end
        end
        
        for i_data = 1:length(data_file)
            info_mat(i_data) = sc_info{i_data}.I
        end
        
        
    end
end

if 0
    ligand_all = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
    info_all ={[0.9594 1.7288, 2.0622];
        [1.2022 1.2095 ,1.2668 ];
        [0.8034 0.95846 ,0.99932];
        [0.8971 1.1914 ,1.4514];% (t = 99.1 sec)
        [0.8544 0.94212  1.0262]};
    
    for i_ligand = 1:length(ligand_all)
        ligand = ligand_all{i_ligand};
        info_mat = info_all{i_ligand};
        
        %data_label = {'exp-data-3doses','sim-data-3doses','sim-data-20doses'};
        
        if 1
            fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
            figure(1)
            paperpos = [0,0,150,80]/1.5;
            papersize = [150,80]/1.5;
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
            %c = categorical(data_label);
            %c = reordercats(c,data_label);
            bar_pct=info_mat;
            %b= bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);
            b= bar( bar_pct, 0.5,'EdgeColor',[0 0 0],'LineWidth',0.5);
            
            ax2= gca;
            % ytickformat(ax2, '%g%%');
            ylim([0,2.1])
            xlim([0.5,3.5])
            b(1).FaceColor = 'flat';
            
            b(1).CData = [[0.7,0.7,0.7];[0.7,0.7,0.7];[0.7,0.7,0.7]];
            
            set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial','XTickLabel',{'','',''})
            
            if i_ligand>0
                %set(gca,'fontsize',7,'XColor','k','FontName','Arial')%'YTickLabel',{},
                ylabel({'Chanel Capacity'},'Color','none')
            else
                ylabel({'Chanel Capacity'})
            end
            saveas(gcf,strcat(fig_save_path,ligand,'_chanel_capacity_codon_3_21doses'),'epsc')
            close
        end
        
    end
end


%% [MI calculation] denoise netowrk
% calculate MI for all red, add noise, vs wt, IkBo 0219, Figure 5
if 1
    
    clear all
    addpath(genpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information'));
    
    local_flag = 0;
    
    min_max_rescale = 1;
    if min_max_rescale
        data_file = {'mutual_info_format_alldata_sc_NFkB_var_red1_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red2_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red3_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red4_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red5_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red6_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red7_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red8_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red9_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red10_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_wt1_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_IkBo2_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_IkBo3_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_IkBo4_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_rcpt_red_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_TAKac_red_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_TAKac_NFkB_red1_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_TAKac_NFkB_red2_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_TAK_noise_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_NFkB_noise_min_max_rescale_0219.mat',...
            'mutual_info_format_alldata_sc_no_noise_min_max_rescale_0219.mat'};
    else
        data_file =  {'mutual_info_format_alldata_sc_NFkB_var_red1.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red2.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red3.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red4.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red5.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red6.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red7.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red8.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red9.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red10.mat',...
            'mutual_info_format_alldata_sc_wt1.mat',...
            'mutual_info_format_alldata_sc_IkBo2.mat',...
            'mutual_info_format_alldata_sc_IkBo3.mat',...
            'mutual_info_format_alldata_sc_IkBo4.mat',...
            'mutual_info_format_alldata_sc_rcpt_red.mat',...
            'mutual_info_format_alldata_sc_TAKac_red.mat'};
    end
    
    data_label = {'NFkB-denoise1',...
        'NFkB-denoise2',...
        'NFkB-denoise3',...
        'NFkB-denoise4',...
        'NFkB-denoise5',...
        'NFkB-denoise6',...
        'NFkB-denoise7',...
        'NFkB-denoise8',...
        'NFkB-denoise9',...
        'NFkB-denoise10',...
        'wt1',...
        'IkBo2',...
        'IkBo3',...
        'IkBo4',...
        'rcp-denoise',...
        'TAKac-denoise',...
        'rcp-noise1',...
        'rcp-noise2',...
        'TAKac-noise',...
        'NFkB-noise',...
        'no-noise'};
    
    if 0 % calculate mutual info or load mutual info mat
        %97 is th time trajectory dimension
        info_mat =[];
        for i_rpt = 1:5
            for i_data = 1:length(data_file)
                data_label{i_data}
                [~,sc_info{i_data}] = run_script_fxn_traj_XG_2024(local_flag, data_file{i_data},97);% current code, 97 is useless
            end
            
            for i_data = 1:length(data_file)
                info_mat(i_data,i_rpt) = sc_info{i_data}.I;
            end
        end
        if min_max_rescale
            save('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0219.mat','data_label','info_mat')
        else
            save('MI_network_flow_alldata_sc_NFkBfeature.mat','data_label','info_mat')
            
        end
    else
        
        % I = 1.0397 for IkBao
        % I = 1.0873 for wt
        % info_mat = [1.0397,1.0873];
    end
    
    
    
    if 1 % add noise different part
        load('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0219.mat')
        %load('MI_network_flow_NFkBfeature.mat')
        
        fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
        figure(2)
        paperpos = [0,0,100,100];
        papersize = [100,100];
        
        bar_means = mean(info_mat, 2);
        bar_stddevs = std(info_mat, 0, 2);
        
        data_label{17} = 'rcp-noise';
        index_data = [21,17,19,20];
        bar_means_plot = bar_means(index_data);
        bar_std_plot = bar_stddevs(index_data);
        
        info_NFkBred_mat = info_mat(17:18,:);
        bar_means_plot(2) = mean(info_NFkBred_mat(:));
        bar_std_plot(2) = std(info_NFkBred_mat(:));
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label(index_data));
        c = reordercats(c,data_label(index_data));
        bar(c,bar_means_plot,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
        
        numGroups = size(bar_means_plot, 1);
        numBars = size(info_mat, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) ; % Adjust the position - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars)
        
        % Add error bars
        errorbar(x, bar_means_plot, bar_std_plot, 'k', 'linestyle', 'none');
        
        x_positions = 1:numGroups; % X positions of bars

        for i = 1:numGroups
            % Extract raw data points

            raw_data = info_mat(index_data(i), :); % Use corresponding data from info_mat


            jitter = linspace(-0.3,0.3,length(raw_data));%(rand(size(raw_data)) - 0.5) * 1; % Add jitter for better visibility
            scatter(x_positions(i) + jitter, raw_data, 3, 'k','filled'); % Plot raw data 'filled', 'MarkerFaceAlpha', 0.6
        end

        hold off;

        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        % ylabel({'Chanel Capacity'})
        xticklabels({});
        yticklabels({});
        
        % title(rmsd_cas)
        fig_save_path = './';
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
        saveas(gcf,strcat(fig_save_path,'CC_Netowrk_addnoise_0219'),'epsc') %_min_max_rescale
        close
    end
    
    if 0 % denoise different part
        load('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0219.mat')
        %load('MI_network_flow_NFkBfeature.mat')
        
        fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
        figure(2)
        paperpos = [0,0,100,100];
        papersize = [100,100];
        
        bar_means = mean(info_mat, 2);
        bar_stddevs = std(info_mat, 0, 2);
        
        data_label{17} = 'rcp-noise';
        index_data = [11,15,16,1];
        bar_means_plot = bar_means(index_data);
        bar_std_plot = bar_stddevs(index_data);
        
        %         info_NFkBred_mat = info_mat(17:18,:);
        %         bar_means_plot(2) = mean(info_NFkBred_mat(:));
        %         bar_std_plot(2) = std(info_NFkBred_mat(:));
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label(index_data));
        c = reordercats(c,data_label(index_data));
        bar(c,bar_means_plot,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
        
        numGroups = size(bar_means_plot, 1);
        numBars = size(info_mat, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars); % Adjust the position
        
        % Add error bars
        errorbar(x, bar_means_plot, bar_std_plot, 'k', 'linestyle', 'none');


        % Plot individual data points

        x_positions = 1:numGroups; % X positions of bars

        for i = 1:numGroups
            % Extract raw data points

            raw_data = info_mat(index_data(i), :); % Use corresponding data from info_mat


            jitter = linspace(-0.3,0.3,length(raw_data));%(rand(size(raw_data)) - 0.5) * 1; % Add jitter for better visibility
            scatter(x_positions(i) + jitter, raw_data, 3, 'k','filled'); % Plot raw data 'filled', 
        end

        hold off;
       
        
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        % ylabel({'Chanel Capacity'})
        xticklabels({});
        yticklabels({});
        
        % title(rmsd_cas)
        fig_save_path = './';
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
        saveas(gcf,strcat(fig_save_path,'CC_Netowrk_wt_denoise_0219'),'epsc') %_min_max_rescale
        close
    end
    
    if 0 % plot information flow
        load('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0219.mat')
        %load('MI_network_flow_NFkBfeature.mat')
        
        fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
        figure(2)
        paperpos = [0,0,100,100];
        papersize = [100,100];
        
        bar_means_all = mean(info_mat, 2);
        bar_stddevs_all = std(info_mat, 0, 2);
        
        clear bar_means bar_stddevs
        info_NFkBred_mat = info_mat(1:10);
        bar_means = [log2(5),bar_means_all(17),mean(info_NFkBred_mat(:)),bar_means_all(11)]; % rcpt & TAK1ac
        bar_stddevs = [0,bar_stddevs_all(17),std(info_NFkBred_mat(:)),bar_stddevs_all(11)];
        
        
        data_label{17} = 'rcp-noise';
        index_data = [21,17,1,11];
        bar_means = bar_means_all(index_data);
        bar_stddevs = bar_stddevs_all(index_data);
        
        
        info_NFkBred_mat = info_mat(17:18,:);
        bar_means(2) = mean(info_NFkBred_mat(:));
        bar_stddevs(2) = std(info_NFkBred_mat(:));
        
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label(index_data));
        c = reordercats(c,data_label(index_data));
        bar(c,bar_means,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
        
        numGroups = length(bar_means);
        numBars = size(info_mat, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars); % Adjust the position

        % Add error bars
        errorbar(x, bar_means, bar_stddevs, 'k', 'linestyle', 'none');



        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        xticklabels({});
        yticklabels({});
        % title(rmsd_cas)
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
        saveas(gcf,strcat(fig_save_path,'CC_Netowrk_level_MI_0219'),'epsc') %_min_max_rescale
        close
    end
    
    if 0 % supp
        load('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0219.mat')
        %load('MI_network_flow_NFkBfeature.mat')
        
        data_label = {'NFkB-denoise1',...
            'NFkB-denoise2',...
            'NFkB-denoise3',...
            'NFkB-denoise4',...
            'NFkB-denoise5',...
            'NFkB-denoise6',...
            'NFkB-denoise7',...
            'NFkB-denoise8',...
            'NFkB-denoise9',...
            'NFkB-denoise10',...
            'wt1',...
            'IkBo2',...
            'IkBo3',...
            'IkBo4',...
            'rcp-denoise',...
            'TAKac-denoise',...
            'rcp-noise1',...
            'rcp-noise2',...
            'TAKac-noise',...
            'NFkB-noise',...
            'no-noise'};
        
        fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
        figure(2)
        paperpos = [0,0,200,100]*1.5;
        papersize = [200,100]*1.5;
        
        clear bar_means bar_stddevs
        
        bar_means_all = mean(info_mat, 2);
        bar_stddevs_all = std(info_mat, 0, 2);
        
        data_label{17} = 'rcp-noise';
        index_data = [11,1:10];
        bar_means = bar_means_all(index_data);
        bar_stddevs = bar_stddevs_all(index_data);
        data_label = data_label(index_data);
        
        
        
        info_NFkBred_mat = info_mat(1:10);
        bar_means(12) = mean(info_NFkBred_mat(:));
        bar_stddevs(12) = std(info_NFkBred_mat(:));
        data_label{12} = 'avg-NFkB-denoise';
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label);
        c = reordercats(c,data_label);
        bar(c,bar_means,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
        
        numGroups = length(bar_means);
        numBars = size(info_mat, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars); % Adjust the position
        
        % Add error bars
        errorbar(x, bar_means, bar_stddevs, 'k', 'linestyle', 'none');
        
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        %         xticklabels({});
        %         yticklabels({});
        % title(rmsd_cas)
        ylabel('Chanel Capacity')
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
        saveas(gcf,strcat(fig_save_path,'CC_Netowrk_level_MI_0219'),'epsc') %_min_max_rescale
        close
    end
    
    
    
end


%% [MI calculation] wt (wild type) vs Ikbo (IkBa mutant)
% calculate MI for all red, add noise, vs wt, IkBo 0401, Figure 5
if 0
    
    clear all
    addpath(genpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information'));
    
    local_flag = 0;
    
    min_max_rescale = 1;
    if min_max_rescale
        data_file = {'mutual_info_format_alldata_sc_NFkB_var_red1_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red2_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red3_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red4_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red5_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red6_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red7_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red8_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red9_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red10_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_wt1_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_IkBo2_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_IkBo3_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_IkBop25_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_rcpt_red_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_TAKac_red_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_TAKac_NFkB_red1_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_TAKac_NFkB_red2_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_TAK_noise_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_NFkB_noise_min_max_rescale_0331.mat',...
            'mutual_info_format_alldata_sc_no_noise_min_max_rescale_0331.mat'};
    else
        data_file =  {'mutual_info_format_alldata_sc_NFkB_var_red1.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red2.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red3.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red4.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red5.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red6.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red7.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red8.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red9.mat',...
            'mutual_info_format_alldata_sc_NFkB_var_red10.mat',...
            'mutual_info_format_alldata_sc_wt1.mat',...
            'mutual_info_format_alldata_sc_IkBo2.mat',...
            'mutual_info_format_alldata_sc_IkBo3.mat',...
            'mutual_info_format_alldata_sc_IkBo4.mat',...
            'mutual_info_format_alldata_sc_rcpt_red.mat',...
            'mutual_info_format_alldata_sc_TAKac_red.mat'};
    end
    
    data_label = {'NFkB-denoise1',...
        'NFkB-denoise2',...
        'NFkB-denoise3',...
        'NFkB-denoise4',...
        'NFkB-denoise5',...
        'NFkB-denoise6',...
        'NFkB-denoise7',...
        'NFkB-denoise8',...
        'NFkB-denoise9',...
        'NFkB-denoise10',...
        'wt1',...
        'IkBo2',...
        'IkBo3',...
        'IkBop25',...
        'rcp-denoise',...
        'TAKac-denoise',...
        'rcp-noise1',...
        'rcp-noise2',...
        'TAKac-noise',...
        'NFkB-noise',...
        'no-noise'};
    
    if 1 % calculate mutual info or load mutual info mat
        %97 is th time trajectory dimension
        info_mat =[];
        for i_rpt = 1:5
            for i_data = 1:length(data_file)
                data_label{i_data}
                [~,sc_info{i_data}] = run_script_fxn_traj_XG_2024(local_flag, data_file{i_data},97);% current code, 97 is useless
            end
            
            for i_data = 1:length(data_file)
                info_mat(i_data,i_rpt) = sc_info{i_data}.I;
            end
        end
        if min_max_rescale
            save('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0401.mat','data_label','info_mat')
        else
            save('MI_network_flow_alldata_sc_NFkBfeature_0401.mat','data_label','info_mat')
            
        end
    else
        % I = 1.0397 for IkBao
        % I = 1.0873 for wt
        % info_mat = [1.0397,1.0873];
    end
    
    if 1 % WT vs IkBas
        load('MI_network_flow_alldata_sc_NFkBfeature_min_max_rescale_0401.mat')
        %load('MI_network_flow_NFkBfeature.mat')
        
        fig_save_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/SubFigures2023/';
        figure(2)
        paperpos = [0,0,100,100];
        papersize = [100,100];
        
        bar_means_all = mean(info_mat, 2);
        bar_stddevs_all = std(info_mat, 0, 2);
        
        %         clear bar_means bar_stddevs
        %         info_NFkBred_mat = info_mat(1:10);
        %         bar_means = [log2(5),bar_means_all(17),mean(info_NFkBred_mat(:)),bar_means_all(11)]; % rcpt & TAK1ac
        %         bar_stddevs = [0,bar_stddevs_all(17),std(info_NFkBred_mat(:)),bar_stddevs_all(11)];
        
        data_label{17} = 'rcp-noise';
        index_data = [11,14];
        bar_means = bar_means_all(index_data);
        bar_stddevs = bar_stddevs_all(index_data);
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label(index_data));
        c = reordercats(c,data_label(index_data));
        bar(c,bar_means,0.4,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
        
        numGroups = length(bar_means);
        numBars = size(info_mat, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars); % Adjust the position
        
        % Add error bars
        errorbar(x, bar_means, bar_stddevs, 'k', 'linestyle', 'none');
        
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        xticklabels({});
        yticklabels({});
        % title(rmsd_cas)
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
        saveas(gcf,strcat(fig_save_path,'CC_WT_IkBass_MI_0401'),'epsc') %_min_max_rescale
        close
    end
    
    
end


%% sub functions drwa_conf_mat


function fig = draw_conf_mat(conf_mat)
conf_mat_plot = NaN(5,5);
i_conf_mat = 1;

for i_ligand_index = 1:size(conf_mat_plot,1)
    for j_ligand_index = (i_ligand_index+1):size(conf_mat_plot,2)
        conf_mat_plot(i_ligand_index,j_ligand_index) = conf_mat(i_conf_mat);
        i_conf_mat = i_conf_mat+1;
    end
end

conf_mat_plot = conf_mat_plot';
conf_mat_plot = conf_mat_plot;

figure(1)
paperpos=[0,0,64,64]*3;

papersize=[64,64]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)


h = heatmap(conf_mat_plot,'Colormap',parula,'GridVisible','off','MissingDataColor',[1 1 1], 'CellLabelFormat', '%.2f');
XLabels = 1:5;
CustomXLabels = string(XLabels);
CustomXLabels(:) = " ";
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
caxis([0,0.8]);
colorbar off
fig = gcf;
end