clear all
load('./example_data_format/mutual_info_cal_data_example.mat')

%nfkb_eg = nfkb;
data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';
data_savepath = '../raw_data2023/2024/';
savepath = '../SubFigures2023/revision_2024/';
load(strcat(data_savepath,'Ade2021_IkBamut_wt_8hr_codon_metric.mat'),'data','metrics','collect_feature_vects');

%index order: TNF-0, Pam-2, CpG-3, LPS-4,PolyIC-5
% need to be revised if change
index_IkBamut = [1:3];
index_WT = [4:6];


data_save_path = './data_signaling_codons_Machine_learning_format/';

% Check if the directory exists
if ~exist(data_save_path, 'dir')
    % If the directory does not exist, create it
    mkdir(data_save_path);
end

codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

data_name = {'IkBamut','WT'};
index_vec = {index_IkBamut,index_WT};


if 1 % plot codon distributions for 2 ligands
    for i_codon = 1:length(codon_list)
        figure(1)
        paperpos=[0,0,200,100]/1.4;
        papersize=[200 100]/1.4;
        draw_pos=[10,10,180,90]/1.4;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        y = collect_feature_vects.(codon_list{i_codon})([3,1,2]); % y{i} size should be cellnum by 1
        z = collect_feature_vects.(codon_list{i_codon})([6,4,5]); % 1000 by 10
        std_cal_y = [];
        for i_y = 1:length(y)
            std_cal_y = [std_cal_y;y{i_y}]
        end
        std_y = std(std_cal_y(:))/2000;

        % subplot(1,length(vis_data_field),i_data_field)

        al_goodplot_pair_RMSD_diff_size(z,[1:3:7],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'bilateral',[],std_y); %left
        al_goodplot_pair_RMSD_diff_size(y,[2:3:8],0.5,ones(size(z,2),1)*[255 0 0]/255,'bilateral',[],std_y);


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

        saveas(gcf,strcat(savepath,'Codon_',codon_list{i_codon},'_wt_IkBmut_signaling_codon_3ligand'),'pdf');
        close
    end
end
