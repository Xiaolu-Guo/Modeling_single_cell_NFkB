function [] = draw_scatter_codon_para_2023_06(data_save_file_path,fig_save_path)
%P3C: speed vs Tak; duration cdeg;
%PIC: OSC vs timed1; EvL vs Cdeg; Peak vs rcp syn
%LPS: tot vs NFkB tot
%
std_thresh = 0.1;
% load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

ligand_vec = {'LPS','PolyIC','Pam3CSK'};% 'TNF','CpG',
% dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
%     {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
%     {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
%     {'10ug/mL';'33ug/mL';'100ug/mL'};
%     {'10ng/mL';'100ng/mL';'1ug/mL'}};
% data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

% fig_save_path = strcat(fig_save_path,'2023_03_');
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];


collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

ligand_vec = {'LPS','PolyIC','Pam3CSK','TNF'};% 'TNF','CpG',
codon_index = {[4],[6,5,2],[1,3],[6]};
para_index = {[7],[2,5,4],[1,5],[5]};% pay attention to NFkB tot!!!
negative_corr = {[0],[0,0,0],[0,1],[1]};
para_str_name = {'TAK1 activation','Time delay','Time delay','Receptor synthesis','Complex degradation','Traffic to Endosome','NFkB tot'};
xticks = {{[0.05,0.1,0.2]},{[1e-2,1e-1,1],[1e-4,1e-3,1e-2],[1e-7,1e-6,1e-5]},{[0.1,1,10],[0.001,0.01,0.1]},{[1e-2,1e-1,1]}};
x_lims = {{[0.03,0.3]},{[9e-3,1],[6e-5,0.01],[1e-7,1e-5]},{[1e-1,10],[5e-4,0.1]},{[1e-2,1.01]}};
%%
for i_sti = 1:length(ligand_vec)
    ligand = ligand_vec{i_sti};
    dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
    for i_dose = length(dose_vec)
        dose  = dose_vec(i_dose);
        para_vals = data.parameters_mode_nan{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose};
        para_vals = para_vals(:,1:end-1);
        sig_codons_model = [collect_feature_vects.(features{1}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{2}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{3}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{4}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{5}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{6}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'}];
        sig_codons_exp = [collect_feature_vects.(features{1}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'},...
            collect_feature_vects.(features{2}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'},...
            collect_feature_vects.(features{3}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'},...
            collect_feature_vects.(features{4}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'},...
            collect_feature_vects.(features{5}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'},...
            collect_feature_vects.(features{6}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'experiments'}];
        
        std_exp = std(sig_codons_exp);
        std_model = std(sig_codons_model);
        mean_exp = mean(sig_codons_exp);
        mean_model = mean(sig_codons_model);
        %% core para
        corr_mat =  corr(para_vals(:,[1:3,size(para_vals,2)]),...
            sig_codons_model,...
            'Type','Spearman');
        index = find(abs(mean_model-mean_exp)>0.1 | abs(std_model-std_exp)>0.15 );
        index_non_wide = find(std_exp<std_thresh);
        
        %         corr_mat(:,abs(mean_model-mean_exp)>0.1) = NaN;
        %         corr_mat(:,abs(std_model-std_exp)>0.15) = NaN;
        %         % comment
        %         corr_mat(:,std_exp<std_thresh) = NaN;
        
        para_names_str = data.para_name{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose}(:,1:end-1);
        para_names_str{end} = 'NFKBInit';
        
        
        %% receptor para
        rcp_mat =  corr(para_vals(:,setdiff(1:size(para_vals,2),[1:3,size(para_vals,2)],'stable')),...
            [collect_feature_vects.(features{1}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{2}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{3}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{4}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{5}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{6}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'}],...
            'Type','Spearman');
        
        
        %         rcp_mat(:,abs(mean_model-mean_exp)>0.1) = NaN;
        %         rcp_mat(:,abs(std_model-std_exp)>0.15) = NaN;
        %         % comment
        %         rcp_mat(:,std_exp<std_thresh) = NaN;
        if 1
            for i_ligand_codon = 1:length(codon_index{i_sti})
                figure(3)
                
                set(gcf, 'PaperUnits','points')
                
                paper_pos = [0,0,120,80];
                paper_size = [120,80];
                set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
                
                para_ind = para_index{i_sti}(i_ligand_codon);
                codon_ind = codon_index{i_sti}(i_ligand_codon);
                
                scatter(para_vals(:,para_ind),sig_codons_model(:,codon_ind),2,'filled');hold on
                x = para_vals(:,para_ind);
                y = sig_codons_model(:,codon_ind);
                if negative_corr{i_sti}(i_ligand_codon)
                    
                    ft = fittype( 'Hill_line_neg(x,n,K,A)' );
                    f = fit( x, y, ft, 'StartPoint', [1, mean(x), 1] );
                else
                    ft = fittype( 'Hill_line_pos(x,n,K,A)' );
                    f = fit( x, y, ft, 'StartPoint', [2, mean(x), 1] );
                    
                end
                
                if (i_sti ==2) && (i_ligand_codon == 3)
                    f = fit( x, y, ft, 'StartPoint', [1, mean(x), 1] );% polyIC peakamp rcp syn
                end
                
                % f=fit(x,y,'smoothingspline');%  lowess
                h = plot(f,x,y);
                h(2).LineWidth = 1;
                hLeg = legend();
                set(hLeg,'visible','off')
 
                xlim(x_lims{i_sti}{i_ligand_codon});
                set(gca,'XScale','log','FontSize',7,'XColor',[0,0,0],'YColor',[0,0,0],'XTick',xticks{i_sti}{i_ligand_codon})
                % xlabel('TAK1 activiation')%TAK1 activiation
                % xlabel('Time delay')
                xlabel(para_str_name{para_ind},'FontWeight','b','FontSize',9,'Color',[1,1,1])%TAK1 activiation
                % xlabel('Traffic to Endosome')%TAK1 activiation
                % xlabel('TLR synthesis')%TAK1 activiation
                
                % ylabel('Oscillation')
                ylabel(features{codon_ind},'FontWeight','b','FontSize',9)
                saveas(gcf,strcat(fig_save_path,ligand,'_',features{codon_ind},'_',para_str_name{para_ind}),'epsc')
                close
                % ylabel('Speed')
            end
        end
    end
end
