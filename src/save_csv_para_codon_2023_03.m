
function [] = save_csv_para_codon_2023_03(data_save_file_path,data_save_path,monolix_data_save_file_path)

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'10ug/mL';'33ug/mL';'100ug/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'}};
dose_symbol = {'Low','Med','High'};
data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

ligand_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};
proj_num_all = [2,3,4,5,6];

para_names = {'rcp syn','C deg','end','TAK ac', 'time d1','time d2','Tot NFkB'};

for i_feature = 1:length(features)
    corr_mat_codon = NaN(7,5*3);
    i_cond = 1;
    for i_sti = 1:length(ligand_vec)
        ligand = ligand_vec{i_sti};
        % dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
        for i_dose = 1:length(dose_vec{i_sti})
            dose  = dose_vec{i_sti}{i_dose};
            dose_sym = dose_symbol{i_dose};
            para_vals = data.parameters_mode_nan{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose};
            para_vals = para_vals(:,1:end-1);
            codon_corr = collect_feature_vects.(features{i_feature}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'};
            para_corr_core = para_vals(:,[1:3,size(para_vals,2)]);
            para_corr_rcp= para_vals(:,setdiff(1:size(para_vals,2),[1:3,size(para_vals,2)],'stable'));
            para_corr = [para_corr_core,para_corr_rcp];
            
            writematrix(para_corr,strcat(data_save_path,'X_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
            writematrix(codon_corr,strcat(data_save_path,'y_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
            
            
            %% other recptor parameter neg control
            para_neg_control = zeros(size(para_vals,1),0);
            proj_num_vec = setdiff(proj_num_all,proj_num_all(i_sti),'stable');
            for i_proj_num = 1:length(proj_num_vec)
                [para_val,est]  = sample_est_parameters_2023(proj_num_vec(i_proj_num),size(para_vals,1),monolix_data_save_file_path,1);
                NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
                non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
                para_rcp_index = 4:length(non_NFkB_index);
                para_neg_control = [para_neg_control,para_val(:,para_rcp_index)];
            end
            writematrix(para_neg_control,strcat(data_save_path,'X_neg_control_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));

            
            %% uniform distribution neg control
            para_neg_control_uni_dist = rand(size(para_vals,1),50000);
            [corr_para,pval_para] = corr(para_corr,para_neg_control_uni_dist);
            I_pval = find(sum(pval_para>0.6)==size(para_vals,2)& sum(abs(corr_para)<0.05)==size(para_vals,2),size(para_neg_control,2));%
            para_neg_control_uni_dis_save = para_neg_control_uni_dist(:,I_pval);
            
            if length(I_pval)<size(para_neg_control,2)
                i_cond,i_feature
            end
            writematrix(para_neg_control_uni_dis_save,strcat(data_save_path,'X_neg_control_unidistr_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
            
            i_cond=i_cond+1;
        end
    end
end

