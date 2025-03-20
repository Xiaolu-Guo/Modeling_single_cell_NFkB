function [] = Module_parameters_sim_04172024(data_save_file_path,module_info)

vers = '_04172024';

sim_data_tbl = SimDataTblInitialize();

features_metric_fields = {'duration','oscpower','max_pos_integral','pos_pk1_time','max_value','time2HalfMaxPosIntegral','time2HalfMaxValue'};%'max_pos_pk1_speed',,'pos_pk1_amp','time2HalfMaxPosIntegral'

for i_fields = 1:length(features_metric_fields)
    sim_data_tbl.(features_metric_fields{i_fields}) = zeros(0);
end

data_info.save_file_path = data_save_file_path;%  '../../NFkB_para_estm_project/NFkB_figures/ParameterScan/data/';
data_info.species_outputname = {'nucNFkB'};%;'TNFR';'IKK'
data_info.species_composition = {{'NFkBn';'IkBaNFkBn'}}; % must ;{'TNFR'};{'IKK'}
% be r x 1, for each cell i must be ri x 1
data_info.flag = '';
data_info.type = 'wt';

if isfield(module_info,'module')
    data_info.save_file_name = strcat('Module_Sens_',module_info.module{1},'_',module_info.ligand{1},vers);
else
    data_info.save_file_name = strcat('Module_Sens_',module_info.ligand{1},'_',module_info.para_type{1},vers);
end



for i_ligand = 1:length(module_info.ligand) 
    clear sim_info
    sim_info.ligand = module_info.ligand(i_ligand);
    sim_info.dose_str = module_info.dose_str(i_ligand);
    sim_info.dose_val = module_info.dose_val(i_ligand);
    
    
    %% parameter setting
    
    switch module_info.para_type{i_ligand}
        case 'deg'
            switch module_info.ligand{i_ligand}
                case 'TNF'
                    sim_info.parameter_name = {'params64','params61','params58'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21); 10.^linspace(-1,1,21);10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.125,0.125,0.125];                    
                case 'LPS'
                    sim_info.parameter_name = {'params44'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.012];                   
                case 'Pam3CSK'
                    sim_info.parameter_name = {'params75'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21); 10.^linspace(-1,1,21);10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.004];                    
                case 'CpG'                    
                    sim_info.parameter_name = {'params93'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [1.60E-03];                    
                case 'polyIC'                    
                    sim_info.parameter_name = {'params83'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [ 0.0007];
                otherwise                   
                    error('wrong type of ligand!')                    
            end
            
        case 'syn'
            switch module_info.ligand{i_ligand}
                case 'TNF'
                    sim_info.parameter_name = {'params54'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [8.224e-06];
                case 'LPS'
                    sim_info.parameter_name = {'params35'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [5.25E-05];
                case 'CpG'                    
                    sim_info.parameter_name = {'params85'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [2e-6];                    
                case 'polyIC'
                    sim_info.parameter_name = {'params77'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [3e-6];
                case 'Pam3CSK'                    
                    sim_info.parameter_name = {'params68'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [1e-6];
                otherwise                    
                    error('wrong type of ligand!')                    
            end
            
        case 'endo'
            switch module_info.ligand{i_ligand}                
                case 'LPS'                    
                    sim_info.parameter_name = {'params40','params36'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21); 10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.065681,0.065681];                    
                case 'CpG'
                    sim_info.parameter_name = {'params88'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.015];                    
                case 'polyIC'
                    sim_info.parameter_name = {'params79'};
                    sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
                    sim_info.parameter_value_vec = [0.04];                    

                otherwise                    
                    error('wrong type of ligand!')
            end
            
        case 'TAKac'
            sim_info.parameter_name = {'params52n2','params65n2'};
            sim_info.fold_change_vec = [10.^linspace(-1,1,21); 10.^linspace(-1,1,21)];
            sim_info.parameter_value_vec = [1,1];  
            
        case 'timed1'
            sim_info.parameter_name = {'params99'};
            sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
            sim_info.parameter_value_vec = [0.4]; 
            
        case 'timed2'
            sim_info.parameter_name = {'params101'};
            sim_info.fold_change_vec = [10.^linspace(-1,1,21)];
            sim_info.parameter_value_vec = [0.4];    
            
        case 'NFkBtot'
            sim_info.species_vec.name = {'NFkB'};
            sim_info.species_vec.val = [10.^linspace(-1,1,21) * 0.08];           
            sim_info.parameter_name = {'params99'};
            sim_info.fold_change_vec = [ones(1,21)];
            sim_info.parameter_value_vec = [0.4];  
         
        otherwise
            error('wrong para_type!')
    end
        


sim_info.parameter_value_vec = diag(sim_info.parameter_value_vec) * sim_info.fold_change_vec;


%%

sim_data_tbl_tmpt = NFkB_signaling_para_conc_sim_2023(sim_info,data_info);

% the parameter numbers has to be the same

% NFkBn + IkBaNFkBn
clear data
data.model_sim{1} = sim_data_tbl_tmpt.trajectory(:,1:5:end);
data.info_ligand{1} = sim_data_tbl_tmpt.ligand(1);
data.info_dose_index{1} = i_ligand;
data.info_dose_str{1} = sim_data_tbl_tmpt.dose_str(1);
data.info_num_cells{1} = size(data.model_sim{1},1);
data.order{1} = (1:data.info_num_cells{1})';
data.exp = data.model_sim;
vis_data_field = {'model_sim'};%,'sample'};
data_label = {'simulation'};%,'sample'};

[~,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter

for i_fields = 1:length(features_metric_fields)
    sim_data_tbl_tmpt.(features_metric_fields{i_fields}) = metrics{1}.(features_metric_fields{i_fields});
end

%clear time2HalfMaxValue_tmp
sim_data_tbl_tmpt.time2HalfMaxValue = get_time2HalfMaxValue(data.model_sim{1});
% = time2HalfMaxValue_tmp;

sim_data_tbl= [sim_data_tbl;sim_data_tbl_tmpt];

end



if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
    save(strcat(data_info.save_file_path,data_info.save_file_name),'sim_data_tbl');
end

end


function time2HalfMaxValue = get_time2HalfMaxValue(time_series_no_base_ded)

% adapted from  get_fold_change_new_20220701 for simulation sensitivity
% analysis
%calculates max fold change in time_series matrix (non-baseline deducted)

baseline = time_series_no_base_ded(:,1)*ones(1,size(time_series_no_base_ded,2));
FramesPerHour = 12;

[output.max_value, idx_max] = nanmax(time_series_no_base_ded - baseline, [], 2);

index_mat = ones(size(time_series_no_base_ded,1),1) * 1:size(time_series_no_base_ded,2) <= idx_max;

halfMaxIntegral = output.max_value/2;

distances = abs(time_series_no_base_ded - baseline - halfMaxIntegral);
distances = distances .* index_mat + 100000* ~index_mat;

[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
time2HalfMaxValue = (idx-1)/FramesPerHour;

end