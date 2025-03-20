% run_me_2023.m for NFkB parameter estimation
% addpath with subfolders: MACKtrack

%All_ligand_codon_2023_t33_cv_filtered_TNF.mat
% data.parameters_mode_nan
% data.exp_mode_filter_nan
% data.pred_mode_filter_nan
% are the TNF cv filtered trajectories

%%
clear all
vers = '202301';

run_sim_example = 0;
rescale_exp_data = 0;
rescale_supriya_exp_data_benchmark = 0;
cal_codon_diff_supriya_exp_data_benchmark = 0;
run_save_all_data_codon = 0;

draw_fig_traj_heatmap = 0;
draw_fig_para_distr = 0;
run_save_codon = 0;
draw_codon = 0;
draw_spearman_corr = 0;
draw_heatmap_para_scatter = 0;
save_para_codon_corr_python = 0;
visualize_para_codon_importance_python = 0;

transfer_codon_data_mutual_info_cal_format = 0;
plot_codon_cc_LPS_TNF = 0;
transfer_codon_data_machine_learning_format = 0;
visualize_stim_classification = 0;

run_sampling = 0;
run_sampling_fitting = 0;
draw_sampling = 0;
draw_sampling_codon = 0;

%_supriya_data_sampling
draw_co_sti_sampling = 0;
draw_co_sti_CpG_polyIC_competetation = 0;

run_multi_sti = 0;
run_5_sti = 0;
run_scramble_dual = 0;
run_5_single_ligand = 0;

% intrinsic
draw_intrinsic_noise = 0;

%% initializing
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
monolix_data_save_file_path = '../raw_data2023/SAEM_proj_2023/';
fig_save_path = '../SubFigures2023/';
data_example_format_path = './example_data_format/';


addpath('./lib/')
addpath('./src/')
addpath('./bin/')
addpath(genpath('./refs_codes/'));


if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

%% [Figure S1] [simulation] & [visualization] parameter sens analysis
% Figure S1, 0507 version: tested 05/09/2024
if 0
    ParameterScan_Sens_traj_04172024()
end

%% [Related to Figure 2] [data preparation] rescaling the data to SI for monolix software input
if 0 % rescale_exp_data and save SAEM format
    % tested 05/15/2024
    data_main
    % cd(pwd)
end

%% [Related to Figure 2] [calculation] calculate signaling codon
% data prepare for all codon visulization
% filter out 2/3 of low variation (CV) of TNF
% 05/10/2024
run_save_codon = 0;
if run_save_codon
    run_all_ligand_codon_2023_03(monolix_data_save_file_path,data_save_file_path)
end

%% [Figure 2] [visualization] Heatmaps of traj, W-dist of signaling codon distri.
% Figure 2, 0507 version: tested 05/10/2024
if 0

    % Figure 2A
    draw_traj_heatmap_2023_03(monolix_data_save_file_path,fig_save_path)

    % Figure 2B-C
    draw_all_ligand_codon_distrib_202404(data_save_file_path,fig_save_path)

end

%% [Related to Figure 3] [data preparation] sampling codon prep for MI calculation
% run till here
if 0% transfer_codon_data_mutual_info_cal_format % transfer the codon data to the mutual information calculation format

    if 1 % tested: 05/15/2024. tansfer sampling data
        data_save_file_path_1 = '../raw_data2023/';
        MI_file_save_path = '../raw_data2023/MI_single_ligand/';
        transfer_data_mutual_info_cal_format(data_save_file_path_1,MI_file_save_path)
    end
    if 1 % tested: 08/13/2024 transfer experiment data
        transfer_exp_data_mutual_info_cal_format_202408
    end
end

if 0 % tested: 08/13/2024 extended doses study: transfer 20doses
    % 'Sim3_codon_r5_metric.mat'; Pam3CSK
    % 'Sim3_codon_r4_metric.mat'; PolyIC
    % 'Sim3_codon_r3_metric.mat'; CpG
    % 'Sim3_codon_r2_metric.mat'; LPS
    % 'Sim3_codon_r1_metric.mat'; TNF

    data_filename = 'Sim3_codon_r5_metric.mat';% please change the name here to run different ligand, the file names are listed above
    save_sampling_20doses_for_channel_capacity_202408(data_filename)
end

%% [related to Figure 2] [data preparation] transfer signaling codons for exp and sim into machine learning format
if 0 % tested: 08/13/2024 transfer_codon_data_machine_learning_format, for Ade's exp, fitting, and sampling dataset
    transfer_data_machine_learning_cal_format_202408
end

%% [Figure 2] [visualization] stimulation classification
% Figure 2D & Figure S2F-G & Figure S3D : 0507 version: : tested 05/30/2024
if 0
    visualize_stim_classification_machine_learning_results
    ML_visual_sens_cross_data
end

%% [Figure S2] [visualization] metrics of good fitting, RMSD, signaling codon distri, and W-dist
% Figure S2: 0507 version: : tested 05/09/2024
if 0

    % Figure S2A RMSD distribution tested 06/18/2024
    draw_traj_RMSD_distrib_2024_06(fig_save_path)

    % Figure S2A
    draw_traj_RMSD_2023_05(fig_save_path)

    % Figure S2B
    draw_traj_RMSD_2024_05(fig_save_path)

    % Figure S2C-E: tested 05/10/2024
    draw_all_ligand_codon_distrib_202309(data_save_file_path,fig_save_path)

end

%% [Related to Figure S3] [simulation] sampling single ligand
if 0
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    % run all conditions

    Initializing_sampling_supriya_data_2023_05
    % single: 1,2,4,6,7
    % dual: 3,5,9,10,14,15,20,26,27
    % tri 21,22

    dual_ligand_length = 10;
    single_ligand_length = 5;

    % input_paras.Num_sample = 1000;

    % the last 5 indices are corresponding to single-ligand stimuli
    % simulation
    index_sim = [3;5;9;10;14;15;20;26;27;29;1;2;4;6;7];%1;2;6;7;11;
    input_paras.proj_num_vec = input_paras.proj_num_vec(index_sim);
    input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_sim);
    input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_sim);
    input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_sim);
    input_paras.Num_sample = [1*ones(dual_ligand_length,1);
        2*ones(single_ligand_length,1)]*300;%
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;

    save_metric_name = 'Sim5_codon_all5dose_metric.mat';
    cal_codon =1;
    % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras)
    para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')

    vers_fig = '2023';
    single_or_dual = 'dual';
end

%% [Figure S3] [visualization] single-ligand sampling results
% Figure S3: 0507 version: : tested 05/10/2024
if 0 % heatmap

    % Figure S3A: 0507 version: : tested 05/10/2024
    draw_traj_heatmap_diff_order_2023_08

    % Figure S3B-C: 0507 version: : tested 05/12/2024
    draw_all_ligand_sampling_codon_distrib_202306(data_save_file_path,fig_save_path)

    if 1 % Figure S3E: plot_codon LPS & TNF : tested 08/13/2024
        scatter_codon_exp_sampling_sim_cc_202408
    end

end

%% [Related to Figure 3] [simulation] doses study
if 0
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    %    run_sampling_simulation(monolix_data_save_file_path,fig_save_path)

    dose_ratio = 10.^[-2.5:0.25:2.5];
    dose_length = length(dose_ratio);
    input_paras.proj_num_vec = [mat2cell(2*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(3*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(4*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(5*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(6*ones(dose_length,1),ones(dose_length,1),[1])];
    input_paras.Num_sample = [1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1)]*600;%600;

    input_paras.proj_ligand_vec = {{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};
        {'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};
        {'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};
        {'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};
        {'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};
        {'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};
        {'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};
        {'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};
        {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};
        {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};};
    input_paras.proj_dose_str_vec = {{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'}};
    a = [mat2cell((1*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((1*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((100*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((10000*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((100*dose_ratio)',ones(dose_length,1),[1])];

    input_paras.proj_dose_val_vec = cell(size(a));
    for i_proj = 1:length(input_paras.proj_dose_val_vec)
        input_paras.proj_dose_val_vec{i_proj} = a(i_proj);
    end

    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;

    cal_codon =1;
    input_paras_filed_names = fieldnames(input_paras);
    for i_ligand = 1:5
        index_ligand = (i_ligand-1)* 21 + (1:21);
        save_metric_name = strcat('Sim3_codon_r',num2str(i_ligand),'_metric.mat');
        for i_field = 1:length(input_paras_filed_names)
            if length(input_paras.(input_paras_filed_names{i_field}))>length(index_ligand)
                input_paras_ligand.(input_paras_filed_names{i_field}) = input_paras.(input_paras_filed_names{i_field})(index_ligand);
            else
                input_paras_ligand.(input_paras_filed_names{i_field}) = input_paras.(input_paras_filed_names{i_field});
            end
        end
        para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,input_paras_ligand,cal_codon,save_metric_name,'alldose')
    end

end

%% [Figure 3] [visualization] for extended doses study of each ligand
% Figure 3A: 0507 version: : tested 05/10/2024
if 0
    ligand_fig_save_path = strcat(fig_save_path,'doses_20/');
    for i_r = 1:5

        data_filename = strcat('Sim3_codon_r',num2str(i_r),'_metric.mat');
        draw_sampling_traj_heatmap_2023_05(ligand_fig_save_path,data_filename)
    end
end

%% [Figure S3] [visualization] extended doses signaling codon distribution
if 0 % Figure S3F : [tested] 06/12/2024 draw 20 doses signaling codon distributions
    % 'Sim3_codon_r5_metric.mat'; Pam3CSK
    % 'Sim3_codon_r4_metric.mat'; PolyIC
    % 'Sim3_codon_r3_metric.mat'; CpG
    % 'Sim3_codon_r2_metric.mat'; LPS
    % 'Sim3_codon_r1_metric.mat'; TNF

    draw_sampling_20doses_codons_distrib()

end

%% [Related to Figure 3] [simulation] & [data preparation] MI loss within the network
% For Figure 4 & S4
% first run simulation to denoise different modules,
% transfer data into MI calculation format
% tested 05/15/2024
if 0

    data_save_file_path = '../raw_data2023/simulation_denoise/';
    MI_file_save_path = '../raw_data2023/MI_denoise/';

    if 1 % [tested] no noise
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params52n2','params65n2',...
            'params99','params101'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            1,1,...
            0.4,0.4]';

        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim24_no_heterogeneity_codon_metric',num2str(i_r),'.mat');

        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);

    end

    if 1 %[tested] only TAK noise Sim22
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params99','params101'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            0.4,0.4]';

        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim22_TAK1_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);

    end

    if 1 % [tested] only core module noise
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params52n2','params65n2'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            1,1]';

        var_input.specienames = {};
        var_input.specievals = [];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim23_NFkB_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);

    end

    if 1 % [tested] only rcp noise
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params52n2','params65n2','params99','params101'};%LPS
        var_input.paravals = [1,1,0.4,0.4]';

        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];

        cal_codon =1;

        for i_r =1:2
            save_metric_name = strcat('Sim21_reduce_TAK1ac_NFkB_heterogeneity_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end


    end

    if 1 % [tested] TAK reduce heterogeneity
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params52n2','params65n2'};%LPS
        var_input.paravals = [1,1]';

        var_input.specienames = {};
        var_input.specievals = [];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim20_reduce_TAK1ac_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);


    end


    if 1 %% [tested] receptor reduce heterogeneity
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35'};%LPS
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05]';

        var_input.specienames = {};
        var_input.specievals = [];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim19_reduce_rcpt_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);


    end

    if 1 %% [tested] receptor reduce heterogeneity for IkBas/s
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample =1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...%LPS
            'params6'};% NFkB regulated IkBa transcriptional rate
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            6e-05*0.25]';

        var_input.specienames = {};
        var_input.specievals = [];


        cal_codon =1;

        i_r =1;
        save_metric_name = strcat('Sim19_IkBas_reduce_rcpt_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);


    end

    if 1 % [tested] control: wt, IkB-/- % only use wt.
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 5;%1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        cal_codon =1;

        para_vec = [6e-05;
            6e-05*0.1;
            6e-05*0.01;
            6e-05*0];

        ncpu=2;
        pc=parcluster('local');
        pc.NumThreads=2;%
        parpool(pc,ncpu)
        var_input = cell(4);

        for i_r =1:size(para_vec,1)
            var_input{i_r}.paranames = {'params6'};
            var_input{i_r}.paravals = para_vec(i_r,:)';
            var_input{i_r}.specienames = {};
            var_input{i_r}.specievals = [];

            save_metric_name = strcat('Sim18_wt_IkBo_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input{i_r}, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end

        delete(gcp)

    end

    if 1 %[tested]  generate data, reducing NFkB var, using representaive, vs 9 random picked selected cells
        % purpose: assuming that all information are transimitted into NFkB

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

        var_input.paranames = {'params99','params101'};
        var_input.specienames = {'NFkB'};

        para_vec = [0.4,0.4;
            0.336004000000000,0.633038000000000;
            0.0673386000000000,0.0619458000000000;
            0.0157939000000000,0.0216838000000000;
            0.217826000000000,0.417812000000000;
            0.424314000000000,0.445427000000000;
            0.0115429000000000,0.0258170000000000;
            0.0802117000000000,0.0685842000000000;
            0.493319000000000,0.669808000000000;
            0.0241970000000000,0.0200540000000000];

        specie_vec = [0.08;
            0.0564181000000000;
            0.242900000000000;
            0.0701397000000000;
            0.0497623000000000;
            0.114479000000000;
            0.0974905000000000;
            0.0421724000000000;
            0.0721400000000000;
            0.0858496000000000];


        cal_codon =1;

        ncpu=2;
        pc=parcluster('local');
        pc.NumThreads=1;%
        parpool(pc,ncpu)
        var_input = cell(10);

        parfor i_r =1:size(para_vec,1)
            var_input{i_r}.paranames = {'params99','params101'};
            var_input{i_r}.specienames = {'NFkB'};
            var_input{i_r}.paravals = para_vec(i_r,:)';
            var_input{i_r}.specievals = specie_vec(i_r,:)';

            save_metric_name = strcat('Sim17_reduce_core_heterogeneity_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input{i_r}, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end

        delete(gcp)

    end

    if 1 % [tested]
        transfer_diff_var_codon_mutual_info_format_0331(data_save_file_path,MI_file_save_path)
    end

end

%% [Figure S4] [visualization] heatmaps and signaling codon distribution for WT and Sjroen syndrome
if 0 % Fgiure S4B-C
    draw_heatmap_codon_SS_p1x_0331
end



%% [Figure 4] [calculation & visualization] coefficient of variance of parameter distribution
if 0% [tested] 06/18 cal para CV

    % purpose: assuming that all information are transimitted into NFkB


    monolix_data_save_file_path = '../SAEM_proj_2023/';
    input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

    data_label = {'TNF-syn','TNF-deg',...
        'LPS-syn','LPS-deg','LPS-endo',...
        'CpG-syn','CpG-deg','CpG-endo',...
        'PolyIC-syn','PolyIC-deg','PolyIC-endo',...
        'Pam3CSK-syn','Pam3CSK-deg','adp','core-timed1','core-timed2','core-NFkBtot'};

    %ade's data highest dose
    input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
    input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
    input_paras.Num_sample = 1000;%5 ligand
    input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

    var_input.paranames = {};%core IKK-NFkB-IKBa
    var_input.paravals = []';

    var_input.specienames = {};
    var_input.specievals = [];

    cal_codon =1;

    i_r =1;
    save_metric_name = strcat('Simxx_sampling_cal_CV','.mat');
    if 0 % new sampling data
        rcp_cvs = [];
        adp_cvs = [];
        core_cvs = [];
        for i=1:10
            % [rcp_cvs_rpc,adp_cvs_rpc,core_cvs_rpc] = para_sampling_cal_CV(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
            [rcp_cvs_rpc,adp_cvs_rpc,core_cvs_rpc] = para_sampling_cal_CV_barplot(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);

            rcp_cvs = [rcp_cvs,rcp_cvs_rpc];
            adp_cvs =[adp_cvs,adp_cvs_rpc];
            core_cvs = [core_cvs,core_cvs_rpc];
        end
        save(strcat('../raw_data2023/2024/',save_metric_name), 'rcp_cvs', 'adp_cvs', 'core_cvs')

    else % load previously saved data

        load(strcat('../raw_data2023/2024/',save_metric_name))
    end


    if 1 % bar plot of paramter cv
        figure(1)
        paperpos=[0,0,300,100]*1.5;
        papersize=[300 100]*1.5;
        draw_pos=[10,10,290,90]*1.5;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        core_cvs_new = [core_cvs(1:3,:),core_cvs(4:6,:),core_cvs(7:9,:),core_cvs(10:12,:),core_cvs(13:15,:)];
        % cvs_all = [rcp_cvs;adp_cvs;core_cvs];
        % rcp_cvs = rcp_cvs(:,[
        bar_means_all = mean(rcp_cvs, 2);
        bar_stddevs_all = std(rcp_cvs, 0, 2);
        bar_means_all(end+1) = mean(adp_cvs(:));
        bar_stddevs_all(end+1) = std(adp_cvs(:));
        bar_means_all(end+1:end+3) = mean(core_cvs_new,2);
        bar_stddevs_all(end+1:end+3) = std(core_cvs_new, 0, 2);



        mean(mean(rcp_cvs, 2))
        mean(mean(adp_cvs(:)))
        mean(mean(core_cvs_new,2))
        %  data_label = {'TNF-syn','TNF-deg',...
        %             'LPS-syn','LPS-deg','LPS-endo',... 3-5
        %             'CpG-syn','CpG-deg','CpG-endo',... 6-8
        %             'PolyIC-syn','PolyIC-deg','PolyIC-endo',... 9-11
        %             'Pam3CSK-syn','Pam3CSK-deg','adp','core-timed1','core-timed2','core-NFkBtot'};% 12, 13, 14 adp, 15-17
        data_label_index =  [1,2,12,13,6,7,8,3,4,5,9,10,11,14,15,16,17];

        bar_means = bar_means_all(data_label_index);
        bar_stddevs = bar_stddevs_all(data_label_index);
        data_label = data_label(data_label_index);

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        c = categorical(data_label);
        c = reordercats(c,data_label);
        b = bar(c,bar_means,0.4, 'FaceColor', 'flat','EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed

        numGroups = length(bar_means);
        numBars = length(bar_means);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        x = (1:numGroups) ; % Adjust the position - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars)

        % Define colors
        colors = [repmat([119, 180, 202]/255, 2, 1);   % for TNF bars
            repmat([229,129,56]/255, 2, 1);   % for Pam bars
            repmat([137,180,66]/255, 3, 1);   % for CpG bars
            repmat([222,78,66]/255, 3, 1); % for LPS bars
            repmat([101,77,123]/255, 3, 1); % for pIC bars
            repmat([0.7 0.7 0.7], 4, 1);]; % for core bars


        % Assign colors to bars
        b.CData = colors;
        % color:


        % Add error bars
        errorbar(x, bar_means, bar_stddevs, 'k', 'linestyle', 'none');

        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        %             xticklabels({});
        %             yticklabels({});
        % title(rmsd_cas)
        set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')

        x_positions = 1:numGroups; % X positions of bars

        for i = 1:numGroups
            % Extract raw data points
            clear raw_data
            if i< 14
                raw_data = rcp_cvs(data_label_index(i), :); % Use corresponding data from info_mat
            elseif i < 15
                raw_data = adp_cvs(:);
            else

                raw_data = core_cvs_new(i-14,:);
            end

            % jitter = linspace(-0.3,0.3,length(raw_data));%(rand(size(raw_data)) - 0.5) * 1; % Add jitter for better visibility
            jitter = (rand(size(raw_data)) - 0.5) * 0.2; % Add jitter for better visibility
            scatter(x_positions(i) + jitter, raw_data, 3, 'k','filled'); % Plot raw data;
        end

        hold off;

        ylim([0,2.5])
        saveas(gcf,strcat(fig_save_path,'cvs_bar_rcp_adp_core'),'epsc');

        close

    end

end

%% [related to Figure 6] [simulation] single cell specificity
if 0

    %% [tested] [simulation] 5 single ligand stim: matching
    if 1  % tested 05/15/2024, Matlab 2020a
        %run_5_single_ligand: to get 'Sim8_5_signle_ligand_codon_metric_r3.mat'
        % doses info:
        % ligands: {{'TNF'};
        %     {'LPS'};
        %     {'CpG'};
        %     {'PolyIC'};
        %     {'Pam3CSK'}};

        % dose_str: {{'10ng/mL'};
        %     {'10ng/mL'};
        %     {'100nM'};
        %     { '100ug/mL'};
        %     {'100ng/mL'}};

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample = 200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        %save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r2.mat';
        for i_r = 3%:7
            save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_sim_single_ligand_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
    end

    %% [tested] [simulation] single cell specificity for IkBa-/- match
    if 1 % tested 05/15/2024, Matlab 2020a

        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample = 200;%5 ligand 200
        input_paras.proj_num_vec ={[2,3,4,5,6]};

        cal_codon =1;

        if 1 % 03/22/2024 params6 = 0.25 * params6_wt;
            for  i_r =1
                save_metric_name = strcat('Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x_r',num2str(i_r),'.mat');
                % params6 = 0.1 * params6_wt;
                para_sample_fitting_match_sim_SRS_IkBao_202401(0.25, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
            end
        end

    end

end

%% [Figure 6 & S6] [visualization] draw single cell specificity
if 0 % Fgiure 6

    % figure 6
    if 1
        draw_single_cell_confusion_20240404
    end
    % figure S6
    if 1
        Figure_S7A_single_cell_specifcity
    end
end



%% [Figures in supplementary notes] benchmarking different methods combine receptor modules
if 0 % testing different matching method
    %% [tested] [supplementary materials] 5 single ligand stim: weighted matching
    if 0 % tested 05/15/2024, Matlab 2020a
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        save_metric_name = 'Sim15_5_signle_ligand_codon_metric.mat';

        % ncpu = 5;
        % pc = parcluster('local');
        % pc.NumThreads = 2;
        % parpool(pc,ncpu)
        % par
        for i_r = 2 :6

            save_metric_name = strcat('Sim15_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_weight_match_sim_single_ligand_var_2023_12(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        % delete(gcp)
    end

    %% [tested] [supplementary materials] 5 single ligand stim: scramble
    if 0 % tested 05/15/2024, Matlab 2020a
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        %    save_metric_name = 'Sim14_5_signle_ligand_codon_metric_r2.mat';

        % ncpu=5;
        % pc=parcluster('local');
        % pc.NumThreads=2;%
        % parpool(pc,ncpu)
        %
        % par
        for  i_r =3:7
            save_metric_name = strcat('Sim14_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_sim_single_ligand_scramble_2024(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        % delete(gcp)

    end

    %% [tested] [supplementary materials] single cell compare draw replicates
    if 0 % please make sure the files Sim15_5_signle_ligand_codon_metric_r3.mat etc is run and saved
        %draw_5_single_ligand_sampling_exp_codon_v202401
        draw_5_single_ligand_sampling_exp_codon_v20240120
    end

end

if 0 % using machine learning to learn parameter signaling codon relationships
    %%  [tested] [supplementary materials] save the parameter signaling codon for python codes: linear regression, RandomForest, XGBoost
    % For Figure supp:
    % save the data of parameters (X), and singnaling codon (y), and random
    % variables as negative control. these dataset will be used to train the
    % regression models in python. the outputs of the python codes have been
    % saved in raw_data2023/singaling_codon_para/, and are visualized in the
    % next two sections
    if 0
        py_data_save_path = strcat(data_save_file_path,'singaling_codon_para_for_python/');

        save_csv_para_codon_2023_03(data_save_file_path,py_data_save_path,monolix_data_save_file_path)
    end

    %% [tested] [supplementary materials] draw R-squared for regression models  and spearman corr and feature importance between para and codon
    if 0 % Figure supp
        % Figure A C: 0507 version: : tested 05/12/2024
        visualize_para_codon_feature_importance_2024_04(data_save_file_path,py_data_save_path,fig_save_path)

        % Fgiure B: 0507 version: : tested 05/12/2024
        draw_spearman_corr_each_codon_2024_04(data_save_file_path,fig_save_path)

    end

    %% [tested] [supplementary materials] draw codon vs parameter: heatmap, spearman corrlation, feature importance
    if 0 % Figure supp
        py_data_save_path = strcat(data_save_file_path,'singaling_codon_para/');

        % Figure A: 0507 version: : tested 05/12/2024
        para_cells_scatter_2023_06();

        % Figure B spearman correlation plot: 0507 version: : tested 05/12/2024
        draw_spearman_corr_each_codon_2023_06(data_save_file_path,fig_save_path)

        % Figure C: 0507 version: : tested 05/12/2024
        visualize_para_codon_feature_importance_2023_03(data_save_file_path,py_data_save_path,fig_save_path)

    end
end


%% Initializing for [Figure 5], [Figure S2H-L], [Figure S3E-G], [Figure 7]
clear all
vers = '202301';

data_save_file_path = '../raw_data2023/2024/';%_fay_parameter/';
monolix_data_save_file_path = '../raw_data2023/SAEM_proj_2023/';
fig_save_path = '../SubFigures2023/revision_2024/';
data_example_format_path = './example_data_format/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')
addpath('./revision/')

if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

%% [Related to Figure 5] [Simulation] IkBamut compare with Ade2021 dataset 
if 0
    % generating new results for IkBa mutant & WT
    % IkBaMut_LPS_100ng vs IkBaWt
    % IkBaMut_PolyIC_50ug vs IkBaWt
    % IkBaMut_TNF10ng vs IkBaWt
    % [tested] control: wt, IkB-/- % only use wt.
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};

    %ade's data highest dose
    input_paras.proj_dose_str_vec = {{'10ng/mL'},{'100ng/mL'},{'333nM'},{'50ug/mL'},{'1000ng/mL'}};
    input_paras.proj_dose_val_vec = {{10},{100},{333},{50000},{1000}};
    input_paras.Num_sample = 1000;%1000;%5 ligand
    input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};

    cal_codon =1;

    para_vec = [6e-05;
        6e-05*0.25];
    vers = '_revision_1115' ;

    ncpu=2;
    pc=parcluster('local');
    pc.NumThreads=2;%
    parpool(pc,ncpu)
    var_input = cell(4);
    name_vec = {'_wt','_ikbamut'}

    for i_r = 1:size(para_vec,1)
        var_input{i_r}.paranames = {'params6'};
        var_input{i_r}.paravals = para_vec(i_r,:)';
        var_input{i_r}.specienames = {};
        var_input{i_r}.specievals = [];

        save_metric_name = strcat('Sim18_wt_IkBo_codon_metric',name_vec{i_r},vers,'.mat');
        para_fitting_sim_alter_var_202402(var_input{i_r}, data_save_file_path,input_paras,cal_codon,save_metric_name);
    end

    delete(gcp)

end

%% [Figure 5 & S5] ML classification: simulation and ade2021 exp IkBamut
if 0
    if 0 % data preparation, for ML classifier run in python
        transfer_data_ML_mutant_sim_20241115 % simulation data
        transfer_data_ML_ade2021_exp_20241101 % exp data (Adelaja 2021)
    end

    if 0 % Figure S5: plot the signaling codon distribution
        visualize_signaling_codon_sim_for_ade2021_exp % simulation data
        visualize_signaling_codon_ade2021_exp % exp data (Adelaja 2021)
    end

    if 0 % Figure S5: plot the trajectories
        visualize_traj_sim_for_ade2021_exp % simulation data
        visualize_ade2021_data_20241218 % visualize Ade2021 exp IkBamut and wt dataset
    end

    if 0 % plot the machine learning classifier results
        ML_visulize_Ade2021_IkBamut_20241218
    end
end



%% [Figure S3] quantitative assessing extend dose sim and experiments
if 0
    % get the dataset:
    Quan_assess_extend_dose_1207
end

%% [Figure S3] run extend dose sim to compare with Adelaja 2021
if 0
    % generating new results for IkBa mutant & WT
    % IkBaMut_LPS_100ng vs IkBaWt
    % IkBaMut_PolyIC_50ug vs IkBaWt
    % IkBaMut_TNF10ng vs IkBaWt
    % [tested] control: wt, IkB-/- % only use wt.
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'PolyIC'}};

    %ade's data highest dose
    input_paras.proj_dose_str_vec = {{'33ng/mL'},{'333ng/mL'},{'50ug/mL'}};
    input_paras.proj_dose_val_vec = {{33},{333},{50000}};
    input_paras.Num_sample = 600;%1000;%5 ligand
    input_paras.proj_num_vec ={[2],[3],[5]};

    cal_codon =1;

    para_vec = [6e-05];
    vers = '_revision_1130' ;

    var_input = cell(1);
    name_vec = {'_ext_dose'};


    var_input{1}.paranames = {'params6'};
    var_input{1}.paravals = [para_vec];
    var_input{1}.specienames = {};
    var_input{1}.specievals = [];

    save_metric_name = strcat('Sim',name_vec{1},'_codon_metric_202402',vers,'.mat');
    para_fitting_sim_alter_var_202402(var_input{1}, data_save_file_path,input_paras,cal_codon,save_metric_name);
    save_metric_name = strcat('Sim',name_vec{1},'_codon_metric_2023_05',vers,'.mat');
    para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')

end
%% [Figure S2] visualize SAEM Pam dataset ( Pam  downsampling)
if 0
    visualize_Pam_downsampling_data
end


%% [Figure 7 & S7] [visualization] draw single cell confusion
% revision for Figure 7 & S7
% to check: IkBa mutant [done]
if 0
    scSRS_corr_heterogeneity_202410
end
