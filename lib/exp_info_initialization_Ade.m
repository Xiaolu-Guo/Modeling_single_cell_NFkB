% should be update if the experiments are new


%filepath=filepath_all{stim};

% filepath_all={'../Model_TNF_v6/'
%     '../Model_CpG_0604/'
%     '../Model_polyIC_0612/'
%     '../Model_P3CSK_0604/'
%     '../Model_LPS_0612/'};

% plot_traj=0;

%num_dose=num_dose_all{stim};%5;

num_dose_all={3 %TNF
    5 %CpG
    3 %P3CSK
    5 %LPS
    3}; %polyIC

%num_cells_each_dose ={stim};
num_cells_each_dose_all = {{420, 744, 645} %TNF
    {258, 247, 381, 308, 394} %CpG
    {507,612,596 } %P3CSK
    {176,245,437,310,327 } %LPS
    {666, 540, 557} }; %polyIC


%
% str_dose= str_dose_all{stim};%
dose_str_all= {{'100pg','1ng','10ng'} %TNF
    {'10nM','33nM','100nM','333nM','1uM'} %CpG
    {'10ng','100ng','1ug'} %P3CSK
    {'1ng','3ng','10ng','33ng','100ng'} %LPS
    {'10ug','33ug','100ug'} }; %polyIC

%stimuli=stimuli_all{stim};%'CpG';
% stimuli_all={'TNF'
%     'CpG'
%     'polyIC'
%     'P3CSK'
%     'LPS'};

%est_name=est_name_all{stim};%'Model_CpG_0604';
% est_name_all = {'Model_TNF_v6'
%     'Model_CpG_0604'
%     'Model_polyIC_0612'
%     'Model_P3CSK_0604'
%     'Model_LPS_0612'};

ligand_all = { 'TNF';
    'CpG';
    'Pam3CSK';
    'LPS'
    'polyIC' };


dose_all = {{'x0_1','x1','x10'} %TNF
    {'x10','x33','x100','x333','x1000'} %CpG
    {'x10','x100','x1000'} %Pam3CSK4
    {'x1','x3_3','x10','x33','x100'} %LPS
    {'x10','x33','x100'} }; %polyIC

dose_scale_all = {5200 %TNF
    1000 %CpG
    1500 % P3CSK
    24000 %LPS
    5e6 }; %polyIC

dose_val_all = {{0.1,1,10} %TNF
    {10,33,100,333,1000} %CpG
    {10,100,1000} % P3CSK
    {1,3,10,33,100} %LPS
    {1000*10,1000*33,1000*100}}; %polyIC


% save('initialzation_info')