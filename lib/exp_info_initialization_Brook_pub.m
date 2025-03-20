% should be update if the experiments are new


%filepath=filepath_all{stim};

% filepath_all={'../Model_TNF_v6/'
%     '../Model_CpG_0604/'
%     '../Model_polyIC_0612/'
%     '../Model_P3CSK_0604/'
%     '../Model_LPS_0612/'};

% plot_traj=0;

%num_dose=num_dose_all{stim};%5;

num_dose_all={5 %TNF
    5 %CpG
    7 %LPS
    5 %P3CSK
    4}; %polyIC

%num_cells_each_dose ={stim};
num_cells_each_dose_all = {{488, 416, 433, 423, 346} %TNF
    {378,347,392,423,296 } %CpG
    {458,301,449,617,2019,546,557 } %LPS
    {287, 375, 354, 467, 448}  %P3CSK
    {360, 370, 423, 359} }; %polyIC


%
% str_dose= str_dose_all{stim};%
% dose_str_all= {{'100pg','1ng','10ng'} %TNF
%     {'10nM','33nM','100nM','333nM','1uM'} %CpG
%     {'1ng','3ng','10ng','33ng','100ng'} %LPS
%     {'10ng','100ng','1ug'} %P3CSK
%     {'10ug','33ug','100ug'} }; %polyIC

dose_str_all= {{'330pg','1ng','3.3ng','10ng','33ng'} %TNF
    {'10nM','33nM','100nM','333nM','1uM'} %CpG
    {'330pg','1ng','3.3ng','10ng','33ng','100ng','330ng'} %LPS
    {'1ng','3.3ng','10ng','33ng','100ng'} %P3CSK
    {'3.3ug','10ug','33ug','100ug'} }; %polyIC
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
    'LPS'
    'Pam3CSK';
    'polyIC' };


dose_all = {{'x0_33','x1','x3_3','x10','x33'} %TNF
    {'x1','x3_3','x10','x33','x100'} %CpG
    {'x0_33','x1','x3_3','x10','x33','x100','x330'} %LPS
    {'x1','x3_3','x10','x33','x100'} %Pam3CSK4
    {'x3_3','x10','x33','x100'} }; %polyIC

dose_scale_all = {5200 %TNF
    1000 %CpG
    24000 %LPS
    1500 % P3CSK
    5e6 }; %polyIC

dose_val_all = {{0.33,1,3.3,10,33} %TNF
    {1,3.3,10,33,100} %CpG
    {0.33,1,3.3,10,33,100,330} %LPS
    {1,3.3,10,33,100} % P3CSK
    {1000*3.3,1000*10,1000*33,1000*100}}; %polyIC


% save('initialzation_info')