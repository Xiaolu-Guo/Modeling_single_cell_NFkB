% TNF SIMULATION
%p_mod = [idx, all_x(:,idx(1))];
% run('../NFkB_common/exp_info_initialization.m')
addpath('./bin/');
addpath('./lib/');
addpath('./src/');


stimuli_info_tbl = get_stimuli_info_tbl();

names = {'C1','IKK','NFkBn'};% 'TAK',
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

sti = 'TNF';
doses = [10];%0.1,1,
dose_str = {'10ng/mL'};%'100pg/mL','1ng/mL',

% sti = 'Pam3CSK';
% doses = [10,100,1000];
% dose_str = {'10ng','100ng','1ug'};

%
% sti = 'polyIC';
% doses = [1000*10,1000*33,1000*100];
% dose_str = {'10ug','33ug','100ug'};

% sti = 'CpG';
% doses = [10,33,100,333,1000];
% dose_str = {'10nM','33nM','100nM','333nM','1uM'};
%
%
% sti = 'LPS';
% doses = [1,3,10,33,100];
% dose_str = {'1ng','3ng','10ng','33ng','100ng'};


% doses = cell2mat(dose_val_all{strcmp(ligand_all,sti)});
%             doses = stimuli_info_tbl.dose_val((stimuli_info_tbl.Ligand == sti));
% dose_scale = 1/dose_scale_all{strcmp(ligand_all,sti)}; % Convert to uM (from nM)
switch sti
    case 'polyIC'
        sti_dose_scale = 'PolyIC';
    otherwise
        sti_dose_scale = sti;
        
end
dose_scale = get_stimuli_info_dose_scale(sti_dose_scale);%             dose_field = dose_all{strcmp(ligand_all,sti)};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added

options.params_size = [99,1;100,1;101,1];
options.params_value = [0.4,0.4,0.4];
[v0.PARAMS, v0.SPECIES] = nfkbInitialize();%_TNF_TRAF2_paraset2
para_vals = [5.03779756734822;0.00999870619270433;4.99747460828447;5001.04796566778;1799.97518484732];
v0.PARAMS(59,1) = para_vals(1);
v0.PARAMS(60,1) = para_vals(2);
v0.PARAMS(62,1) = para_vals(3);
v0.PARAMS(63,1) = para_vals(4);
v0.PARAMS(65,1) = para_vals(5);

v0.PARAMS(52,2) = 0.636376;
v0.PARAMS(99,1) = 0.362951;
v0.PARAMS(101,1) = 0.532915;
v0.PARAMS(54,1) = 8.22649e-7;
v0.PARAMS(55,1) = 2.42833e-15;
v0.PARAMS(56,1) = 20000;
v0.PARAMS(57,1) = 1;
v0.PARAMS(58,1) = 0.000125009;
v0.PARAMS(59,1) = 1.58318;
v0.PARAMS(60,1) = 2.01374e-15;
v0.PARAMS(62,1) = 49.9958; 
v0.PARAMS(63,1) = 5.7768;
v0.PARAMS(65,1) = 0.0219448;
v0.PARAMS(65,2) = 0.636376;

v0.PARAMS(52,2)= 0.331;
v0.PARAMS(65,2)= 0.331;
v0.PARAMS(99,1)=   0.391;
v0.PARAMS(101,1)=  0.278;
v0.PARAMS(54,1)=   2.94e-6;
v0.PARAMS(55,1)=   0.00153;
v0.PARAMS(56,1)=   2e4;
v0.PARAMS(57,1)=  1.21e-8;
v0.PARAMS(58,1)=  0.494;
v0.PARAMS(59,1)=  217;
v0.PARAMS(60,1)=  0.124;
v0.PARAMS(62,1)=  2.01;
v0.PARAMS(63,1)=  1.12e3;
v0.PARAMS(65,1)=  1.27e3;

% cell 1698 fit peak and trough
v0.PARAMS(52,2) =0.190971605302866;
v0.PARAMS(99,1) =0.999981501324547;
v0.PARAMS(101,1) =0.424308530705832;
v0.PARAMS(54,1) =2.56142543633536e-6;
v0.PARAMS(55,1) =0.00407136814760358;
v0.PARAMS(56,1) =19999.4074685683;
v0.PARAMS(57,1) =0.00226217100589618;
v0.PARAMS(58,1) =0.137815834106994;
v0.PARAMS(59,1) =29.0689153836765;
v0.PARAMS(60,1) =0.0398500238297566;
v0.PARAMS(62,1) =2.76960730927809;
v0.PARAMS(63,1) =452.155762244646;
v0.PARAMS(65,1) =987.280867696068;
NFkB_cyto_init =0.0400000700857019;

v0.PARAMS(65,2) = v0.PARAMS(52,2);

v0.PARAMS(61,1) = 0;
v0.PARAMS(64,1) = 0;


options.v.PARAMS = v0.PARAMS;
options.v.SPECIES = v0.SPECIES;
options.v.init_val.NFkB = 0.2;
options.v.init_val.NFkB = NFkB_cyto_init;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate all doses (only need to equilibrate on first iteration)
output = [];

for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options);% _TNF_TRAF2
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [t,x] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options); % _TNF_TRAF2
    end
    output = cat(3,output,x);
end

C1_curves = squeeze(output(:,strcmp(names,'C1'),:));
ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

for ii = 1:length(doses)
    figure(1)
    plot(1:size(nfkb_curves,1),nfkb_curves(:,ii),'LineWidth',2);hold on
    
    figure(2)
    plot(1:size(C1_curves,1),C1_curves(:,ii),'LineWidth',2);hold on
    
    figure(3)
    plot(1:size(ikk_curves,1),ikk_curves(:,ii),'LineWidth',2);hold on 
end

for ii =1:length(doses)
    legend_str{ii} = num2str(doses(ii));
end



figure(1)
legend(dose_str)
title(sti)
ylabel('NFkBn')
xlabel('Time')
set(gca,'FontSize',14,'FontWeight','b')
xlim([0,480])

figure(2)
legend(dose_str)
title(sti)
ylabel('C1')
xlabel('Time')
set(gca,'FontSize',14,'FontWeight','b')
xlim([0,480])
%draw_opt.save_filename ='./fig';

figure(3)
legend(dose_str)
title(sti)
ylabel('IKK')
xlabel('Time')
set(gca,'FontSize',14,'FontWeight','b')
xlim([0,480])

a = nfkb_curves';

% for i =1:3
%     figure(i)
% 
% legend_str = {'NFkBtot = 0.2','0.1','0.05'}
% legend(legend_str)
% end