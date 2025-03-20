% TNF SIMULATION
%p_mod = [idx, all_x(:,idx(1))];
% run('../NFkB_common/exp_info_initialization.m')
addpath('../NFkB_common/');
    stimuli_info_tbl = get_stimuli_info_tbl();

addpath('../NFkB_Sensitive_module_parameter/');
names = {'IKK','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

sti = 'TNF';

% doses = [1 3.3 10 33 100];
            % doses = cell2mat(dose_val_all{strcmp(ligand_all,sti)});
            doses = stimuli_info_tbl.dose_val((stimuli_info_tbl.Ligand == sti));
            % dose_scale = 1/dose_scale_all{strcmp(ligand_all,sti)}; % Convert to uM (from nM)
            dose_scale = 1/get_stimuli_info_dose_scale(sti);
%             dose_field = dose_all{strcmp(ligand_all,sti)};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added

options.params_size = [99,1;100,1;101,1];
options.params_value = [0.4,0.4,0.4];
[v0.PARAMS, v0.SPECIES] = nfkbInitialize();
options.v.PARAMS = v0.PARAMS;
options.v.SPECIES = v0.SPECIES;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate all doses (only need to equilibrate on first iteration)
output = [];

for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [t,x] = nfkbSimulate({sti,doses(i)*dose_scale},names, [], {},options);
    end
    output = cat(3,output,x);
end

ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));
for ii = 1:length(doses)
    plot(1:size(nfkb_curves,1),nfkb_curves(:,ii),'LineWidth',2);hold on
end
for ii =1:length(doses)
    legend_str{ii} = num2str(doses(ii));
end

legend(legend_str)
ylabel(sti)
xlabel('Time')
set(gca,'FontSize',14,'FontWeight','b')
%draw_opt.save_filename ='./fig';

a = nfkb_curves';
%
% %% - - - - - - - LPS SIMULATION (match doses)- - - - - - - -
% p_mod = [];
% doses = [0.33 3.3 33];
% dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa
% names = {'TLR4','TLR4LPS','TLR4LPSen','TRIF','MyD88','TRAF6','IKK','IkBat','NFkBn'};
% options = struct;
% options.DEBUG = 1;
% options.SIM_TIME = 48*60;
% % Simulate all doses (only need to equilibrate on first iteration)
% output = [];
% for i = 1:length(doses)
%     if isempty(output)
%         [t,x,simdata] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
%     else
%         options.STEADY_STATE = simdata.STEADY_STATE;
%         [~,x] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
%     end
%     output = cat(3,output,x);
% end
%
% %plotSpecies(output,names,doses);
%
% ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
% nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));
%
% %% - - - - - - - CpG SIMULATION (match doses)- - - - - - - -
% p_mod = [
% %88 1 4e-4 % Degradation of TLR9_N (sets excess amt based on ratio w/ C-terminius degrdataion, 4e-4)
% %92 1 3 % Degradation rate via TLR9_N
% %93 1 1.6e-3 % Constitutive degradation (unbound = 4e-4)
% %89 1 0.028
% ];
% doses = [10 33 100 330];
% dose_scale = 1/1000; % Convert to uM (from nM)
% names = {'CpG','CpG_en','TLR9','TLR9_CpG','TLR9_N','MyD88','TRAF6','IKK','NFkBn'};
% options = struct;
% options.DEBUG = 1;
% options.SIM_TIME = 8*60;
% % Simulate all doses (only need to equilibrate on first iteration)
% output = [];
% for i = 1:length(doses)
%     if isempty(output)
%         [t,x,simdata] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
%     else
%         options.STEADY_STATE = simdata.STEADY_STATE;
%         [~,x] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
%     end
%     output = cat(3,output,x);
% end
%
%
% ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
% nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));
%
% %% - - - - - - - poly(I:C) SIMULATION (match doses)- - - - - - - -
% p_mod = [];
% doses = 1000*[3.3 10 33 100];
% dose_scale = 1/5e6; % Convert to uM. PolyI:C molecular weight: 1000KDa(+)
% names = {'polyIC','polyIC_en','TLR3','TLR3_polyIC','TRIF','TRAF6','IKK','NFkBn'};
% options = struct;
% options.DEBUG = 1;
% options.SIM_TIME = 8*60;
% % Simulate all doses (only need to equilibrate on first iteration)
% output = [];
% for i = 1:length(doses)
%     if isempty(output)
%         [t,x,simdata] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
%     else
%         options.STEADY_STATE = simdata.STEADY_STATE;
%         [~,x] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
%     end
%     output = cat(3,output,x);
% end
%
% ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
% nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));
%
% %% - - - - - - - Pam3CSK SIMULATION (match doses)- - - - - - - -
% p_mod = [];
% doses = [1 3.3 10 33];
% dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
% names = {'CD14',  'Pam3CSK', 'CD14_P3CSK','TLR2','TLR2_P3CSK','MyD88','TRAF6','TAK1','IKK','NFkBn'};
% options = struct;
% options.DEBUG = 1;
% options.SIM_TIME = 8*60;
%
% % Simulate all doses (only need to equilibrate on first iteration)
% output = [];
% for i = 1:length(doses)
%     if isempty(output)
%         [t,x,simdata] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
%     else
%         options.STEADY_STATE = simdata.STEADY_STATE;
%         [~,x] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
%     end
%     output = cat(3,output,x);
% end
%
% ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
% nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

