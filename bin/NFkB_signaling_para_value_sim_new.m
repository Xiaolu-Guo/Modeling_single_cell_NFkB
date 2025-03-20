function [sim_data_tbl] = NFkB_signaling_para_value_sim(sim_info,data_info)

%% add this part for p100o
% for i_regulator = 1:size(R,1)
%         
%         options.v.init_val.IKK_off = R(i_regulator,1);
%         options.v.init_val.MyD88_off = R(i_regulator,2);
%         options.v.init_val.TRIF_off = R(i_regulator,3);
%         options.v.init_val.TRAF6_off = R(i_regulator,4);
%         options.v.init_val.TAK1_off = R(i_regulator,5);
%         
% %         if isempty(nfkb_dynamics_ligand)
%             [t,x] = nfkbSimulate_distr_regulators({sim_info.ligand{i},sim_info.doses{i}/sim_info.dose_scale{i}},names, [], {},options);


%% read parameters, check input, initialize

% sim_info.parameter_name = {'p1';'p2'};
% sim_info.parameter_value_vec = [xxxx; xxxxx];

if isfield(sim_info,'sim_time')
else
    sim_info.sim_time = 8*60;
end

sim_info_fields = {'ligand';'dose_str';'dose_val'};
% sim_info.ligand = {'sti'};
% sim_info.dose_str = {'dose1','dose2'};
% sim_info.dose_val = {[],[]};
for i_sim_info_fields = 1:length(sim_info_fields)
    if ~isfield(sim_info,sim_info_fields{i_sim_info_fields})
        error(strcat(sim_info_fields{i_sim_info_fields},'is not specified!'))
    end
end

data_info_fields = {'species_outputname';'species_composition'};
% data_info.species_outputname = {'nucNFkB','TNFR'};
% data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'}}; % must
% be r x 1, for each cell i must be ri x 1
for i_data_info_fields = 1:length(data_info_fields)
    if ~isfield(data_info,data_info_fields{i_data_info_fields})
        error(strcat(data_info_fields{i_data_info_fields},'is not specified!'))
    end
end
names = unique(vertcat(data_info.species_composition{:}));


% data_info.flag =
% data_info.type =
if ~isfield(data_info,'flag')
    data_info.flag = '';
end

if ~isfield(data_info,'type')
    data_info.type = '';
end


% initialize the sim_data
sim_data_tbl = SimDataTblInitialize();

reac_num = zeros(1,length(sim_info.parameter_name));
para_num = zeros(1,length(sim_info.parameter_name));
parameter_module = cell(1,length(sim_info.parameter_name));
for i_para_name = 1:length(sim_info.parameter_name)
    [reac_num(i_para_name),para_num(i_para_name)] = paraname2rpnum(sim_info.parameter_name{i_para_name});
    parameter_module{i_para_name} = get_para_module(reac_num(i_para_name));
    
end

dose_scale = get_stimuli_info_dose_scale(sim_info.ligand); % Convert to uM (from nM)


options0 = struct;
options0.DEBUG = 0;
options0.SIM_TIME = sim_info.sim_time;
[v0.PARAMS, v0.SPECIES] = nfkbInitialize();
options0.v.PARAMS = v0.PARAMS;
options0.v.SPECIES = v0.SPECIES;

if isfield(sim_info,'species_vec')
    for i_species_vec = 1:length(sim_info.species_vec.name)
        options0.v.init_val.(sim_info.species_vec.name{i_species_vec}) = sim_info.species_vec.val{i_species_vec};
    end
end


i_sim_data=1;

for i_para_val = 1:size(sim_info.parameter_value_vec,2)% 10%
    
    %options.v.PARAMS(reac_num,para_k)
    
    % Simulate all doses (only need to equilibrate on first iteration)
    
    
    for i_para_name = 1:length(sim_info.parameter_name)
        options0.v.PARAMS(reac_num(i_para_name),para_num(i_para_name)) = sim_info.parameter_value_vec(i_para_name, i_para_val);
    end
    
    output = [];
    
    
    for i_dose = 1:length(sim_info.dose_val)%
        if isempty(output)
            options=options0;
            [~,x,simdata] = nfkbSimulate({sim_info.ligand,sim_info.dose_val{i_dose}*dose_scale},names, [], {},options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({sim_info.ligand,sim_info.dose_val{i_dose}*dose_scale},names, [], {},options);
        end
        output = cat(3,output,x);
    end
    
    for i_species_outputname = 1:length(names)
        species = names{i_species_outputname};
        Toclear_curve.(species) = squeeze(output(:,strcmp(names,species),:));
        Toclear_curve.(species) = Toclear_curve.(species)';
    end
    
    for i_dose = 1:length(sim_info.dose_val)
        for i_species_outputname = 1:length(data_info.species_outputname)
            
            for i_para_name = 1:length(sim_info.parameter_name)
                sim_data_tbl.parameter_module(i_sim_data,i_para_name) = parameter_module{i_para_name};%has t obe change!!!!
                sim_data_tbl.parameter_reac_num(i_sim_data,i_para_name) = reac_num(i_para_name);
                sim_data_tbl.parameter_para_num(i_sim_data,i_para_name) = para_num(i_para_name);
                
                parameter_name_tmp = rpnum2paraname(reac_num(i_para_name), para_num(i_para_name));
                
                if i_para_name == 1
                    parameter_name_full =  parameter_name_tmp;
                elseif i_para_name >1
                    parameter_name_full = strcat(parameter_name_full,'-',...
                        parameter_name_tmp);
                end
                
                sim_data_tbl.parameter_fold_change(i_sim_data,i_para_name) = sim_info.fold_change_vec(i_para_name, i_para_val);
                sim_data_tbl.parameter_value(i_sim_data,i_para_name) = sim_info.parameter_value_vec(i_para_name, i_para_val);
            end
            sim_data_tbl.parameter_name(i_sim_data) = parameter_name_full;
            %  sim_data.parameter_i_para{i_sim_data} = i_para_vec; %delete
            sim_data_tbl.ligand(i_sim_data,1) = sim_info.ligand;
            sim_data_tbl.dose_str(i_sim_data,1) = sim_info.dose_str{i_dose};
            sim_data_tbl.dose_val(i_sim_data,1) = sim_info.dose_val{i_dose};
            % might be changed for different genotype
            sim_data_tbl.flag(i_sim_data,1) = data_info.flag;
            sim_data_tbl.type(i_sim_data,1) = data_info.type;
            
            % traj
            sim_data_tbl.species(i_sim_data,1) = data_info.species_outputname{i_species_outputname};
            sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) = 0; % initialize
            for i_species_composition = 1:length(data_info.species_composition{i_species_outputname})
                
                species = data_info.species_composition{i_species_outputname}{i_species_composition};
                sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) = sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) + Toclear_curve.(species)(i_dose,:);

            end
            
            i_sim_data = i_sim_data+1;
        end
        
    end
    
    
end


end
