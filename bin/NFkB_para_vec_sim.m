function sim_data = NFkB_para_vec_sim(sim_info,data_info,sim_data)
% reac_num,para_num,i_para_vec,para_vec(i_para_vec)

reac_num =          sim_info.reac_num;
para_num =          sim_info.para_num;
para_vec =          sim_info.para_vec;
sti =               sim_info.stimuli;
dose_scale =        sim_info.dose_scale;
dose_field =        sim_info.dose_field;
parameter_module =  sim_info.parameter_module;
doses =             sim_info.doses;
fold_change_vec =   sim_info.fold_change_vec;
% sim_info.dose_val_all
% dose_vec = sim_info.dose_vec;

% data_info.species;
% data_info.save_file_path
% data_info.save_file_name

% Modified from sim_expamples

% updated on 04/19/2022

% sim_info.stimuli
% sim_info.parameter
% data_info.species
% data_info.save_file_path
% data_info.save_file_name

%% initialize
% if nargin <1
%     Module_sens = 'TNF';
% end


%% run the simulation


i_sim_data = length(sim_data.trajectory)+1;

names = data_info.species;

options0 = struct;
options0.DEBUG = 0;
options0.SIM_TIME = 8*60;
%
[v0.PARAMS, v0.SPECIES] = nfkbInitialize();

options0.v.PARAMS = v0.PARAMS;
options0.v.SPECIES = v0.SPECIES;

for i_para_vec = 1:length(para_vec)% 10%
    
    options0.v.PARAMS(reac_num,para_num) = para_vec(i_para_vec);
    
    
    %options.v.PARAMS(reac_num,para_k)
    
    % Simulate all doses (only need to equilibrate on first iteration)
    
        output = [];

        
        for i_dose = 1:length(doses)%
            if isempty(output)
                options=options0;
                [~,x,simdata] = nfkbSimulate({sti,doses(i_dose)*dose_scale},names, [], {},options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({sti,doses(i_dose)*dose_scale},names, [], {},options);
            end
            output = cat(3,output,x);
        end
        
        for i_species = 1:length(names)
            species = names{i_species};
            Toclear_curve.(species) = squeeze(output(:,strcmp(names,species),:));
            Toclear_curve.(species) = Toclear_curve.(species)';
        end
        
        for i_dose = 1:length(dose_field)
            for i_species = 1:length(names)
                species = names{i_species};
                sim_data.parameter_module{i_sim_data,1} = parameter_module;
                sim_data.parameter_reac_num{i_sim_data,1} = reac_num;
                sim_data.parameter_para_num{i_sim_data,1} = para_num;
                if para_num ==1
                    sim_data.parameter_name{i_sim_data,1} = strcat('params',num2str(reac_num));
                else
                    sim_data.parameter_name{i_sim_data,1} = strcat('params',num2str(reac_num),'n',num2str(para_num));
                end
                
                sim_data.parameter_fold_change{i_sim_data,1} = fold_change_vec(i_para_vec);
                sim_data.parameter_value{i_sim_data,1} = para_vec(i_para_vec);
               %  sim_data.parameter_i_para{i_sim_data} = i_para_vec; %delete
                sim_data.ligand{i_sim_data,1} = sti;
                sim_data.dose_str{i_sim_data,1} = dose_field{i_dose};
                sim_data.dose_val{i_sim_data,1} = dose_field{i_dose};
                sim_data.species{i_sim_data,1} = species;
                % might be changed for different genotype
                sim_data.flag{i_sim_data,1} = 'none';
                sim_data.type{i_sim_data,1} = 0;
                
                sim_data.trajectory{i_sim_data,1} = Toclear_curve.(species)(i_dose,:);
                i_sim_data = i_sim_data+1;
            end
            
        end
    
    
end



