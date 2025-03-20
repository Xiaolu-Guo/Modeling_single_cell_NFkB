%plot distribution
%
function [parameters,para_name,estimates] = read_est_individual_tbl(proj_num,data_parent_path)


[proj_path,~] = subfunc_get_proj_path(proj_num,data_parent_path);

% opts = detectImportOptions('airlinesmall.csv')
%  opts.DataLines=1;

result_filepath= strcat(proj_path,'results/');


estimates_cell=readtable(strcat(result_filepath,'summary.txt'));%,opts);
% para_simul_cell=readtable(strcat(result_filepath,'IndividualParameters/simulatedIndividualParameters.txt'));
para_pred_cell=readtable(strcat(result_filepath,'IndividualParameters/estimatedIndividualParameters.txt'));

%range=[0,1];

para_num=find(strcmp(estimates_cell.Effects,'Deviation'))-1;

minmax_range=find(strcmp(estimates_cell.Fixed,'min'));
minmax_range=minmax_range+1:minmax_range+para_num;

estimates.min=estimates_cell.x____________________________(minmax_range);
estimates.max=estimates_cell.Var7(minmax_range);
estimates.mean=estimates_cell.x____________________________(1:para_num);
estimates.std=estimates_cell.x____________________________(para_num+2:para_num+para_num+1);
%estimates.params=estimates_cell.Fixed(1:para_num);
estimates.params = strrep(estimates_cell.Fixed(1:para_num),'_pop','');
estimates.params = strrep(estimates.params,'arams','');

%x_mean=0;


draw_para_num=length(estimates.params);

for ii=1:draw_para_num


    para_name{ii} = estimates.params{ii};% para_pred_cell.Properties.VariableNames{1+para_num*2+ii};
    %  histogram(para_pred_cell.(para_pred_cell.Properties.VariableNames{1+para_num*2+ii}),100,'Normalization','pdf')
    parameters(:,ii) = para_pred_cell.(para_pred_cell.Properties.VariableNames{1+para_num*2+ii});
    
end

sim_data_tbl = SimDataTblInitialize();

    sim_data_tbl.parameter_module(i_sim_data,i_para_name) = parameter_module{i_para_name};%has t obe change!!!!
sim_data_tbl.parameter_reac_num(i_sim_data,i_para_name) = reac_num(i_para_name);
sim_data_tbl.parameter_para_num(i_sim_data,i_para_name) = para_num(i_para_name);

sim_data_tbl.parameter_value(i_sim_data,i_para_name) = sim_info.parameter_value_vec(i_para_name, i_para_val);
%  sim_data.parameter_i_para{i_sim_data} = i_para_vec; %delete
sim_data_tbl.ligand(i_sim_data,1) = sim_info.ligand;
sim_data_tbl.dose_str(i_sim_data,1) = sim_info.dose_str{i_dose};
sim_data_tbl.dose_val(i_sim_data,1) = sim_info.dose_val{i_dose};
% might be changed for different genotype

% traj
sim_data_tbl.species(i_sim_data,1) = data_info.species_outputname{i_species_outputname};
sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) = 0; % initialize
sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) = sim_data_tbl.trajectory(i_sim_data,1:size(Toclear_curve.(species),2)) + Toclear_curve.(species)(i_dose,:);
    
