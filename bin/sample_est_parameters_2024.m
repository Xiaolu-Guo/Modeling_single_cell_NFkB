function [para_val,estimates_params_sample] = sample_est_parameters_2024(proj_num,Num_sample,data_parent_path,var_fold)

if nargin <4
    var_fold =1; % change the variance, change the heterogeneity? 
end

[proj_path,~] = subfunc_get_proj_path_2(proj_num,data_parent_path,'XGESN2023');
% subfunc_get_proj_path_2(proj_num,data_parent_path);

% opts = detectImportOptions('airlinesmall.csv')
%  opts.DataLines=1;

result_filepath= strcat(proj_path,'results/');


estimates_cell=readtable(strcat(result_filepath,'summary.txt'));%,opts);
estimates_pop = readtable(strcat(result_filepath,'populationParameters.txt'));%,opts);
para_num=find(strcmp(estimates_cell.Effects,'Deviation'))-1;

minmax_range=find(strcmp(estimates_cell.Fixed,'min'));
minmax_range=minmax_range+1:minmax_range+para_num;

estimates.min=estimates_cell.x____________________________(minmax_range);
estimates.max=estimates_cell.Var7(minmax_range);
estimates.mean=estimates_cell.x____________________________(1:para_num);

estimates.mean=estimates_pop.value(1:para_num);
estimates.std = estimates_pop.value(para_num+1:2*para_num);
% estimates.std=estimates_cell.x____________________________(para_num+2:para_num+para_num+1);
%estimates.params=estimates_cell.Fixed(1:para_num);
estimates.params = strrep(estimates_pop.parameter(1:para_num),'_pop','');

opts1 = detectImportOptions('parameter_setting.xlsx','Sheet','param_setting_2023');
for jj = 1: length(opts1.VariableNames)
    opts1 = setvartype(opts1, opts1.VariableNames{jj}, 'char');
end

parameter_setting = readtable('parameter_setting.xlsx',opts1);
opts2 = detectImportOptions('parameter_setting.xlsx','Sheet','para_default');
for jj = 1: length(opts2.VariableNames)
    opts2 = setvartype(opts2, opts2.VariableNames{jj}, 'char');
end
default_parameter_setting = readtable('parameter_setting.xlsx',opts2);



estimates_params_sample.name = estimates.params(~strcmp(estimates.params,'shift'))';
estimates_params_sample.mean = estimates.mean(~strcmp(estimates.params,'shift'))';
estimates_params_sample.std = estimates.std(~strcmp(estimates.params,'shift'))';
estimates_params_sample.corr = ones(length(estimates_params_sample.name));

for i_para = 1:length(estimates_params_sample.name)
    
    index_para = strcmp(estimates_params_sample.name{i_para},parameter_setting.parameter);
    %     index_para_default = strcmp(est_info_proj.model.est_para{ii},default_parameter_setting.parameter);
    if ~sum(index_para)
        index_para_default = strcmp(estimates_params_sample.name{i_para},default_parameter_setting.parameter);
        if ~sum(index_para_default)
            error(strcat('no such parameter ',estimates_params_sample.name{i_para}));
        end
        min_val = str2double( default_parameter_setting.min_val{index_para_default});
        max_val = str2double( default_parameter_setting.max_val{index_para_default});
    else
        
        min_val = parameter_setting.min_val{index_para};
        
        if isempty(min_val)
            min_val = str2double(parameter_setting.Value{index_para})*0.1;
        else
            min_val = str2double(min_val);
        end
  
        max_val = parameter_setting.max_val{index_para};
        
        if isempty(max_val)
            max_val = str2double(parameter_setting.Value{index_para})*10;
        else
            max_val = str2double(max_val);
        end
    end
    
    estimates_params_sample.min(i_para) =  min_val;
    estimates_params_sample.max(i_para) = max_val;
    
    for j_para = 1:length(estimates_params_sample.name)
        if i_para~=j_para
            corr_name1 = strcat('corr_',estimates_params_sample.name(i_para),'_',estimates_params_sample.name(j_para));
            corr_name2 = strcat('corr_',estimates_params_sample.name(j_para),'_',estimates_params_sample.name(i_para));
            index_corr = strcmp(estimates_pop.parameter,corr_name1)|strcmp(estimates_pop.parameter,corr_name2);
            if sum(index_corr)~=1
                error('there is no corresponding correlationship')
            else
                estimates_params_sample.corr(i_para,j_para) = estimates_pop.value(index_corr);
            end
        end
    end
end

%% sampling

mu =       estimates_params_sample.mean; % mean value
sigma =    estimates_params_sample.std * var_fold; %standard deviation of the ramdom effects
corr_mat = estimates_params_sample.corr;
para_min = estimates_params_sample.min;
para_max = estimates_params_sample.max;

sigma_mat = diag(sigma)*corr_mat*diag(sigma);
mu_log = log((mu - para_min)./(para_max-mu));

N_cells = Num_sample;

R = mvnrnd(mu_log,sigma_mat,N_cells);

para_val = (exp(R).*(ones(N_cells,1)*para_max) + ones(N_cells,1)*para_min)./(1+exp(R));
