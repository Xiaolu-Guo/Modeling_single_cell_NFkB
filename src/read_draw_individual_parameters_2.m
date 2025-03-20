%plot distribution
%
function [parameters,para_name,estimates] = read_draw_individual_parameters_2(proj_num,data_parent_path)


    [proj_path,~] = subfunc_get_proj_path_2(proj_num,data_parent_path);

% opts = detectImportOptions('airlinesmall.csv')
%  opts.DataLines=1;

    result_filepath= strcat(proj_path,'results/');


estimates_cell=readtable(strcat(result_filepath,'summary.txt'));%,opts);
para_simul_cell=readtable(strcat(result_filepath,'IndividualParameters/simulatedIndividualParameters.txt'));
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
