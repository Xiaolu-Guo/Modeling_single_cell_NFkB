function [para_val,estimates_params_sample] = sample_est_parameters_multi_receptor_2023(proj_num,Num_sample,data_parent_path,var_fold)

if length(Num_sample)<1
    error('wrong input: length of proj_num should be the same as Num_sample')
end

if nargin <4
    var_fold =1; % change the variance, change the heterogeneity?
end


if length(proj_num)~= length(Num_sample)
    Num_sample = Num_sample*ones(size(proj_num));
end

% Num_sample_1 = round(Num_sample/2);
% Num_smaple_2 = Num_sample - Num_sample_1;

for i_ligand =1:length(proj_num)
    [para_sample_ligand_large{i_ligand},~] = sample_est_parameters_2023(proj_num(i_ligand),max(Num_sample)*1000,data_parent_path,var_fold);
    [para_sample_ligand{i_ligand},estimates_ligand{i_ligand}] = sample_est_parameters_2023(proj_num(i_ligand),Num_sample(i_ligand),data_parent_path,var_fold);
    
    index_core_ligand{i_ligand} = [find(strcmp(estimates_ligand{i_ligand}.name,'params66')), ...
        find(strcmp(estimates_ligand{i_ligand}.name,'params99')),...
        find(strcmp(estimates_ligand{i_ligand}.name,'params100')), ...
        find(strcmp(estimates_ligand{i_ligand}.name,'params101') ), ...
        find(strcmp(estimates_ligand{i_ligand}.name,'NFkB_cyto_init'))];
    index_rcp_ligand{i_ligand} = setdiff(1:size(para_sample_ligand{i_ligand},2),index_core_ligand{i_ligand},'stable');
    
end


clear estimates_params_sample

i_ligand =1;
index_core = 1:length(index_core_ligand{i_ligand});
estimates_params_sample.name(index_core) = estimates_ligand{i_ligand}.name(index_core_ligand{i_ligand});

for i_ligand =1 :length(proj_num)
    index_rcp{i_ligand} = (length(estimates_params_sample.name)+1):(length(estimates_params_sample.name)+length(index_rcp_ligand{i_ligand}));
    estimates_params_sample.name(index_rcp{i_ligand}) =...
        estimates_ligand{i_ligand}.name(index_rcp_ligand{i_ligand});
end

para_val = zeros(sum(Num_sample),length(estimates_params_sample.name));

j_para = 1;
for i_ligand = 1:length(proj_num)
    
    for i_para = 1:size(para_sample_ligand{i_ligand},1)
        para_val(j_para,index_core) = para_sample_ligand{i_ligand}(i_para,index_core_ligand{i_ligand});
        
        para_val(j_para,index_rcp{i_ligand}) = ...
            para_sample_ligand{i_ligand}(i_para,index_rcp_ligand{i_ligand});
        
        for j_ligand = 1:length(proj_num)
            if i_ligand ~= j_ligand
                diff = sum((para_sample_ligand{i_ligand}(i_para,index_core_ligand{i_ligand}) - para_sample_ligand_large{j_ligand}(:,index_core_ligand{j_ligand})).^2,2);
                para_val(j_para,index_rcp{j_ligand}) = ...
                    para_sample_ligand_large{j_ligand}(diff == min(diff),index_rcp_ligand{j_ligand});
            end
        end
        j_para = j_para+1;
    end
end

estimates_params_sample.NFkB_index = find(strcmp(estimates_params_sample.name,'NFkB_cyto_init'));
estimates_params_sample.non_NFkB_index = setdiff(1:length(estimates_params_sample.name),estimates_params_sample.NFkB_index,'stable');



%
% [para_sample_ligand1_large,~] = sample_est_parameters_2(proj_num(1),Num_smaple_2*1000,data_parent_path,var_fold);
% [para_sample_ligand1,estimates_ligand1] = sample_est_parameters_2(proj_num(1),Num_sample_1,data_parent_path,var_fold);
%
% [para_sample_ligand2_large,~] = sample_est_parameters_2(proj_num(2),Num_sample_1*1000,data_parent_path,var_fold);
% [para_sample_ligand2,estimates_ligand2] = sample_est_parameters_2(proj_num(2),Num_smaple_2,data_parent_path,var_fold);
%
%
% index_core_ligand2 = [find(strcmp(estimates_ligand2.name,'params66')), ...
%     find(strcmp(estimates_ligand2.name,'params99')),...
%     find(strcmp(estimates_ligand2.name,'params100')), ...
%     find(strcmp(estimates_ligand2.name,'params101') ), ...
%     find(strcmp(estimates_ligand2.name,'NFkB_cyto_init'))];
% index_rcp_ligand2 = setdiff(1:size(para_sample_ligand2,2),index_core_ligand2,'stable');
%
% index_core_ligand1 = [find(strcmp(estimates_ligand1.name,'params66')), ...
%     find(strcmp(estimates_ligand1.name,'params99')),...
%     find(strcmp(estimates_ligand1.name,'params100')), ...
%     find(strcmp(estimates_ligand1.name,'params101') ), ...
%     find(strcmp(estimates_ligand1.name,'NFkB_cyto_init'))];
% index_rcp_ligand1 = setdiff(1:size(para_sample_ligand1,2),index_core_ligand1,'stable');
%
% estimates_params_sample.name((length(index_core_ligand2)+length(index_rcp_ligand1)+1):(length(index_core_ligand2)+length(index_rcp_ligand1)+length(index_rcp_ligand2))) =...
%     estimates_ligand2.name(index_rcp_ligand2);

% 
% for i_para = 1:size(para_sample_ligand2,1)
%     j_para = i_para + size(para_sample_ligand1,1);
%     
%     para_val(j_para,1:length(index_core_ligand2)) = para_sample_ligand2(i_para,index_core_ligand2);
%     
%     diff = sum((para_sample_ligand1(i_para,index_core_ligand1) - para_sample_ligand2_large(:,index_core_ligand2)).^2,2);
%     para_val(j_para,(length(index_core_ligand2)+1):(length(index_core_ligand2)+length(index_rcp_ligand1))) = ...
%         para_sample_ligand1_large(diff == min(diff),index_rcp_ligand1);
%     
%     para_val(j_para,(length(index_core_ligand2)+length(index_rcp_ligand1)+1):(length(index_core_ligand2)+length(index_rcp_ligand1)+length(index_rcp_ligand2))) = ...
%         para_sample_ligand2(i_para,index_rcp_ligand2);
%     
% end
