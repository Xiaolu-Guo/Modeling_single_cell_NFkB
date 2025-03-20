function [para_sample_multi_ligand,estimates_2ligand]= sample_fitting_multi_receptor_2023_11(proj_num,Num_sample,data,data_all,sim_info,fitting_dose)

%3 ligands or more 5 ligands
old_str = {'p','p','p','p','p','p','p','p','p','p'};
new_str = {'params','params','params','params','params','params','params','params','params','params'};


data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

i_data = [];
for i_sti = 1:length(proj_num)
    
    switch fitting_dose
        %         case 'middose'
        %             i_data = find(strcmp(data.info_ligand,sim_info.ligand),1,'last')-1;
        %             para_val_ele = data.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
        
        %         case 'lowdose'
        %             i_data = find(strcmp(data.info_ligand,sim_info.ligand),1);
        %             para_val_ele = data.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
        
        %         case 'highdose'
        %             i_data = find(strcmp(data.info_ligand,sim_info.ligand),1,'last');
        %             para_val_ele = data.parameters_mode_nan{i_data}  ; % [cellnumber x parameter]
        
        case 'alldose'
            i_data_0 = find(strcmp(data_all.info_ligand,sim_info.ligand{i_sti}));

%             i_data = find(strcmp(data.info_ligand,sim_info.ligand));
            para_val_ele{i_sti} = [];
            for i_i_data= 1:length(i_data_0)
                para_val_ele{i_sti}  = [para_val_ele{i_sti} ;data_all.parameters_mode_nan{i_data_0(i_i_data)} ] ; % [cellnumber x parameter]
            end
            i_data(i_sti) = find(strcmp(data.info_ligand,sim_info.ligand{i_sti}),1)
        case 'originaldose'
            %  info_dose_str_cmp = cellfun(@replace,data.info_dose_str,...
            %     old_str_dose,...
            %     new_str_dose,...
            %     'UniformOutput',false);
            i_data(i_sti) = find(strcmp(data.info_ligand,sim_info.ligand{i_sti}) &...
                strcmp(data_dose_str, sim_info.dose_str{i_sti})) ;
            
            para_val_ele{i_sti} = data.parameters_mode_nan{i_data(i_sti)} ; % [cellnumber x parameter]
        otherwise
            error('wrong input for fitting_dose value!')
    end
    
    
    rpt_time = ceil(Num_sample*5/size(para_val_ele{i_sti},1));
    para_sample_ligand_total{i_sti} = [];
    for i_rpt = 1:rpt_time
        para_sample_ligand_total{i_sti} = [para_sample_ligand_total{i_sti};para_val_ele{i_sti}];
    end
    para_sample_ligand{i_sti} = para_sample_ligand_total{i_sti}(randperm(size(para_sample_ligand_total{i_sti},1),Num_sample),:);
    estimates_ligand{i_sti}.name = cellfun(@replace,data.para_name{i_data(i_sti)},old_str(1:length(data.para_name{i_data(i_sti)})),new_str(1:length(data.para_name{i_data(i_sti)})),'UniformOutput',false);
    NFkB_index = find(strcmp(estimates_ligand{i_sti}.name,'NFkB_cyto_init'));
    shift_index = find(strcmp(estimates_ligand{i_sti}.name,'shift'));
    parameter_index = setdiff(setdiff(1:length(estimates_ligand{i_sti}.name),NFkB_index,'stable'),shift_index,'stable');
    
    
end

para_sample_ligand_large = para_sample_ligand_total;

for i_ligand =1:length(proj_num)
    
    index_core_ligand{i_ligand} = [find(strcmp(estimates_ligand{i_ligand}.name,'params52n2')), ...
        find(strcmp(estimates_ligand{i_ligand}.name,'params99')),...
        find(strcmp(estimates_ligand{i_ligand}.name,'params101')), ...
        find(strcmp(estimates_ligand{i_ligand}.name,'NFkB_cyto_init'))];
    index_rcp_ligand{i_ligand} = setdiff(1:size(para_sample_ligand{i_ligand},2),index_core_ligand{i_ligand},'stable');
    
end


i_ligand =1;
index_core = 1:length(index_core_ligand{i_ligand});
estimates_2ligand.name(index_core) = estimates_ligand{i_ligand}.name(index_core_ligand{i_ligand});

for i_ligand =1 :length(proj_num)
    index_rcp{i_ligand} = (length(estimates_2ligand.name)+1):(length(estimates_2ligand.name)+length(index_rcp_ligand{i_ligand}));
    estimates_2ligand.name(index_rcp{i_ligand}) =...
        estimates_ligand{i_ligand}.name(index_rcp_ligand{i_ligand});
end

para_val = zeros(sum(Num_sample),length(estimates_2ligand.name));

j_para = 1;
for i_ligand = 1:length(proj_num)
    
    for i_para = 1:size(para_sample_ligand{i_ligand},1)
        para_val(j_para,index_core) = para_sample_ligand{i_ligand}(i_para,index_core_ligand{i_ligand});
        
        para_val(j_para,index_rcp{i_ligand}) = ...
            para_sample_ligand{i_ligand}(i_para,index_rcp_ligand{i_ligand});
        
        for j_ligand = 1:length(proj_num)
            if i_ligand ~= j_ligand
                diff = sum((para_sample_ligand{i_ligand}(i_para,index_core_ligand{i_ligand}) - para_sample_ligand_large{j_ligand}(:,index_core_ligand{j_ligand})).^2,2);
                min_index = find(diff == min(diff));
                min_index_1 = min_index(randperm(length(min_index),1));
                para_val(j_para,index_rcp{j_ligand}) = ...
                    para_sample_ligand_large{j_ligand}(min_index_1,index_rcp_ligand{j_ligand});
            end
        end
        j_para = j_para+1;
    end
end

estimates_2ligand.NFkB_index = find(strcmp(estimates_2ligand.name,'NFkB_cyto_init'));
shift_index = find(strcmp(estimates_2ligand.name,'shift'));
estimates_2ligand.non_NFkB_index = setdiff(setdiff(1:length(estimates_2ligand.name),estimates_2ligand.NFkB_index,'stable'),shift_index,'stable');

para_sample_multi_ligand = para_val;