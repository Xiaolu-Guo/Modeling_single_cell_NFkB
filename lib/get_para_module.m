function module = get_para_module(para_name_num)
% based on 'parameter_setting.xlsx' sheet 'param_setting'

if ~isnumeric(para_name_num)
    [para_name_num,~] = paraname2rpnum(para_name_num);
end

module_index = sum([26,46,52,65,67,76,84,94,101]<para_name_num)+1;

module_vec = {'core','LPS','coreup','TNF','core','Pam3CSK','polyIC','CpG','core'};

module = module_vec{module_index};

end