function [reac_num, para_num] = paraname2rpnum(parameter_name)

para_name = replace(parameter_name, 'params','');
pnum = str2num(replace(para_name, 'n',' '));
switch length(pnum)
    case 1
        reac_num = pnum(1);
        para_num = 1;
    case 2
        reac_num = pnum(1);
        para_num = pnum(2);
    otherwise
        error('wrong input for paramter name')
end

end