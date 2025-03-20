function [proj_path,proj_name] = subfunc_get_proj_path_2(proj_num,proj_parent_path,proj_type)
% relative path to Postdoc projects


if nargin<1 || nargin>3
    error('wrong input number!')
end

switch nargin
    case 1
        proj_type = 'XGESN';
        proj_path = './NFkB_para_estm_project/SAEM_proj_2022_2/';
        
    case 2
        proj_type = 'XGESN';
        proj_path = proj_parent_path;
    case 3
        proj_path = proj_parent_path;

end


proj_name = strcat(proj_type,sprintf( '%03d', proj_num));


proj_path = strcat(proj_path,proj_name,'/');

end

