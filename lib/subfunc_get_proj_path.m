function [proj_path,proj_name] = subfunc_get_proj_path(proj_num,proj_parent_path,proj_type)
% relative path to Postdoc projects


if nargin<1 || nargin>3
    error('wrong input number!')
end

switch nargin
    case 1
        proj_type = 'XGES';
        proj_path = './NFkB_para_estm_project/SAEM_proj_2022/';
        
    case 2
        proj_type = 'XGES';
        proj_path = proj_parent_path;
end


proj_name = strcat(proj_type,sprintf( '%04d', proj_num));


proj_path = strcat(proj_path,proj_name,'/');

end

