function parameter_name = rpnum2paraname(reac_num, para_num)


if para_num ==1
    parameter_name = strcat('params',num2str(reac_num));
else
    parameter_name = strcat('params',num2str(reac_num),'n',num2str(para_num));
end

end