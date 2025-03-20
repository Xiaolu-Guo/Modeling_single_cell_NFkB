function sim_data_tbl = SimDataTblInitialize()
sim_data_tbl = table();
sim_data_tbl.parameter_module = categorical();
sim_data_tbl.parameter_reac_num = zeros(0);
sim_data_tbl.parameter_para_num = zeros(0);
sim_data_tbl.parameter_name = categorical();
sim_data_tbl.parameter_fold_change = zeros(0);
sim_data_tbl.parameter_value = zeros(0);
sim_data_tbl.ligand = categorical();
sim_data_tbl.dose_str = categorical();
sim_data_tbl.dose_val = zeros(0);
sim_data_tbl.species = categorical();
sim_data_tbl.flag = categorical();
sim_data_tbl.type = categorical();
sim_data_tbl.trajectory = zeros(0);

end