function stimuli_info_tbl = get_stimuli_info_tbl(sheetname)

if nargin<1
opts1 = detectImportOptions('stimuli_info.xlsx','Sheet','Ade_exp_data_2023');
else
    opts1 = detectImportOptions('stimuli_info.xlsx','Sheet',sheetname);

end

for jj = 1: length(opts1.VariableNames)
    opts1 = setvartype(opts1, opts1.VariableNames{jj}, 'char');
end
stimuli_info_tbl = readtable('stimuli_info.xlsx',opts1);

categorical_fields = {'Ligand';'dose_fold';'dose_str';'dose_str_SAEMData'};
double_fields = {'dose_val';'dose_scale';'num_cells';'ADM';'dose_index'};

for i_categorical_fields = 1:length(categorical_fields)
    stimuli_info_tbl.(categorical_fields{i_categorical_fields}) = ...
        categorical(stimuli_info_tbl.(categorical_fields{i_categorical_fields}) );
end

for i_double_fields = 1:length(double_fields)
    stimuli_info_tbl.(double_fields{i_double_fields}) = ...
        cellfun(@str2num,stimuli_info_tbl.(double_fields{i_double_fields}));

end
