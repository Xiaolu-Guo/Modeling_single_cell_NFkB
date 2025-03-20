function cellstr_vec = categorical2cellstr(categorical_vec)

cellstr_vec = cell(size(categorical_vec));
for i_row = 1:size(categorical_vec,1)
    for i_column = 1:size(categorical_vec,2)
    cellstr_vec{i_row,i_column} = char(categorical_vec(i_row,i_column));

    end
end