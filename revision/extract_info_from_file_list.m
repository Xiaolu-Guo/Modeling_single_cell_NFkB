% Given data_file_list, extract ligand, dose, genotype, and test info into separate variables

data_file_list = {'1000nMCpG_test_trajectories';
    '1000ngP3C4_test_trajectories';
    '100ngLPS_test_trajectories';
    '100ugPIC_test_trajectories';
    '10ngTNF_test_trajectories';
    'CpG_100nM_trajectories';
    'CpG_1uM_trajectories';
    'CpG_333nM_trajectories';
    'CpG_33nM_trajectories';
    'CpG_TNF_trajectories';
    'Ctrl_trajectories';
    'IkBaKO_TNF_3ng_trajectories';
    'IkBaMut_LPS_100ng_r1_trajectories';
    'IkBaMut_PolyIC_50ug_r1_trajectories';
    'IkBaMut_TNF_10ng_r1_trajectories';
    'IkBaWt_LPS_100ng_r1_trajectories';
    'IkBaWt_PolyIC_50ug_r1_trajectories';
    'IkBaWt_TNF_10ng_r1_trajectories';
    'IkBbeKO_LPS_10ng_trajectories';
    'IkBbeKO_TNF_3ng_trajectories';
    'LPS_0_3ng_trajectories';
    'LPS_100ng_trajectories';
    'LPS_10ng_trajectories';
    'LPS_1ng_trajectories';
    'LPS_333ng_trajectories';
    'LPS_33ng_r1_trajectories';
    'LPS_33ng_r2_trajectories';
    'LPS_33ng_r3_trajectories';
    'LPS_33ng_trajectories';
    'LPS_3ng_trajectories';
    'Pam3CSK4_100ng_trajectories';
    'Pam3CSK4_10ng_trajectories';
    'Pam3CSK4_1ng_trajectories';
    'Pam3CSK4_33ng_trajectories';
    'Pam3CSK4_3ng_trajectories';
    'PolyIC_100ug_trajectories';
    'PolyIC_10ug_trajectories';
    'PolyIC_33ug_trajectories';
    'PolyIC_3ug_trajectories';
    'TNF_0_3ng_trajectories';
    'TNF_10ng_trajectories';
    'TNF_1ng_trajectories';
    'TNF_33ng_trajectories';
    'TNF_3ng_trajectories'};

n_files = length(data_file_list);
ligand_vec = cell(n_files, 1);
dose_vec = cell(n_files, 1);
genotype_vec = cell(n_files, 1);
test_or_not_vec = cell(n_files, 1);
replicate_vec = cell(n_files, 1);

for i = 1:n_files
    data_file = data_file_list{i};
    
    % Extract if it's a test dataset or not
    if contains(data_file, 'test')
        test_or_not_vec{i} = 'test';
    else
        test_or_not_vec{i} = 'not_test';
    end
    
    % Extract genotype information
    if contains(data_file, 'IkBaKO')
        genotype_vec{i} = 'IkBaKO';
    elseif contains(data_file, 'IkBaMut')
        genotype_vec{i} = 'IkBaMut';
    elseif contains(data_file, 'IkBbeKO')
        genotype_vec{i} = 'IkBbeKO';
    else
        genotype_vec{i} = 'WT'; % Default Wild Type
    end
    
    % Extract ligand information
    ligands = {'CpG_TNF','CpG', 'P3C4', 'LPS', 'PIC', 'PolyIC', 'TNF', 'Pam3CSK4', 'Ctrl'};
    ligand_found = false;
    for j = 1:length(ligands)
        if contains(data_file, ligands{j})
            ligand_vec{i} = ligands{j};
            ligand_found = true;
            break;
        end
    end
    if ~ligand_found
        ligand_vec{i} = 'unknown'; % In case the ligand is not found in the list
    end
    
    % Extract dose information
    dose_pattern = '\d+[unm]g|\d+_\d+ng|\d+ng|\d+uM|\d+nM';
    dose_match = regexp(data_file, dose_pattern, 'match');
    if ~isempty(dose_match)
        dose_vec{i} = dose_match{1};
    else
        dose_vec{i} = 'unknown'; % Default if no dose is found
    end

    % Extract replicate information
    replicate_pattern = 'r\d+';
    replicate_match = regexp(data_file, replicate_pattern, 'match');
    if ~isempty(replicate_match)
        replicate_vec{i} = replicate_match{1};
    else
        replicate_vec{i} = 'r0'; % Default if no replicate information is found
    end
end

% Displaying the vectors (optional)
disp('Ligand Vector:'), disp(ligand_vec)
disp('Dose Vector:'), disp(dose_vec)
disp('Genotype Vector:'), disp(genotype_vec)
disp('Test or Not Vector:'), disp(test_or_not_vec)
