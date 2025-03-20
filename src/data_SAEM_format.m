function [] = data_SAEM_format(data_info,params)

% Input ligand index vector, dose vector, cell numbers vector, timepts
% vector; must be the same length


% sti_vec.ligand_index_vec = [];

% for ii = 1:length(data_info.ligand_vec)
%     switch data_info.ligand_vec{ii}
%         case 'TNF'
%             sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 1];
%         case 'CpG'
%             sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 2];
%         case 'LPS'
%             sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 3];
%         case 'Pam3CSK'
%             sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 4];
%         case 'PolyIC'
%             sti_vec.ligand_index_vec = [sti_vec.ligand_index_vec 5];
%     end
% end

% sti_vec.dose_index_vec = cell2mat( data_info.dose_index_vec);
sti_vec.SAEM_ending_time_mins_vec  = cell2mat(data_info.SAEM_ending_time_mins_vec);

%% set the default value for

if (~isfield(sti_vec,'cell_num_vec')) || ((isfield(sti_vec,'cell_all_vec') && (sum(sti_vec.cell_all_vec == 1) == length(sti_vec.cell_all_vec ))))
    
    sti_vec.cell_num_vec=-ones(length(data_info.ligand_vec),1);
    for sti_index= 1: length(data_info.ligand_vec)
        if isfield(params,'sheetname')
            stimuli_info_tbl = get_stimuli_info_tbl(params.sheetname);
        else
            stimuli_info_tbl = get_stimuli_info_tbl();
        end    
        index = (stimuli_info_tbl.Ligand == char(data_info.ligand_vec{sti_index})) & (stimuli_info_tbl.dose_index == data_info.dose_index_vec{sti_index});
        sti_vec.cell_num_vec(sti_index) =  stimuli_info_tbl.num_cells(index);
    end
end

if ~isfield(sti_vec,'SAEM_ending_time_mins_vec')
    sti_vec.SAEM_ending_time_mins_vec = params.end_time_mins * ones(length(data_info.ligand_vec),1);
end

if isfield(params,'TimePtsWeight')
else
    params.TimePtsWeight = ones(max(sti_vec.SAEM_ending_time_mins_vec)/params.time_each_frame+1,sum(sti_vec.cell_num_vec));
end


%%



if isfile(params.write_SAEM_data_filename)
    delete(params.write_SAEM_data_filename);
end

file_to_write=fopen(params.write_SAEM_data_filename,'w'); %create the txt data file

fprintf(file_to_write,'%s \t %s \t %12s \t %12s \t %12s \t %12s \r\n','ID','TIME','AMT','ADM','Y','TimePtWeight'); %write the txt data file

cell_id=1;

time_stimuli=0.15;

for sti_index= 1: length(data_info.ligand_vec)
    
    if isfield(params,'sheetname')
        stimuli_info_tbl = get_stimuli_info_tbl(params.sheetname);
    else
        stimuli_info_tbl = get_stimuli_info_tbl();
    end
    
    index = (stimuli_info_tbl.Ligand == char(data_info.ligand_vec{sti_index})) & (stimuli_info_tbl.dose_index == data_info.dose_index_vec{sti_index});
    
        cell_num = stimuli_info_tbl.num_cells(index);

    
    stim_info.ADM = stimuli_info_tbl.ADM(index);
    stim_info.dose_val = stimuli_info_tbl.dose_val(index);
    stim_info.dose_str_SAEM = char(stimuli_info_tbl.dose_str_SAEMData(index));
    stim_info.data_name_rescale = strcat(data_info.ligand_vec{sti_index},'_',stim_info.dose_str_SAEM ,...
        '.xls');
    
    
    exp_data=readmatrix(strcat(params.save_path,stim_info.data_name_rescale));
    
    time_matrix=ones(cell_num,1) * (0:params.time_each_frame: sti_vec.SAEM_ending_time_mins_vec(sti_index));
    
    exp_data_wr=exp_data';
    
    time_matrix_wr=time_matrix';
    
    for i_cell=1:cell_num
        
        fprintf(file_to_write,'%d \t %d \t %12s \t %12s \t %12.5g \t %12.5g \r\n',cell_id,time_matrix_wr(1,i_cell),'.','.',exp_data_wr(1,i_cell),1); %write the txt data file
        
        fprintf(file_to_write,'%d \t %.2f \t %12.6g \t %12d \t %12s \t %12s \r\n',cell_id,time_stimuli,stim_info.dose_val,stim_info.ADM,'.','.'); %write the txt data file
        
        
        for k=2: (sti_vec.SAEM_ending_time_mins_vec(sti_index)/params.time_each_frame+1)
            
            fprintf(file_to_write,'%d \t %d \t %12s \t %12s \t %12.5g \t %12.5g \r\n',cell_id,time_matrix_wr(k,i_cell),'.','.',exp_data_wr(k,i_cell),params.TimePtsWeight(k,i_cell));
        end
        
        cell_id=cell_id+1;
    end
    
end

fclose(file_to_write);

% data_info.dose_str_vec = cell(1,1);
% for i_sti = 1:length(data_info.ligand_vec)
% dose_str_all = get_stimuli_info_dose_str_all(data_info.ligand_vec{i_sti});
% %         info_dose_str_vec{i_sti} = string( dose_str_all);
%       data_info.dose_str_vec{i_sti} =  char(dose_str_all(data_info.dose_index_vec{i_sti}));
% end
%
% if nargin<2
%     data_info = subfunc_get_data_info(data_info);
% else
%     data_info = subfunc_get_data_info(data_info,params.sheetname);
%
% end


% params.TimePtsWeight;
% params.write_filename
% params.time_mins

% output: write SAEM format file

% example:


% % sti_vec.ligand_index_vec = [1,1,1];
% sti_vec.dose_index_vec = [1,2,3];
%
% sti_vec.ligand_index_vec = [1];
% sti_vec.dose_index_vec = [3];
% sti_vec.cell_num_vec = [30]; %or not define
% sti_vec.SAEM_ending_time_mins_vec = [480];

% sti_vec.ligand_index_vec = [1,1,2];
% sti_vec.dose_index_vec = [2,3,1];
% sti_vec.cell_num_vec = [30,30,30]; %or not define
% sti_vec.SAEM_ending_time_mins_vec = [480,480,480];

% calledin data_main.m

%% check whether the same length
%if

