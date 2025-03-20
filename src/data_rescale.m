function [] = data_rescale(rescale_factor,ligand_index,dose_index,params)
% updated on 11/24/2021
% rescale the data using the rescaling factor

load(params.exp_original_data_file);

sti_info = get_stimulus_info(ligand_index,dose_index,params.exp_info);

%% recale all the data

index = ((dataTbl.Ligand == sti_info.ligand_name) & (dataTbl.Dose == sti_info.dose));
data = dataTbl.time_series(index,:);

% for debugging
%histogram(max(data,[],2));hold on

data_rescale = data * rescale_factor;

% for debugging
% histogram(max(data_rescale,[],2));hold on


if isfolder(params.save_path)
else
    mkdir(params.save_path)
end

if isfile(strcat(params.save_path,sti_info.data_name_rescale))
     delete(strcat(params.save_path,sti_info.data_name_rescale));
end

writematrix(data_rescale,strcat(params.save_path,sti_info.data_name_rescale));


