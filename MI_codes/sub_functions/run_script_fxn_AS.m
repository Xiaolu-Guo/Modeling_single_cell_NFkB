%run script for information theory on veqC
function [info, sc_info, nfkb, all_dims, names_1D, sc_dims, names_sc_1D] = run_script_fxn_AS(local_flag, polarization, save_name)
    [nfkb, all_dims, sc_dims, names_1D, names_sc_1D] = loadnfkb_AS(polarization); % run this locally, to be able to load files, then transfer nfkb.mat to server, then run this script to load into right variable names etc
    disp('data loaded')
    if ~local_flag
        addpath(genpath('/home/apeksha/channel_capacity'))
        cd(['/home/apeksha/channel_capacity/M', polarization])
        [info, sc_info] = inforun_AS(all_dims, names_1D, sc_dims, names_sc_1D);
        save(strcat(save_name,'.mat'));
    end
end