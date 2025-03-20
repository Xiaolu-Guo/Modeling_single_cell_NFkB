%run script for information theory for signaling codons
function [info, nfkb, all_dims, names_1D, sc_dims, names_sc_1D] = run_cc_traj_202306(local_flag, data_file)
P = mfilename('fullpath');
P2 = mfilename;
nfkb_name = [P(1:(length(P)-length(P2))) , data_file];


[nfkb, all_dims, sc_dims, names_1D, names_sc_1D] = loadnfkb_traj(nfkb_name); % run this locally, to be able to load files, then transfer nfkb.mat to server, then run this script to load into right variable names etc
disp('data loaded')
if ~local_flag
    % addpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information')
    % cd('/home/apeksha/channel_capacity/')
    % mkdir('all_stim')
    % cd('/home/apeksha/channel_capacity/all_stim')
    %cd('./data_for_example/all_stim/')
    [info] = inforun_traj(all_dims,  names_sc_1D);
    % save(strcat(save_name,'.mat'));
    % cd('../../')
end
end