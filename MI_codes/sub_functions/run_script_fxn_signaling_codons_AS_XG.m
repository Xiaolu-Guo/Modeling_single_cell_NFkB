%run script for information theory for signaling codons
function [info, sc_info, nfkb, all_dims, names_1D, sc_dims, names_sc_1D] = run_script_fxn_signaling_codons_AS_XG(local_flag, data_file)
P = mfilename('fullpath');
P2 = mfilename;
nfkb_name = [P(1:(length(P)-length(P2))) , data_file];


[nfkb, all_dims, sc_dims, names_1D, names_sc_1D] = loadnfkb_signaling_codons_AS_XG(nfkb_name); % run this locally, to be able to load files, then transfer nfkb.mat to server, then run this script to load into right variable names etc
disp('data loaded')
if ~local_flag
    % addpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information')
    % cd('/home/apeksha/channel_capacity/')
    % mkdir('all_stim')
    % cd('/home/apeksha/channel_capacity/all_stim')
    %cd('./data_for_example/all_stim/')
    [info, sc_info] = inforun_singaling_codons_AS(all_dims, sc_dims, names_sc_1D);
    % save(strcat(save_name,'.mat'));
    % cd('../../')
end
end