%% Find qualifying replicates based off of codeword distribution distances
% Apeksha Singh 07/2020, can now select different cutoff value for distance
% betweeen distributions for each codeword
% By default, one replicate is kept at minimum for each experimental condition
% The column reference replicate in results tbl has a 1 for the replicate 
% that was used as the standard for comparison for its experimental condition
%% Change path to folder containing distance calculation files
cd 'G:\Shared drives\Signaling Systems Lab\Apeksha\violin_plot_dist_11_04'
%% Select the cutoff value for each codeword
% make sure to list codewords and cutoff values in same order
codewords = {'Duration', 'EarlyVsLate', 'OscVsNonOsc', 'PeakAmplitude', 'Speed', 'TotalActivity'};
cutoffs = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
%% Select distribution distance calculation method and chose to save results
% can view results in Workspace without saving
% options for distance method: 'Bhattacharyya','DistinctArea_nbinsFD',
% 'Earth_Mover_Wasserstein_clustered','Hellinger', 'Jensen_Shannon',
% 'Kolmgorov_Smirnov', 'Kullback_Leibler'
distance_method = 'Jensen_Shannon';
save_results = false;
results = get_qual_reps(codewords, cutoffs, distance_method, save_results);
%% Plotting Functionality
% allows user to see how many replicates qualify while one cutoff value is
% varied and the others are fixed to the values in cutoffs above
create_plot = false;
% from codewords list above, enter index of codeword to vary (ex:
% 'Duration' = 1) and select range of values to explore
select_codeword = 1;
vary_codeword_cutoff = [0:0.01:1];
if create_plot
    num_qual_reps = zeros(size(vary_codeword_cutoff));
    for i = 1:length(vary_codeword_cutoff)
        cutoffs(select_codeword) = vary_codeword_cutoff(i);
        results = get_qual_reps(codewords, cutoffs, distance_method, false);
        num_qual_reps(i) = size(results.tbl, 1);
    end
    plot(vary_codeword_cutoff, num_qual_reps)
    xlabel(['value of ' codewords{select_codeword} ' cutoff value'])
    ylabel('number of qualifying replicates')
end
%% Code for finding qualified replicates as described above
function results= get_qual_reps(codewords, cutoffs, distance_method, save_results)
%get qualified replicates
    preconditions = {'control_', 'preIFNB(100U)_', 'preIFNg(10ng)_', 'preIL10(20ng)_', 'preIL4(10ng)_', 'preIL13(50ng)_'};
    stimulations = {'TNF', 'R848', 'Pam3', 'FSL1', 'PolyIC', 'CpG', 'LPS', 'FLA'};
    num_tot_reps = 0;
    qualified_reps_labels = [];
    qualified_reps_id = [];
    ref_index = [];
    conditions_list = [];
    for i = 1:length(stimulations)
        switch stimulations{i}
            case 'TNF'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('1ngTNF_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'R848'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('1ugR848_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'Pam3'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('100ngPam3CSK_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'FSL1'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('3ngFSL1_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'PolyIC'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('50ugHMWpolyIC_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'CpG'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('CpG100nM_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'LPS'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('LPS10ng_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
            case 'FLA'
                preconditions_lst = preconditions;
                [qualified_reps_id2add, qualified_reps_labels2add, ref_index2add, conditions_list2add, num_reps2add] = get_qualified_reps_for_stim('10ngFLA_distribution_distances', preconditions_lst, codewords, cutoffs, distance_method);
                num_tot_reps = num_tot_reps+num_reps2add;
                qualified_reps_labels = [qualified_reps_labels; qualified_reps_labels2add];
                qualified_reps_id = [qualified_reps_id; qualified_reps_id2add];
                ref_index = [ref_index; ref_index2add];
                conditions_list = [conditions_list; conditions_list2add];
        end
    end
    results = struct;
    results.number_total_replicates= num_tot_reps;
    results.number_experimental_conditions=numel(unique(conditions_list));
    results.tbl = table(conditions_list, qualified_reps_id, qualified_reps_labels, ref_index, 'VariableNames', {'experimental_condition', 'experimental_id', 'experimental_label', 'reference_replicate'});
    if save_results
        save('results', 'results')
    end
end
function [qualified_reps_id, qualified_reps_labels, ref_index, conditions_list, num_tot_reps] = get_qualified_reps_for_stim(dist_file_name, preconditions_list, codewords, cutoffs, distance_method)
    num_tot_reps = 0;
    qualified_reps_id = [];
    qualified_reps_labels = [];
    ref_index = [];
    conditions_list = [];
    for j = 1: length(preconditions_list)            
        condition = [string([preconditions_list{j} dist_file_name])];
        load(condition(1), 'distance_calculation')
        expts = distance_calculation.expt_info;
        num_tot_reps = num_tot_reps + size(expts, 1);
        distances = distance_calculation.(distance_method);
        [ref_rep, kept_reps] = find_best_ref(distances, size(expts, 1), codewords, cutoffs);
        id = expts.experimental_id;
        labels = expts.experimental_label;
        qualified_reps_id = [qualified_reps_id; id(ref_rep); id(kept_reps)];
        qualified_reps_labels = [qualified_reps_labels; labels(ref_rep); labels(kept_reps)];
        ref_index = [ref_index; 1; zeros(length(kept_reps),1)];
        conditions_list = [conditions_list; repmat(condition, length(kept_reps) +1, 1)];
    end
end
function [ref_num, kept_reps] = find_best_ref(distances, num_refs, codewords, cutoffs)
    num_rep_kept_per_ref = zeros(1, num_refs);
    kept_reps_per_ref = {};
    row_info = distances.(codewords{1})(:, 1:2);
    for i=1:length(num_refs)
        [r, ~] = find(row_info==i);
        ref_rows = row_info(r, :);
        for j=1:length(codewords)
            codeword_distances = distances.(codewords{j});
            ref_dists = codeword_distances(r, 3);
            [delete_r, ~] = find(ref_dists >cutoffs(j));
            if ~isempty(delete_r)
                delete_r = unique(delete_r);
                ref_rows(delete_r, :) = repmat([0 0], length(delete_r), 1);
            end
        end
        [r, ~] = find(ref_rows == 0);
        ref_rows(unique(r), :) = [];
        num_rep_kept_per_ref(i) = size(ref_rows, 1);
        kept_reps_per_ref = [kept_reps_per_ref, ref_rows];
    end
    [~, ref_num] = max(num_rep_kept_per_ref);
    if ~isempty(kept_reps_per_ref)
        kept_reps = kept_reps_per_ref{ref_num};
        kept_reps = kept_reps(kept_reps~=ref_num);
    else
        kept_reps = [];
    end
end           