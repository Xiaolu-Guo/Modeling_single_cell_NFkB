load('example_data.mat')

%% data_example is an example of the data file.

% data_example.exp and data_example.pred stores the experimental data and
% corresponding predicted data (from modeling simulation).

% data_example.info_* stores the corresponding stimuli info, including
% ligand, doses, cell numbers.

% data_example.order is the decend order of experimental data based on the
% Amplitude, which can be changed by hand to based on other codon/dynamical
% features.



% which data type to be calculated for the codon
vis_data_field = {'sample'};%,'pred'};

% labels for the data
data_label = {'sample'};%,'sample'};

% calculate the codon and dynamical features
[collect_feature_vects,metrics] = calculate_codon(data_example,vis_data_field,data_label);%,  parameter

% collect_feature_vects.speed (Osc, Duration, etc) stores the corresponding
% codon values;
% metrics{i_sti}.* stores the correspoding dynmical feature values.
