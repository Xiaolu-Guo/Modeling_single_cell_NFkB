function [collect_feature_vects,metric_structs] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metric_structs) %, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% codeword distribution for different replicates/experimental id's
% utilizes the following to generate violn plots: [] = violin(vects, places, varargin)
% code can also generate distance metric between generated distribuitions
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% INPUTS (required)
% conditions = ["..." "..."] %either expt condition or date

% INPUT PARAMETERS (optional; specify with name-value pairs)
% saveFeatureData: is true if you want to save codeword data, these
%   calculations only hold for comparisons within experimental condition
%   (replicates). Once final dataset produced, can get true codeword values
% plotViolin: is true if want to generate plot
% saveVPlt: is true if wat to save violin plots
% date: is true or false depending if want to grab experiments by expt date or
%   expt condition respectively
% calcResponders: is true or false depending if want to only calculate
%   features for responders or all trajectories respectively]
% calcDistance: is true if want to calculate distance between distributions
% saveDistance: is true if want to save distributions distances
% idFile: path to excel sheet with experimental ids
% CodewordListFile: path to excel sheet with codewords to be calculated
% (describes which features compose them too)
% CodewordListFileSheet: specifies which sheet of excel file to read
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% optional parameters
% p = inputParser;
% addParameter(p,'plotViolin',true, @islogical);
% addParameter(p, 'saveVPlt', true, @islogical);
% addParameter(p, 'saveFeatureData', false, @islogical);
% addParameter(p,'date',false, @islogical);
% addParameter(p,'calcResponders',true, @islogical);
% addParameter(p,'calcDistance',true,@islogical);
% addParameter(p,'saveDistance',true,@islogical);
% addParameter(p, 'idFile', 'C:\Users\apeks\MACKtrack\new_code_2020\Polarization_Stimuli.xlsx');
% addParameter(p, 'CodewordListFile', 'C:\Users\apeks\MACKtrack\new_code_2020\FeatureList.xlsx');
% addParameter(p, 'CodewordListFileSheet', 'Codewords');
% parse(p, varargin{:})
%% get revelant experimental id's and labels according to input condition
% idFile = parameters.idFile;
% idTable = tableread(idFile, 'Sheet', expt_excel_sheet);
% idList = table2cell(idTable(:, 1));
% labelList = table2cell(idTable(:, 2));
% true_labels = table2cell(idTable(:, 3));
%find id number of pertinent replicates
% IDs={1;2;3;4;5;6};

% switch nargin
%     case 2
%         data_label=vis_data_field;
%     case 1
%         vis_data_field=fieldnames(data);
% end


ids = 1:length(metric_structs);

% labels=[];
% for i_data_label =1:length(data_label)
%     labels = [labels string(data_label{i_data_label})];
% end

%% calculate features of replicates
pred_prameters = set_pred_parameters();
parameters = set_parameters();
% metric_structs = cell(1, length(ids));

collect_feature_vects = struct;

collect_feature_vects.info_ligand = cell(2,1);
collect_feature_vects.info_dose_str = cell(2,1);
collect_feature_vects.info_data_type = cell(2,1);

for i_ids=1:length(metric_structs)
        % metric_structs{ i_ids} = computeFeatures_pred_data(data.(vis_data_field{i_data_field}){i_sti}, pred_prameters);
        collect_feature_vects.info_ligand{i_ids} = data_info.info_ligand{i_ids};
    collect_feature_vects.info_dose_str{i_ids} = data_info.info_dose_str{i_ids};
    collect_feature_vects.info_data_type{i_ids} = data_info.data_label{i_ids};

end

FeatureListFile = parameters.CodewordListFile;
FeatureListTable = readtable(FeatureListFile, 'Sheet', parameters.CodewordListFileSheet);
FeatureList = table2cell(FeatureListTable(:,1));
FeatureSpecifiers = table2cell(FeatureListTable(:,2));
CodewordList = table2cell(FeatureListTable(:,3));

%% generate violin plots of codewords
nanzscore = @(x)(x-nanmean(x, 1))./nanstd(x, 0, 1);

codewords = unique(CodewordList);


for i = 1:length(codewords)
    vects = cell(1, length(ids));
    for j= 1:length(FeatureList)
        if strcmp(CodewordList{j}, codewords{i})
            for k=1:length(ids)
                metric_struct = metric_structs{k};
                feature = metric_struct.(FeatureList{j});
                if abs(FeatureSpecifiers{j}) ~= 0
                    feature = feature(:, abs(FeatureSpecifiers{j}));
                end
                if FeatureSpecifiers{j} < 0
                    feature = feature*-1;
                end
                if parameters.calcResponders
                    feature = feature(logical(metric_struct.responder_index));
                end
                vects{k} = [vects{k} feature];
            end
        end
    end
    if size(vects{1}, 2) > 1  %%% for codewords that are composed of more than one metric-->take z score and then take average
        vects_zscore=vects(1:2:end);
        all_data_zscore = cell2mat(vects_zscore(:));
        mean_zscore= nanmean(all_data_zscore,1);
        std_zscore=nanstd(all_data_zscore);
        all_data = cell2mat(vects(:));
        all_data = (all_data-mean_zscore)./std_zscore;
        
        count = 1;
        for k = 1:length(ids)
            vects{k} = nanmean(all_data(count:(size(vects{k}, 1)+count-1), :), 2);
            count = count + size(vects{k}, 1);
        end
    end
    all_data = cell2mat(vects(:)); %%% normalize each codeword
    all_data = (all_data - nanmean(all_data,1))./nanstd(all_data);
    % all_data = (all_data-prctile(all_data,0.5))/(prctile(all_data,99.5)-prctile(all_data,0.5)); %normalized except those extrema
    count = 1;
    for k = 1:length(ids)
        vects{k} = all_data(count:(size(vects{k}, 1)+count-1), :);
        count = count + size(vects{k}, 1);
    end
    collect_feature_vects.(codewords{i}) = vects';
end

end

function pred_parameters = set_pred_parameters()
pred_parameters.baseline=0;
pred_parameters.FramesPerHour=12;
pred_parameters.thresholdFrames=5; %number of consecutive frames cell needs to pass activity threshold to be considered a responder

pred_parameters.OnThresh=3;%addParameter(p,'onThresh', 3); %% changed from foldThresh/inputThresh; old default 1.25 -AS 6/18/2020
pred_parameters.offPad=12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)

pred_parameters.MinLifetime=100;   %   final frame used to filter for long-lived cells (default = 100)
pred_parameters.Delay=1;%no delay. the delayed frame, which frame to start?
%[metrics,aux, graph] = nfkbmetrics_pred_data(pred_data_SAEM,pred_parameters)    ;

pred_parameters.FeatureListFile=strcat('','FeatureList_v2.xlsx');%
%'/Users/admin/Documents/my document/Postdoc projects/NFkB_para_estm_project/SAEM_cas/TNF_ExpNoblock/Untitled/FeatureList.xlsx';
pred_parameters.additional_features=cell.empty;
end

function parameters = set_parameters()
parameters.plotViolin=true;
parameters.saveVPlt=true;
parameters.saveFeatureData = true;
parameters.date = false;
parameters.calcResponders = true;
parameters.calcDistance = true;
parameters.saveDistance = true;
%parameters.idFile = 'C:\Users\apeks\MACKtrack\new_code_2020\Polarization_Stimuli.xlsx';

parameters.CodewordListFile=strcat('','FeatureList_v2.xlsx');%./
parameters.CodewordListFileSheet = 'Codewords';
parameters.calcDistance = false;
end