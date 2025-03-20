function [features, metrics]= computeFeatures_pred_data(pred_data,pred_parameters)
%function that takes experimental ID number and computes features (optionally given 'FeatureList' table) or in varargin 
%from basic metrics function and additional functions
%based of Stefanie's code (derived from Ade's computeFeatures function, Brooks' nfkbmetrics function, nfkbpeaks function, etc.
%% INPUT PARSING
% Create input parser object, add required params from function input
% p = inputParser;
% % Required: ID input
% valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
%     'ID input must be spreadsheet ID or full file path');
% addRequired(p,'id',valid_id);
% % Optional parameters to be passed to metrics function
% expectedFlags = {'on','off'};
% addParameter(p,'Verbose','on', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
% valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
%     'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
% addParameter(p,'ConvectionShift',1, valid_conv); %allows adjustment of convection shift (?)
% addParameter(p,'MinLifetime',97, @isnumeric); %allows adjustment of minimum lifetime (?); old default 130-AS 6/18/2020
% addParameter(p,'TrimFrame',100, @isnumeric);
% addParameter (p, 'onThresh', 3, @isnumeric); %sigma threshold for determining responders
% addParameter(p, 'Delay', 3 , @isnumeric)
% addParameter(p,'Frames',5); %old default 4 =AS 6/18/2020
% addParameter(p, 'FeatureListFile', 'C:\Users\apeks\MACKtrack\new_code_2020\FeatureList.xlsx') %provide file path for Excel table with list of feature to be computed
% addParameter(p, 'FeatureListFileSheet', 'Codewords') %sheet from previous excel to pull features from
% addParameter(p, 'additional_features', {}) %add in more features to calculate beyond what is listed on excel spreadsheet {'feat_name_1';....'feat_name_n'}
% parse(p,id, varargin{:})

%% get list of features to be calculated
if ~isempty(pred_parameters.FeatureListFile)
    FeatureListFile = pred_parameters.FeatureListFile;
    FeatureListTable = readtable(FeatureListFile, 'Sheet', 'FeatureList');%
    FeatureList = table2cell(FeatureListTable(:,1));
    FeatureList = [FeatureList; pred_parameters.additional_features];
else
    FeatureList = pred_parameters.additional_features;
end
%% get names of peak and signal features
%FeatureDescriptionFile = 'C:\Users\apeks\MACKtrack\new_code_2020\FeatureList.xlsx';
%FeatureDescriptionFile = '/Users/admin/Documents/my document/Postdoc projects/NFkB_para_estm_project/SAEM_cas/TNF_ExpNoblock/Untitled/FeatureList.xlsx';
FeatureDescriptionFile = pred_parameters.FeatureListFile;

PeakStatsTable = readtable(FeatureDescriptionFile, 'Sheet', 'PeakStats');
PeakStatsList = table2cell(PeakStatsTable(:,1));
SignalStatsTable = readtable(FeatureDescriptionFile, 'Sheet', 'SignalStats');% FeatureListFile
SignalStatsList = table2cell(SignalStatsTable(:,1));

%% call function for basic metrics
[metrics,aux, graph] = nfkbmetrics_pred_data(pred_data,pred_parameters);
%% calculate new features that rely on external functions
FramesPerHour = 12;
Delay = pred_parameters.Delay;
MinLifetime = pred_parameters.MinLifetime;
if Delay < 1
    error('delay must be greater than zero, to compute with delay = 0 use old computeFeatures.m')
end
%% get Peak Stat, Signal Stat, positive integrals, positive peaks, and fold change features
if any(ismember(FeatureList, PeakStatsList))
    peak_stats =get_peak_stats_new(metrics.time_series_no_base_ded, metrics.baseline, 'FramesPerHour', FramesPerHour);
end

if any(ismember(FeatureList, SignalStatsList))
    % sig_stats =get_sig_stats_new(metrics.time_series, 'FS', FramesPerHour);
    sig_stats =get_sig_stats_v2023(metrics.time_series, 'FS', FramesPerHour);
end

if any(ismember(FeatureList, {'integrals_pos','time2HalfMaxPosIntegral','max_pos_integral'}))
    endFrame = min(96, size(metrics.integrals,2));
    pos_integral_features = get_pos_integrals(metrics.integrals, FramesPerHour, endFrame);
end

if any(ismember(FeatureList, {'pos_pk1_amp', 'pos_pk1_time', 'pos_pk1_width', 'pos_pk1_prom', 'pos_pk1_height',...
        'pos_pk2_amp', 'pos_pk2_time', 'pos_pk2_width', 'pos_pk2_prom', 'pos_pk2_height', 'max_pos_pk1_speed'...
        'pos_pk2_ratio', 'pos_pk2_ratio_prom'}))
    pos_peak_features = get_positive_peaks(metrics.time_series_no_base_ded, metrics.baseline, FramesPerHour, MinLifetime);
end

if any(ismember(FeatureList, {'fold_change', 'max_fold_change', 'max_value','time2HalfMaxValue'}))
    % fold_change_features = get_fold_change_new(metrics.time_series_no_base_ded, metrics.baseline);
    fold_change_features = get_fold_change_new_20220701(metrics.time_series_no_base_ded, metrics.baseline,FramesPerHour);

end
%% intialize feature structure to return
features =struct;

%% get basic metrics, including filtered time_series/trajectories from metrics function
for j = 1:length(FeatureList)
    featName = FeatureList{j};
    if ~isfield(features, featName)
        switch featName
%get basic metrics, including filtered time_series/trajectories from metrics function
            case {'derivatives', 'integrals', 'time_series'}
                features.(featName)     = metrics.(featName);
            case{'intwin1','intwin3'}
                features.(featName)     = metrics.(featName);
            case{'baseline','responder_index','off_times', 'responders_fraction'}
                features.(featName) = metrics.(featName);
            case{'max_amplitude','max_integral','max_derivative','min_derivative'}
                features.(featName)     = metrics.(featName);
            case{'peakfreq'}
                features.(featName)     = metrics.(featName);
         %todo osc_frac from basic metrics function, pick threshold or replace entirely with osc_cat?
            case{'oscfrac'}
                features.(featName)     = metrics.(featName);
            case{'pk1_time','pk1_amp','pk2_time','pk2_amp', 'pk1_width', 'pk1_prom', 'pk1_height', 'pk2_width', 'pk2_prom', 'pk2_height',}  %% note these are calculated in nkbmetrics all now, not in peak_stats
                features.(featName)     = metrics.(featName);
         %todo proper picking of duration and envelope threshold and method,
            case{'envelope','duration'}
                features.(featName)     = metrics.(featName);
            case{'envelope_sigma','duration_sigma'}
                features.(featName) = metrics.(featName);
%% get metrics from Ade's computeFeature function
% metrics that are directly calculated
             case {'pk2_ratio'}
                features.(featName)     =(metrics.pk2_amp+metrics.baseline)./(metrics.pk1_amp+metrics.baseline); %% AS 6/19/2020 changed to add baseline
            case {'pk2_ratio_prom'}
                features.(featName)     =(metrics.pk2_prom)./(metrics.pk1_prom);
             case {'median_derivative'}
                features.(featName)     = nanmedian(metrics.derivatives(:,1:end),2);
             case {'mean_derivative'}
                features.(featName)     = nanmean(metrics.derivatives(:,1:end),2);           
% metrics that require additional function calls
            case{'fold_change'}
                features.(featName)     = fold_change_features.fold_change;
            case{'max_fold_change'}
                features.(featName)     = fold_change_features.max_fold_change;
            case{'max_value'}
                features.(featName)     = fold_change_features.max_value;
            case{'time2HalfMaxValue'}
                features.(featName)     = fold_change_features.time2HalfMaxValue;
            %time to half max integral within first 8 hours
            case {'time2HalfMaxIntegral'}
                endFrame = min(96, size(metrics.integrals,2));
                features.(featName)     = get_time_to_half_max_integral_new(metrics.integrals(:,1:endFrame), FramesPerHour);
            case{'max_pk1_speed'}
                features.(featName)     = get_max_pk1_speed_new(metrics.pk1_time, metrics.derivatives, FramesPerHour);
            case {'osc_cats'}
                feature_data            =  get_osc_cats_new(metrics.peakfreq,metrics.off_times,'cutoff_fq', 0.42);
                feature_data            = removecats(feature_data, {'off'});
                feature_data            = grp2idx(reordercats(feature_data, {'non_osc', 'osc'}));
                features_data           = feature_data -1;
                features.(featName)     = features_data + 0.1*rand(length(feature_data), 1); %%% adding a little noise so dist is not binary
            case{'integrals_pos','time2HalfMaxPosIntegral','max_pos_integral'}
                features.(featName)     = pos_integral_features.(featName);
            case{'last_falltime'}
                features.(featName)     = computeDuration_new(metrics.time_series, metrics.responder_index);
            case{'pos_pk1_amp', 'pos_pk1_time', 'pos_pk1_width', 'pos_pk1_prom', 'pos_pk1_height',...
        'pos_pk2_amp', 'pos_pk2_time', 'pos_pk2_width', 'pos_pk2_prom', 'pos_pk2_height'}
                features.(featName)     = pos_peak_features.(featName);
            case{'max_pos_pk1_speed'}
                features.(featName)     = get_max_pk1_speed_new(pos_peak_features.pos_pk1_time, metrics.derivatives, FramesPerHour);
            case {'pos_pk2_ratio'}
                features.(featName)     =(pos_peak_features.pos_pk2_amp)./(pos_peak_features.pos_pk1_amp);
            case {'pos_pk2_ratio_prom'}
                features.(featName)     =(pos_peak_features.pos_pk2_prom)./(pos_peak_features.pos_pk1_prom);
%PeakStats            
            case PeakStatsList
                    if contains(featName, 'ipt')
                        all_cells= nan(size(peak_stats.kept, 1),size(peak_stats.(featName),2));
                        all_cells(peak_stats.kept,:)= peak_stats.(featName);
                        features.(featName) = all_cells;
                    else
                        features.(featName)=peak_stats.(featName);
                    end

%SignalStats
            case SignalStatsList
                features.(featName)=sig_stats.(featName);  
                   
        end
    end
end


features.period_peak_times = diff(features.peak_times,1,2)*5;%mins
features.period_valley_times = diff(features.valley_times,1,2)*5;%mins

features.median_peak_period = median(features.period_peak_times,2,'omitnan');
features.mean_peak_period = mean(features.period_peak_times,2,'omitnan');
features.std_peak_period = std(features.period_peak_times,0,2,'omitnan');

features.median_valley_period = median(features.period_valley_times,2,'omitnan');
features.mean_valley_period = mean(features.period_valley_times,2,'omitnan');
features.std_valley_period = std(features.period_valley_times,0,2,'omitnan');
