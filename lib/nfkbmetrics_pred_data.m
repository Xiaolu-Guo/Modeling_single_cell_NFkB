function [metrics,aux, graph] = nfkbmetrics_pred_data(pred_data,pred_parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% metrics = nfkbmetrics(id)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMETRICS_NEW uses the filter_nfkb function to filter and preprocess NFkB trajectories,
% then calculates related metrics regarding activation. Metric types include:
% 
% 1) time series (base NFkB dynamics, resampled to 12 frames/hr)
% 2) integrated activity
% 3) differentiated activity
% 4) calculated metrics: measuring aspects of oscillation, duration, timing ,and amplitude
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'          'on' or 'off' - show verbose output
% 'MinLifetime'      final frame used to filter for long-lived cells (default = 100)
% 'TrimFrame'        trim sets to common length (default = 254 timepoints) 
% 'ConvectionShift'  max allowed time shift between scenes (to correct for poor mixing - default is no shift allowed)
%
% OUTPUT: 
% metrics   structure with output fields
% aux       Extra data (e.g. fourier information (FFT, power, frequencies), thresholds used in envelope/duration)
% graph     main structure output from see_nfkb_native
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% INPUT PARSING
% Create input parser object, add required params from function input
            % p = inputParser;
            % % Required: ID input
            % valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
            %     'ID input must be spreadsheet ID or full file path');
            % addRequired(p,'id',valid_id);
            % % Optional parameters
            % addParameter(p,'MinLifetime',97, @isnumeric); %old defualt 130 -AS 6/18/2020
            % addParameter(p,'Verbose','off');
            % addParameter(p, 'GraphLimits', [-0.25, 4]);
            % addParameter(p,'TrimFrame',100, @isnumeric);
            % addParameter(p,'Delay',3);
            % addParameter(p,'onThresh', 3); %% changed from foldThresh/inputThresh; old default 1.25 -AS 6/18/2020
            % addParameter(p,'Frames',5); %old default 4-AS 6/18/2020
            % addParameter(p,'offPad',12) % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)
            % valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
            %     'Convection correction parameter must be single integer >= 0');
            % addParameter(p,'ConvectionShift',1, valid_conv);
            % parse(p,id, varargin{:})
            % 
            % OnThresh =p.Results.onThresh; %value set as the standard deviation theshold above baseline reqd to be considered a "responder" (also called sigmaThresh)
            % thresholdFrames = p.Results.Frames; 
            % %% INITIALIZATION. Load and process data. Interpolate time series, calculate deriv/integral approximations
            % 
            % MinLifetime = p.Results.MinLifetime; 
            % endFrame = p.Results.TrimFrame;
            % Delay = p.Results.Delay;
            % 
            % [graph, info, measure] = filter_nfkb(id,'MinLifetime',MinLifetime,'ConvectionShift',p.Results.ConvectionShift,...
            %     'Verbose', p.Results.Verbose, 'GraphLimits', p.Results.GraphLimits, 'Delay',Delay, 'TrimFrame', endFrame);
Delay=pred_parameters.Delay;
graph.var = pred_data;
graph.var_nfkb_no_base_ded = pred_data-pred_parameters.baseline;
graph.t = (0:5:5*(size(pred_data,2)-1))/60; % 1:time_points or 5*[1:time_points]
metrics.baseline = pred_parameters.baseline*ones(size(graph.var,1),1);%****************************************************************************************************??
FramesPerHour = pred_parameters.FramesPerHour; %should be 12?
%%
% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
t = min(graph.t):1/12:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series = nan(size(graph.var,1),length(t));
    for i = 1:size(graph.var,1)
        metrics.time_series(i,:) = interp1(graph.t,graph.var(i,:),t);
    end
    metrics.time_series_no_base_ded = nan(size(graph.var_nfkb_no_base_ded,1),length(t));%var_nfkb_no_base_ded****************************************************************************************************??
    for i = 1:size(graph.var_nfkb_no_base_ded,1)
        metrics.time_series_no_base_ded(i,:) = interp1(graph.t,graph.var_nfkb_no_base_ded(i,:),t);
    end
else
    metrics.time_series = graph.var;
    metrics.time_series_no_base_ded = graph.var_nfkb_no_base_ded;
end
baseline_stdv = nanstd(metrics.time_series(:,1:Delay),0,2);
metrics.time_series = metrics.time_series(:, Delay:end);
metrics.time_series_no_base_ded = metrics.time_series_no_base_ded(:, Delay:end);
t = t(Delay:end);

% 2) integrated activity
metrics.integrals = nan(size(metrics.time_series));
nan_removed = metrics.time_series;
nan_removed(isnan(nan_removed)) = 0;
for i = 1:size(metrics.integrals,1)
    metrics.integrals(i,:) = cumtrapz(t,nan_removed(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = medfilt1(metrics.time_series,3,[],2);
metrics.derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);

%% TRIM EVERYBODY to a common length (of "good" sets)
endFrame = size(pred_data,2);
try
    metrics.time_series = metrics.time_series(:,1:endFrame);
    metrics.integrals = metrics.integrals(:,1:endFrame);
    metrics.derivatives = metrics.derivatives(:,1:(endFrame-2));
    smoothed = smoothed(:,1:endFrame);
    t = t(1:endFrame);
catch
    disp(['Note: vectors too short to cap @ ',num2str(endFrame),' frames'])
end
%% MISC METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows (0-3, 1-4, etc) of activity
max_hr = floor(max(t));
metrics.intwin1 = nan(size(metrics.time_series,1),max_hr);
metrics.intwin3 = nan(size(metrics.time_series,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i);
    metrics.intwin1(:,i) = trapz(t(win),metrics.time_series(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2);
        metrics.intwin3(:,i) = trapz(t(win),metrics.time_series(:,win),2);
    end
end

% 5) amplitude/peak/on-vs-off metrics
% MAX/MIN metrics
metrics.max_amplitude = nanmax(metrics.time_series(:,1:end),[],2);
metrics.max_integral = nanmax(metrics.integrals(:,1:end),[],2);
metrics.max_derivative = nanmax(metrics.derivatives(:,1:end),[],2);
metrics.min_derivative = nanmin(metrics.derivatives(:,1:end),[],2);

% 6) ACTIVITY Determine which cells are responders and compute off-times
Wliml = 1; %first/lower time point of window to check for activity
Wlimu = 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = pred_parameters.thresholdFrames; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
smoothed_by_sigma = smoothed./baseline_stdv;
OnThresh=pred_parameters.OnThresh;
[metrics.responder_index, metrics.responders_fraction, metrics.off_times] = get_activity_metrics(smoothed_by_sigma, Wliml, Wlimu, OnThresh, blockLengthThresh);
metrics.off_times = metrics.off_times/FramesPerHour;
metrics.off_times(metrics.off_times<0) = 0;

%% METRICS OF OSCILLATION
% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series,1)
    if(metrics.off_times(i)>0)
        y = metrics.time_series(i,1:(depth));
        off_frame = min([length(y), metrics.off_times(i)*FramesPerHour+1+pred_parameters.offPad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq(i) = 3600*aux.freq(1);
    end
end

% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series,1)
        metrics.oscfrac(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac(i,j))
            metrics.oscfrac(i,j) = 0;
        end
    end
end

%% METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude/prominence/width/height
pk_feats = {'pk1_amp', 'pk1_time', 'pk1_width', 'pk1_prom', 'pk1_height',...
        'pk2_amp', 'pk2_time', 'pk2_width', 'pk2_prom', 'pk2_height'};
for i=1:length(pk_feats)
    metrics.(pk_feats{i}) = nan(size(metrics.time_series,1),1);
end
for i = 1:size(metrics.pk1_time,1)    
    [pks, locs, width, prom, heights] = globalpeaks_new(metrics.time_series_no_base_ded(i,1:min([90,pred_parameters.MinLifetime])),metrics.baseline(i),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order); width = width(order); prom = prom(order); heights = heights(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = []; locs(tmp) = []; width(tmp) = []; prom(tmp) = []; heights(tmp) = [];
    end
    pks(locs< 3) = [];  %SL has + 1, but here our delay is smaller
    locs(locs< 3) = [];
    width(locs< 3) = [];
    prom(locs< 3) = [];
    heights(locs< 3) = [];
    if ~isempty(locs)
        metrics.pk1_time(i) = locs(1);
        metrics.pk1_amp(i) = pks(1);
        metrics.pk1_width(i) = width(1);
        metrics.pk1_prom(i) = prom(1);
        metrics.pk1_height(i) = heights(1);
    end
    if length(locs)>1
        metrics.pk2_time(i) = locs(2);
        metrics.pk2_amp(i) = pks(2);
        metrics.pk2_width(i) = width(2);
        metrics.pk2_prom(i) = prom(2);
        metrics.pk2_height(i) = heights(2);
    end
end
    metrics.pk1_time = (metrics.pk1_time-1)/FramesPerHour;
    metrics.pk2_time = (metrics.pk2_time-1)/FramesPerHour;

%% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = medfilt1(metrics.time_series,5,[],2);
lowerThresh = 0;
upperThresh = 0.4;%7.1
aux.thresholds = linspace(lowerThresh,upperThresh,9); %%% need to think about these cutoffs
metrics.envelope = zeros(size(metrics.time_series,1),9);
for j = 1:length(aux.thresholds)
    thresholded = smoothed2(:, 1:end)>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope(i,j)
                    metrics.envelope(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope = metrics.envelope/FramesPerHour;

% Number of frames above a given threshold
metrics.duration = zeros(size(metrics.time_series,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration(:,i) = nansum(smoothed(:, 1:end)>aux.thresholds(i),2)/FramesPerHour;
end

%% NFkB METRICS OF DURATION Using sigma of baseline
%
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs after stimulation)
smoothed2_by_sigma= smoothed2./baseline_stdv;
upperThresh = 41.7;
aux.thresholds = linspace(0, upperThresh, 25);
metrics.envelope_sigma = zeros(size(metrics.time_series,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2_by_sigma(:,1:end)>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope_sigma(i,j)
                    metrics.envelope_sigma(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope_sigma = metrics.envelope_sigma/FramesPerHour;


% Number of frames above a given threshold
metrics.duration_sigma = zeros(size(metrics.time_series,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration_sigma(:,i) = nansum(smoothed_by_sigma(:,1:end)>aux.thresholds(i),2)/FramesPerHour;
end
end