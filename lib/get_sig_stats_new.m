function sig_stats = get_sig_stats_new(time_series, varargin)
% Calculates signal statistics for NFkB trajectories
%-----------------------------------------------------------------------
% INPUTS
%   REQUIRED
%       Data:         nxt array 
%   OPTIONAL
%   'FqRange',  frequency range of interest (2-element vector)
%   'Fs'         sample rate
%   'FillMethod' method for filling missing Data 
%                     {'previous',next', 'linear', 'spline', 'pchip'}
%   'SmoothMethod' {'lowess', 'sgolay', 'none', 'wavelet'}
% OUTPUT
%   stats:    struct of wavelet, and spectral statistics

valid_FillMethod  = {'previous', 'next','linear', 'spline', 'pchip'};
isvalidfill =@(x)ismember (x, valid_FillMethod);
p=inputParser;
addRequired(p, 'time_series',@isnumeric);
addParameter(p, 'FqRange',[0.33 1] , @isnumeric);
addParameter(p, 'Fs',12, @isnumeric);
addParameter(p, 'FillMethod','linear', isvalidfill);
addParameter(p, 'SmoothMethod',"sgolay", @istext);  %% can also be lowess --edited by AS on 6/17/2020

parse(p,time_series, varargin{:});

sig_stats=struct;

Data = time_series(:, 1:end);
Data=fillmissing(Data,p.Results.FillMethod, 2, 'EndValues','extrap');
smoothData = zeros(size(Data));
for i = 1:size(smoothData, 1)
    smoothData(i, :) = smooth(Data(i, :), p.Results.SmoothMethod);
end

Fs=p.Results.Fs; freq_range =p.Results.FqRange;

sig_stats.powerbw = powerbw(smoothData',Fs, freq_range)';
sig_stats.medfreq      = medfreq(smoothData', Fs,freq_range)';
sig_stats.meanfreq     = meanfreq(smoothData', Fs, freq_range)'; 
sig_stats.peak2rms     =peak2rms(smoothData,2);     
sig_stats.rms          =rms(smoothData')';
sig_stats.peak2peak    = peak2peak(smoothData,2);
sig_stats.mean_movmad  =mean( movmad(smoothData,3,2),2);
sig_stats.mean_movstd       =mean( movstd(smoothData,3,[],2),2);
sig_stats.mean_movvar       =mean( movvar(smoothData,3,[],2),2);

%psd and power (scales PSD by equiv noise bandwidth of window)
[psd,fq]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
sig_stats.fq=fq';
sig_stats.psd=transpose(psd./sum(psd,1));

[pwr,~]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','power');
%sig_stats.fq_2=fq';
sig_stats.power=transpose(pwr./sum(pwr,1));

%oscpower aka bandpower
psd = transpose(sig_stats.psd) ; fq = transpose(sig_stats.fq);% 
bp = bandpower(psd,fq,freq_range, 'psd')';
sig_stats.oscpower =bp;

%oscillation frequency--find peaks within the frequency range
    bandfilter= @(x) x<= max(freq_range) & x>= min(freq_range);normalize =@(x) x/sum(x);
    ix =bandfilter(sig_stats.fq);
    peakFun =@(a) arrayfun(@(j) findpeaks(a(j,:), sig_stats.fq(ix),...
                    'SortStr', 'descend', 'MinPeakProminence', 0.0055), 1:size(a,1), 'UniformOutput',false);
    [peaks,locs] = peakFun(sig_stats.psd(:,ix)) ; %peaks = psd, locs = frequency
    freq =zeros(size(peaks)); 
        for j = 1:numel(peaks)
            if numel(peaks{j}) > 1
            % more than one peak within the range, take weighted
            % sum of frequency 
                wgts = normalize(peaks{j}); 
                freq(j) = sum(locs{j}.*wgts); 
            elseif ~isempty(peaks{j})
                freq(j) = locs{j}; 
            end
        end
sig_stats.oscfreq = freq';
%oscillatory bandwidth       
sig_stats.oscbandwidth     =obw(smoothData',Fs)';
%max entropy
max_entropy= zeros(size(smoothData,1),1);    
time_pts = max_entropy;
for j =1:numel(max_entropy)
  [sig_entropy, tp]= pentropy(smoothData(j,:), Fs/3600,'FrequencyLimit',freq_range./3600);   
   [max_entropy(j), ix] = max(sig_entropy);
   time_pts(j) = tp(ix)/3600;
end
sig_stats.max_entropy = max_entropy; 
sig_stats.noise_est =  wnoisest(Data)';
end
