function [output] = get_positive_peaks(time_series_no_base_ded, baselines, FramesPerHour, MinLifetime)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% base off of global_peaks_new and code in nfkbmetrics_new 6/22/2020
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GLOBALPEAKS seeks to find the "dominant" peaks in an input vector - output will be sorted
% from most-to-least dominant.
%
% INPUT:
% vect            input vector--needs to be not baseline subtracted so all
% input is positive (can more easily get width)
% num_peaks       number of desired "dominant" peaks
%
% OUTPUT:
% peaks          peak values - are outputed as baseline subtracted
% locs           peak locations (index)
% heights        peak height (above nearest troughs (defined by nearest dominant peak neigbor only))
% HMW            peak width
% prom           peak prominence (similar to height, except now includes
% height above nearest troughs where all neighboring peaks are shorter)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% METRICS OF Positive AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude/prominence/width/height
pk_feats = {'pos_pk1_amp', 'pos_pk1_time', 'pos_pk1_width', 'pos_pk1_prom', 'pos_pk1_height',...
        'pos_pk2_amp', 'pos_pk2_time', 'pos_pk2_width', 'pos_pk2_prom', 'pos_pk2_height'};
for j=1:length(pk_feats)
    output.(pk_feats{j}) = nan(size(time_series_no_base_ded,1),1);
end
for j = 1:size(output.pos_pk1_time,1)    
    [pks, locs, width, prom, heights] = global_pos_peaks(time_series_no_base_ded(j,1:min([90,MinLifetime])),baselines(j),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order); width = width(order); prom = prom(order); heights = heights(order);
    while min(diff(locs))<6
        tmplst = find(diff(locs)==min(diff(locs)),1,'first');
        tmplst = tmplst + (pks(tmplst)>=pks(tmplst+1));
        pks(tmplst) = []; locs(tmplst) = []; width(tmplst) = []; prom(tmplst) = []; heights(tmplst) = [];
    end
    pks(locs< 3) = [];  %SL has + 1, but here our delay is smaller
    locs(locs< 3) = [];
    width(locs< 3) = [];
    prom(locs< 3) = [];
    heights(locs< 3) = [];
    if ~isempty(locs)
        output.pos_pk1_time(j) = locs(1);
        output.pos_pk1_amp(j) = pks(1);
        output.pos_pk1_width(j) = width(1);
        output.pos_pk1_prom(j) = prom(1);
        output.pos_pk1_height(j) = heights(1);
    end
    if length(locs)>1
        output.pos_pk2_time(j) = locs(2);
        output.pos_pk2_amp(j) = pks(2);
        output.pos_pk2_width(j) = width(2);
        output.pos_pk2_prom(j) = prom(2);
        output.pos_pk2_height(j) = heights(2);
    end
end
    output.pos_pk1_time = (output.pos_pk1_time-1)/FramesPerHour;
    output.pos_pk2_time = (output.pos_pk2_time-1)/FramesPerHour;

function [peaks, locs, HMW, prom, heights] = global_pos_peaks(vect, baseline, num_peaks)
% Find all peaks in vector; 1st global peak is maximum value overall
[all_peaks, all_locs, all_HMW, all_prom] = findpeaks(vect, 'WidthReference', 'halfheight');

all_peaks(((all_locs==1)) | (all_locs==length(vect))) = [];
all_HMW(((all_locs==1)) | (all_locs==length(vect))) = [];
all_prom(((all_locs==1)) | (all_locs==length(vect))) = [];
all_locs(((all_locs==1)) | (all_locs==length(vect))) = [];

getclosest = @(idx,vect) vect(find(abs(vect-idx)==min(abs(vect-idx)),1,'first'));
peaks = [];
locs  = [];
HMW = [];
prom = [];
heights = [];

while length(peaks) < num_peaks
    % Eliminate peaks that have been found already
    tmp_peaks = all_peaks;
    tmp_locs = all_locs;
    tmp_HMW = all_HMW;
    tmp_prom = all_prom;
    tmp_peaks(ismember(tmp_locs,locs)) = [];
    tmp_locs(ismember(tmp_locs,locs)) = [];
    tmp_HMW(ismember(tmp_locs,locs)) = [];
    tmp_prom(ismember(tmp_locs,locs)) = [];
    if isempty(tmp_peaks)
      break
    end
    
    % For each candidate, identify nearest peaks - maximize difference btw candidate and two nearest troughs.  
    diffs = zeros(size(tmp_peaks));
    loc_compare = [1 locs length(vect)];
    for i = 1:length(tmp_locs)
        tmp = loc_compare; tmp(tmp>=tmp_locs(i)) = inf;
        trough1 = min(vect(getclosest(tmp_locs(i),tmp):tmp_locs(i)));
        tmp = loc_compare; tmp(tmp<=tmp_locs(i)) = inf;
        trough2 = min(vect(tmp_locs(i):getclosest(tmp_locs(i),tmp)));
        diffs(i) = tmp_peaks(i) - max([trough1, trough2]);
    end
    
    all_HMW(tmp_peaks-baseline<0) = [];
    all_prom(tmp_peaks-baseline<0) = [];
    all_locs(tmp_peaks-baseline<0) = [];
    diffs(tmp_peaks-baseline<0) = [];
    all_peaks(tmp_peaks-baseline<0) = [];
    
    peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
    locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
    HMW = [HMW, tmp_HMW(find(diffs==max(diffs),1,'first'))];
    prom = [prom, tmp_prom(find(diffs==max(diffs),1,'first'))];
    heights = [heights, max(diffs)];
end
peaks = peaks - baseline;
end
end