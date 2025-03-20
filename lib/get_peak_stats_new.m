function output= get_peak_stats_new(time_series_no_base_ded,baseline,varargin)

%Calculates peak statistics
%--------------------------------------------------------------------------
p=inputParser;
addRequired(p,'time_series_no_base_ded');

addParameter (p,'min_pks',2,@isnumeric);
addParameter (p,'max_pk_diff',35,@isnumeric);
addParameter(p, 'FramesPerHour', 12, @isnumeric);
parse (p,time_series_no_base_ded,varargin{:});

min_pks =p.Results.min_pks;
max_pk_diff=p.Results.max_pk_diff;
FramesPerHour = p.Results.FramesPerHour;

[output.peak_times,output.peak_amps, output.valley_times, output.valley_amps]=nfkbpeaks_new(time_series_no_base_ded(:, 1:end), baseline, 'BeginFrame',3,'MinHeight',0.75,'MinDist',6);
ipt=diff(output.peak_times,1,2);

%%
ipt(ipt>max_pk_diff)=nan;

tot_pks = sum(~isnan(ipt),2)+1;
output.kept = true([size(ipt,1),1]);

ipt(tot_pks<min_pks,:) = [];
output.kept(tot_pks<min_pks)=false;
ipt = ipt.*60/FramesPerHour;%convert to minutes

%Brooks' analysis excludes ipt > 35, excludes cells with <5peaks

output.mean_ipt    =nanmean(ipt,2);
output.median_ipt  =nanmedian(ipt,2);
output.std_ipt     =nanstd(ipt,[],2);
% output.var_ipt        =nanvar(ipt,[],2);
output.max_ipt     = nanmax(ipt,[],2);
output.min_ipt     = nanmin(ipt,[],2);
output.cv_ipt      = output.std_ipt./output.mean_ipt;
output.ipt         =ipt;
output.num_peaks   =tot_pks;
%%
output.mean_peak_amp   = nanmean(output.peak_amps, 2); 
output.median_peak_amp = nanmedian(output.peak_amps, 2); 
% output.var_peak_amp  = nanvar(output.peak_amps, [],2); 
output.std_peak_amp    = nanstd(output.peak_amps, [],2); 
output.cv_peak_amp     = output.std_peak_amp./output.mean_peak_amp; 

output.peak2trough = output.peak_amps(:,1:end-1)-output.valley_amps;
minVals = zeros(size(time_series_no_base_ded,1),1);
time_series = time_series_no_base_ded - baseline;
for row=1:size(minVals,1)
    pk_frame =output.peak_times(row, 1);
    if isnan(pk_frame)
        continue;
    end
    minVals(row) = min(time_series(row, 1:pk_frame));
end
troughs = [minVals, output.valley_amps];

output.trough2peak     = troughs-output.peak_amps;

% output.peak2
output.max_peak2trough     =nanmax(output.peak2trough, [],2);
output.max_trough2peak    =nanmax(output.trough2peak, [],2);

output.mean_peak2trough    =nanmean(output.peak2trough, 2);
output.mean_trough2peak    =nanmean(output.trough2peak,2);

output.median_peak2trough  =nanmedian(output.peak2trough,2); 
output.median_trough2peak  =nanmedian(output.trough2peak,2);

output.min_peak2trough     =nanmin(output.peak2trough,[],2); 
output.min_trough2peak     =nanmin(output.trough2peak, [],2);

output.std_trough2peak    =nanstd(output.trough2peak,[],2);
output.std_peak2trough    =nanstd(output.peak2trough,[],2); 
output.cv_trough2peak     =output.std_trough2peak./output.mean_trough2peak;
output.cv_peak2trough     =output.std_peak2trough./output.mean_peak2trough; 

end