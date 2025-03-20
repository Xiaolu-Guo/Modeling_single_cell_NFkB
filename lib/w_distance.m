function w_dis = w_distance(u_samples, v_samples, p)
% W_DISTANCE 1- and 2- Wasserstein distance between two probability
% measures estimated by matlab ksdensity() runction
%
%   wsd = WS_DISTANCE(u_samples, v_samples) returns the 1-Wasserstein
%   distance between the probability measures u and v
%   corresponding to the sample vectors u_samples and v_samples
%
%   wsd = WS_DISTANCE(u_samples, v_samples, p) returns the p-Wasserstein
%   distance between the probability measures u and v
%   corresponding to the sample vectors u_samples and v_samples.
%   p must be 1 or 2.
%

if ~exist('p', 'var')
    p = 1;
end

pts = linspace(min(min(u_samples),min(v_samples)),max(max(u_samples),max(v_samples)),101);

[~,~,bw_u_samples] = ksdensity(u_samples,pts,...
    'Function','pdf');
[~,~,bw_v_samples] = ksdensity(v_samples,pts,...
    'Function','pdf');%
bw = min(bw_u_samples,bw_v_samples);


if p == 1
    
    pts = linspace(min(min(u_samples),min(v_samples)),max(max(u_samples),max(v_samples)),101);
    [u_cdf,~] = ksdensity(u_samples,pts,...
        'Function','cdf','Bandwidth',bw);%,'Bandwidth',bw
    [v_cdf,~] = ksdensity(v_samples,pts,...
        'Function','cdf','Bandwidth',bw);%,'Bandwidth',bw
    
    w_dis = sum(abs(u_cdf - v_cdf) .* diff(pts));
    
elseif p == 2
    
    
    pts = linspace(0.01,0.99,99);
    [u_icdf,~] = ksdensity(u_samples,pts,...
        'Function','icdf','Bandwidth',bw);%,'Bandwidth',bw
    [v_icdf,~] = ksdensity(v_samples,pts,...
        'Function','icdf','Bandwidth',bw);%,'Bandwidth',bw
    
    pts_0 = [0,pts];
    w_dis = sqrt(sum((u_icdf-v_icdf).^2 .* diff(pts_0)));
    
else
    
    error('Only p=1 or p=2 allowed.')
    
end
end
