tab_colors = [31, 119, 180;
    255, 127, 14;
    44, 160, 44;
    214, 39, 40;
    148, 103, 189;
    140, 86, 75
    ]/255;
%%
%create PCA plot using all features

path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["PIC"];%, "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
load(strcat(path, stimuli,"_nfkb.mat"));
%build matrix and labels list for PCA

polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];
labels = [];
for i=1:6
    labels = [labels; repmat(polarizations(i), length(nfkb(i).sc_metrics.Duration), 1)];
end

colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"];
scatter_colors = [];
for i=1:6
    scatter_colors = [scatter_colors; repmat(colors(i), length(nfkb(i).sc_metrics.Duration), 1)];
end

%data = zeros(length(labels), 6 + length(union_stim));

metrics_lst = fieldnames(nfkb(1).metrics);
data = zeros(length(labels), length(metrics_lst));

for i=1:length(metrics_lst)
    feat_name = metrics_lst{i};
    feat_data = [];
    for j=1:6
        feat_data = [feat_data; nfkb(j).metrics.(feat_name)];
    end
    data(:, i)=feat_data;
end

data = array2table(data,'VariableNames', metrics_lst); 

%signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};
%for i=1:length(signaling_codons )
%    feat_name = signaling_codons {i};
%    feat_data = [];
%    for j=1:6
%        feat_data = [feat_data; nfkb(j).sc_metrics.(feat_name)];
%    end
%    data(:, length(union_stim)+i)=feat_data;
%end

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

nan_feats = sum(isnan(table2array(data)))>0.5*length(metrics_lst);
low_var_feats = abs(coeff_var)<1;

data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(data)));


gscatter(score(:, 1), score(:, 2), labels, lines(6))

scatter(score(:, 1), score(:, 2), [], scatter_colors)

%trying UMAP
%(https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

%umap = run_umap(data);
%lbls_omit_nan = labels(sum(isnan(data), 2) == 0);

%%
%create PCA plot with axes shared for all stimuli

path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];

stim = [];
pol = [];
for s=1:length(stimuli)
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    %build matrix and labels list for PCA
    labels = [];
    for i=1:6
        labels = [labels; repmat(polarizations(i), length(nfkb(i).sc_metrics.Duration), 1)];
    end
    length(labels)
    pol = [pol; labels];
    stim = [stim; repmat(stimuli(s), length(labels), 1)];
end

metrics_lst = fieldnames(nfkb(1).metrics);
data = [];

for s=1:8
    stim_data = zeros(sum(stim==stimuli(s)), length(metrics_lst));
    size(stim_data)
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    for i=1:length(metrics_lst)
        feat_name = metrics_lst{i};
        feat_data = [];
        for j=1:6
            feat_data = [feat_data; nfkb(j).metrics.(feat_name)];
        end
        stim_data(:, i)=feat_data;
    end
    data = [data; stim_data];
end

data = array2table(data,'VariableNames', metrics_lst); 

%signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};
%for i=1:length(signaling_codons )
%    feat_name = signaling_codons {i};
%    feat_data = [];
%    for j=1:6
%        feat_data = [feat_data; nfkb(j).sc_metrics.(feat_name)];
%    end
%    data(:, length(union_stim)+i)=feat_data;
%end

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
low_var_feats = coeff_var==0;

pca_data = data;
pca_data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(pca_data)));

for s=1:8
    select_stim = stim==stimuli(s);
    subplot(2, 4, s)
    gscatter(score(select_stim, 1), score(select_stim, 2), pol(select_stim), lines(6))
    xlim([-50, 100])
    ylim([-40, 80])
    title(stimuli(s))
end

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 1)), 10);
pca_data.Properties.VariableNames{index}

%%
% PCA analysis using info theory results
%collect important features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
features = cell(0);
for i=1:length(stimuli)
    try
        load(strcat(path, stimuli(i), "\","cc_sc_plus_", stimuli(i), ".mat"));
        stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
    catch
        load(strcat(path, stimuli(i), "\","info.mat"));
    end
    [max_I, index] = max(info(4).I(:));
    I_size = size(info(4).I);
    feats_2_add = info(4).names(index);
    feats_2_add = strsplit(feats_2_add{1}, '+');
    features(end+1:end+length(feats_2_add)) = feats_2_add;
end
unique_feat = unique(features);

%collect labels and features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];

stim = [];
pol = [];
for s=1:length(stimuli)
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    %build matrix and labels list for PCA
    labels = [];
    for i=1:6
        labels = [labels; repmat(polarizations(i), length(nfkb(i).sc_metrics.Duration), 1)];
    end
    length(labels)
    pol = [pol; labels];
    stim = [stim; repmat(stimuli(s), length(labels), 1)];
end

data = [];
signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};

for s=1:8
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    
    stim_data = zeros(sum(stim==stimuli(s)), length(unique_feat)+length(signaling_codons));
    for i=1:length(unique_feat)
        feat_name = unique_feat{i}(1:end-4);
        feat_data = [];
        for j=1:6
            feat_data = [feat_data; nfkb(j).metrics.(feat_name)];
        end
        stim_data(:, i)=feat_data;
    end
    
    for i=1:length(signaling_codons )
        feat_name = signaling_codons{i};
        feat_data = [];
        for j=1:6
            feat_data = [feat_data; nfkb(j).sc_metrics.(feat_name)];
        end
        stim_data(:, length(unique_feat)+i)=feat_data;
    end
    data = [data; stim_data];
end

data = array2table(data,'VariableNames', [unique_feat, signaling_codons]);

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

%nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
%low_var_feats = coeff_var==0;

pca_data = data;
%pca_data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(pca_data)));

for s=1:8
    select_stim = stim==stimuli(s);
    subplot(2, 4, s)
    gscatter(score(select_stim, 1), score(select_stim, 2), pol(select_stim), lines(6))
    xlim([-10, 15])
    ylim([-20, 10])
    title(stimuli(s))
end

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 2)), 10);
pca_data.Properties.VariableNames{index}

%%
%using just signaling codons

%collect labels and features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "R84", "PIC", "P3K", "FLA", "CpG", "FSL", "LPS"];
polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];

stim = [];
pol = [];
for s=1:length(stimuli)
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    %build matrix and labels list for PCA
    labels = [];
    for i=1:6
        labels = [labels; repmat(polarizations(i), length(nfkb(i).sc_metrics.Duration), 1)];
    end
    length(labels)
    pol = [pol; labels];
    stim = [stim; repmat(stimuli(s), length(labels), 1)];
end

data = [];
signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};

for s=1:8
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    
    stim_data = zeros(sum(stim==stimuli(s)), length(signaling_codons));
    
    for i=1:length(signaling_codons)
        feat_name = signaling_codons{i};
        feat_data = [];
        for j=1:6
            feat_data = [feat_data; nfkb(j).sc_metrics.(feat_name)];
        end
        stim_data(:, i)=feat_data;
    end
    data = [data; stim_data];
end

data = array2table(data,'VariableNames', [signaling_codons]);

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

%nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
%low_var_feats = coeff_var==0;

pca_data = data;
%pca_data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(pca_data)));

for s=1:8
    select_stim = stim==stimuli(s);
    subplot(2, 4, s)
    gscatter(score(select_stim, 1), score(select_stim, 2), pol(select_stim), tab_colors, [], [] , 'off')
    xlim([-5, 7])
    ylim([-5, 7])
    title(stimuli(s))
end

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 1)), 10);
pca_data.Properties.VariableNames{index}

b = bar(coeff(:, (1:2)))
b(1).FaceColor=[87, 86, 86]/255
b(2).FaceColor= [150, 149, 149]/255
%%
%mean/median responses
%collect labels and features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];

data = zeros(6, 48);
range_data = zeros(12, 48);
signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};

for s=1:8
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    
    for j = 1:6
        pol_data = zeros(1, length(signaling_codons));
        pol_range_data = zeros(2, length(signaling_codons));
        for i=1:length(signaling_codons)
            feat_name = signaling_codons{i};
            pol_data(i)=nanmean(nfkb(j).sc_metrics.(feat_name));
            pol_range_data(1, i) = prctile(nfkb(j).sc_metrics.(feat_name), 25);
            pol_range_data(2, i) = prctile(nfkb(j).sc_metrics.(feat_name), 75);
        end
        data(j, (s-1)*6+1:s*6)=pol_data;
        range_data((j-1)*2+1:j*2, (s-1)*6+1:s*6)=pol_range_data;
    end
end

var_names = [];
for s=1:8
    var_names = [var_names, strcat(stimuli(s), '_',signaling_codons)];
end

data = array2table(data,'VariableNames', var_names);

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

%nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
%low_var_feats = coeff_var==0;

pca_data = data;
%pca_data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(data)));

gscatter(score(:, 1), score(:, 2), polarizations, tab_colors, [], [300] , 'off')
ylim([-8, 8])
xlim([-8, 8])
%plot 25 and 75% percentiles
proj_range = range_data * coeff;
gscatter(proj_range(:, 1), proj_range(:, 2), repelem(polarizations, 2), tab_colors, [], [], 'off')

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 1)), 5);
var_names{index}
table(var_names(index)', coeff(index, 2))

bar(coeff(:, (1:2)))
set(gca, 'XTickLabel', var_names) 

%%
%add theoretical values?
%collect labels and features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
polarizations = ["M0", "IFNb", "IFNg", "Il10", "Il13", "Il4"];

data = zeros(6, 48);
range_data = zeros(12, 48);
signaling_codons = {'Duration', 'EarlyVsLate', 'OscVsNon', 'PeakAmplitude', 'Speed', 'TotalActivity'};

for s=1:8
    load(strcat(path, stimuli(s),"_nfkb.mat"));
    
    for j = 1:6
        pol_data = zeros(1, length(signaling_codons));
        pol_range_data = zeros(2, length(signaling_codons));
        for i=1:length(signaling_codons)
            feat_name = signaling_codons{i};
            pol_data(i)=nanmedian(nfkb(j).sc_metrics.(feat_name));
            pol_range_data(1, i) = prctile(nfkb(j).sc_metrics.(feat_name), 25);
            pol_range_data(2, i) = prctile(nfkb(j).sc_metrics.(feat_name), 75);
        end
        data(j, (s-1)*6+1:s*6)=pol_data;
        range_data((j-1)*2+1:j*2, (s-1)*6+1:s*6)=pol_range_data;
    end
end

var_names = [];
for s=1:8
    var_names = [var_names, strcat(stimuli(s), '_',signaling_codons)];
end

data = array2table(data,'VariableNames', var_names);

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);
histogram(abs(coeff_var))

%remove variables with mostly nan values (or low variance?)

%nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
%low_var_feats = coeff_var==0;

pca_data = data;
%pca_data(:, nan_feats|low_var_feats)=[];

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(data)));

gscatter(score(:, 1), score(:, 2), polarizations, tab_colors, [], [] , 'off')

%plot 25 and 75% percentiles
proj_range = range_data * coeff;
gscatter(proj_range(:, 1), proj_range(:, 2), repelem(polarizations, 2), tab_colors, [], [], 'off')

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 2)), 5);
var_names{index}
table(var_names(index)', coeff(index, 2))

bar(coeff(:, (1:2)))
set(gca, 'XTickLabel', var_names) 