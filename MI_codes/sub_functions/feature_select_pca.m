%%
%using F-score from ANOVA for feature selection

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
    pol = [pol; labels];
    stim = [stim; repmat(stimuli(s), length(labels), 1)];
end

metrics_lst = fieldnames(nfkb(1).metrics);
data = [];

for s=1:8
    stim_data = zeros(sum(stim==stimuli(s)), length(metrics_lst));
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

feat_vars = std(table2array(data), 'omitnan');
feat_means = mean(table2array(data), 'omitnan');
coeff_var = feat_vars;
coeff_var(feat_vars>0) = feat_vars(feat_vars>0)./feat_means(feat_vars>0);

%remove variables with mostly nan values (or low variance?)

nan_feats = sum(isnan(table2array(data)))>0.5*length(stim);
low_var_feats = coeff_var==0;
data(:, nan_feats|low_var_feats)=[];

F_scores = zeros(8, size(data, 2));

for m=1:size(data, 2)
    for s=1:8
        [~, tbl] = anova1(table2array(data(stim==stimuli(s), m)), pol(stim==stimuli(s)), 'off');
        F_scores(s, m) = tbl{2, 5};
    end
end

%now could take mean of F-scores across stimuli to select features or just
%take top values

[max_F, index_F] = maxk(mean(F_scores), 100);

[max_F, index_F] = sort(F_scores(:), 'descend');
index_tbl = repmat(1:length(metrics_lst), 8, 1);
index_tbl = index_tbl(:);
index_tbl = index_tbl(index_F);
index_tbl = unique(index_tbl, 'stable');

pca_data = data;
pca_data = pca_data(:, index_tbl(1:150));

nanzscore= @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

[coeff,score,latent,tsquared,explained,mu] = pca(nanzscore(table2array(pca_data)));

for s=1:8
    select_stim = stim==stimuli(s);
    subplot(2, 4, s)
    gscatter(score(select_stim, 1), score(select_stim, 2), pol(select_stim), lines(6), [], [], 'off')
    %gscatter(score(select_stim, 1), score(select_stim, 2), pol(select_stim), lines(6), [], repmat(1,6,1), 'off')
    %xlim([-20, 55])
    %ylim([-10, 20])
    title(stimuli(s))
end

%determine max loadings
[coeff_vals, index] = maxk(abs(coeff(:, 2)), 10);
pca_data.Properties.VariableNames{index}