%Information theory exploring output from veqC

path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\M';
polarization = ["", "ib", "ig", "i0", "i3", "i4"];
I = zeros(5, 6);
feat = cell(4, 6);
for i=1:length(polarization)
    load(strcat(path, polarization(i), "\","sc_info.mat"));
    I(1, i) = sc_info.I;
    
    load(strcat(path, polarization(i), "\","info.mat"));
    [I(2, i), index] = nanmax(info(1).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(1).I), index);
    feat(1, i) = info(1).names(r, c);
    [I(3, i), index] = nanmax(info(2).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(2).I), index);
    feat(2, i) = info(2).names(r, c);
    [I(4, i), index] = nanmax(info(3).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(3).I), index);
    feat(3, i) = info(3).names(r, c);
    [I(5, i), index] = nanmax(info(4).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(4).I), index);
    feat(4, i) = info(4).names(r, c);
    
    p = plot([0:4],I(:, i), '-o', 'DisplayName',polarization(i));
    hold on
    %ylim([0.4,1.1])
end
xlabel('dimension')
ylabel('max information')

path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
colors = ['#1b9e77', '#d95f02', "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"];
I = zeros(5, 8);
feat = cell(4, 8);
for i=1:length(stimuli)
    load(strcat(path, stimuli(i), "\","sc_info.mat"));
    I(1, i) = sc_info.I;
    
    try
        load(strcat(path, stimuli(i), "\","cc_sc_plus_", stimuli(i), ".mat"));
        stimuli = ["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
    catch
        load(strcat(path, stimuli(i), "\","info.mat"));
    end
    
    [I(2, i), index] = nanmax(info(1).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(1).I), index);
    feat(1, i) = info(1).names(r, c);
    [I(3, i), index] = nanmax(info(2).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(2).I), index);
    feat(2, i) = info(2).names(r, c);
    [I(4, i), index] = nanmax(info(3).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(3).I), index);
    feat(3, i) = info(3).names(r, c);
    [I(5, i), index] = nanmax(info(4).I, [], 'all', 'linear');
    [r, c] = ind2sub(size(info(4).I), index);
    feat(4, i) = info(4).names(r, c);
    
    p = plot([0:4],I(:, i), '-o', 'DisplayName',stimuli(i), 'Color', colors(i));
    hold on
    %ylim([0.4,1.1])
end
xlabel('dimension')
ylabel('max information')

%collect important features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\M';
polarization = ["", "ib", "ig", "i0", "i3", "i4"];
feat_by_pol = cell(1, 6);
union_pol = cell(0, 0);
for i=1:length(polarization)
    load(strcat(path, polarization(i), "\","info.mat"));
    [~, index] = sort(info(3).I(:), 'descend');
    I_size = size(info(3).I);
    features = info(3).names(index(1:round(0.01*I_size(2))));
    feat_lst = cell(0, 0);
    for j=1:length(features)
        feats_2_add = strsplit(features{j}, '+');
        feat_lst(end+1:end+length(feats_2_add)) = feats_2_add;
    end
    feat_by_pol{i} = unique(feat_lst);
    union_pol(end+1:end+length(feat_by_pol{i})) = feat_by_pol{i};
end
union_pol = unique(union_pol);

%collect important features
path = 'C:\Users\apeks\Box\Hoffmann\Spring 2021\information_theory-master\channel_capacity\';
stimuli = ["P3K"]; %["TNF", "PIC", "R84", "P3K", "FLA", "CpG", "FSL", "LPS"];
feat_by_stim = cell(1, 8);
union_stim = cell(0, 0);
for i=1:length(stimuli)
    load(strcat(path, stimuli(i), "\","info.mat"));
    [~, index] = sort(info(2).I(:), 'descend');
    I_size = size(info(2).I);
    features = info(2).names(index(1:round(0.01*I_size(2))));
    feat_lst = cell(0, 0);
    for j=1:length(features)
        feats_2_add = strsplit(features{j}, '+');
        feat_lst(end+1:end+length(feats_2_add)) = feats_2_add;
    end
    feat_by_stim{i} = unique(feat_lst);
    union_stim(end+1:end+length(feat_by_stim{i})) = feat_by_stim{i};
end
union_stim = unique(union_stim);

i = 1;
stimuli="LPS";
load(strcat(path, stimuli, "\","info.mat"));
[~, index] = sort(info(i).I(:), 'descend');
info_values = info(i).I(index);
feat_names = info(i).names(index);