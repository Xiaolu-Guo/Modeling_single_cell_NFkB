function [nfkb, all_dims, sc_dims, names_1D, names_sc_1D] = loadnfkb_signaling_codons_AS_XG(nfkb_name)
%%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [nfkb, all_dims, sc_dims, names_1D] = loadnfkb()
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% LOADNFKB loads the master set of NFkB runs used for (most) BMDM analyses
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% If nfkb.mat doesn't exist, load, measure, and save the individual experiments
stimuli = {'TNF', 'PIC', 'LPS', 'P3K', 'CpG'};
%doses = {};

ids = cell(1, 5);
count = 1;
for i = 1:length(stimuli)
    ids{count} = stimuli{i}; %strcat(stimuli{i},doses{j});
    count = count + 1;
end

% Define name of nfkb file

if nargin < 1
    P = mfilename('fullpath');
P2 = mfilename;
    nfkb_name = [P(1:(length(P)-length(P2))) , 'all_stim', '_nfkb.mat'];
    
    nfkb_name = [P(1:(length(P)-length(P2))) , 'example_2bits_nfkb.mat'];
end
% nfkb_name = [P(1:(length(P)-length(P2))) , 'example_1bit_nfkb.mat'];
% nfkb_name = [P(1:(length(P)-length(P2))) , 'example_0bit_nfkb.mat']; % not working
% nfkb_name = [P(1:(length(P)-length(P2))) , 'example_0bit_v2_nfkb.mat']; % not working
% nfkb_name = [P(1:(length(P)-length(P2))) , 'example_s1bit_nfkb.mat']; %not working
% nfkb_name = [P(1:(length(P)-length(P2))) , 'example_s1bit_nfkb.mat'];%not working




% Combine processed data into a single structure
if ~exist(nfkb_name,'file')
    nfkb = struct;
    for i =1:length(ids)
        %nfkb(i).metrics = table2struct(readtable(['F:\Shared drives\Signaling Systems Lab\Apeksha\information theory code for Xiaolu\channel_capacity\data_for_example\', ids{i}, '_metrics.xlsx'], 'Sheet', 'all_metrics'), 'ToScalar', true);
        nfkb(i).sc_metrics = table2struct(readtable(['/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_Mutual_Information/data_for_example/', ids{i}, '_metrics.xlsx'], 'Sheet', 'signaling_codons'), 'ToScalar', true);
        nfkb(i).id = ids(i);
        nfkb(i).ids = ids;
    end
    save(nfkb_name, 'nfkb')
else
    load(nfkb_name);
end

%metric_names = fieldnames(nfkb(1).metrics);
%disp(['number of features = ', num2str(length(metric_names))])
%for i = 1:length(metric_names)
%    metric_name = metric_names{i};
%    pct_nan = zeros(1, length(ids));
%    for j =1:length(ids)
%        pct_nan(j) = sum(isnan(nfkb(j).metrics.(metric_name)))/length(nfkb(j).metrics.(metric_name));
%    end
%    if max(pct_nan)>0.2 %if more than 20% of feature values missing, remove that metric field
%        for j = 1:length(ids)
%            nfkb(j).metrics = rmfield(nfkb(j).metrics, metric_name);
%        end
%    end
%end
%metric_names = fieldnames(nfkb(1).metrics);
%disp(['number of features after removing missing values = ', num2str(length(metric_names))])

% If extra output arguments are defined, combine all metrics/measurements into one big matrix
% (record dimension names as well)
if nargout~=1
    %mfields = fieldnames(nfkb(1).metrics);
    sc_fields = fieldnames(nfkb(1).sc_metrics);
    names_1D = {};
    names_sc_1D = {};
    
    for i = 1:length(nfkb)
        % Downsample time series, integrals, and derivatives slightly
        % nfkb(i).metrics.time_series = nfkb(i).metrics.time_series(:,[1:36,37:2:180]);
        % nfkb(i).metrics.integrals = nfkb(i).metrics.integrals(:,1:2:240);
        % nfkb(i).metrics.derivatives = nfkb(i).metrics.derivatives(:,[1:36,37:2:180]);
        
        % Combine all dimensions into a single gigantic cell thing
        all_dims{i} = [];
        %for j = 1:length(mfields)
        %    all_dims{i} = cat(2,all_dims{i},nfkb(i).metrics.(mfields{j}));
        %    if i==1
        %        for k = 1:size(nfkb(i).metrics.(mfields{j}),2)
        %            names_1D = cat(1,names_1D,[mfields{j},'_',numseq(k,3)]);
        %        end
        %    end
        %end
        sc_dims{i} = [];
        for j = 1:length(sc_fields)
            sc_dims{i} = cat(2,sc_dims{i},nfkb(i).sc_metrics.(sc_fields{j}));
            if i==1
                for k = 1:size(nfkb(i).sc_metrics.(sc_fields{j}),2)
                    names_sc_1D = cat(1,names_sc_1D,[sc_fields{j},'_',numseq(k,3)]);
                end
            end
        end
    end
    
    % Reload full nfkb.mat
    load(nfkb_name)
end

% If called without output, assign everything into base workspace
if nargout<1
    assignin('base','nfkb',nfkb)
    assignin('base','all_dims',all_dims)
    assignin('base','names_1D',names_1D)
    assignin('base','sc_dims',sc_dims)
    assignin('base','names_sc_1D',names_sc_1D)
end

