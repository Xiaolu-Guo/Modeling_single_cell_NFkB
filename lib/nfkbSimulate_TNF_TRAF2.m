function [t,x, simdata, x_end] = nfkbSimulate_TNF_TRAF2(stimulus, names, params_mod, init_mod, options)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% nfkbSimulate(stimulus, names, params_mod, init_mod, options)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% NFKBSIMULATE initializes the model, makes requested modifications (to parameters and
% initial values), finds the basal model steady state, then simulates the response to the
% input stimulus.
%
% INPUTS:
% stimulus     (n x 2) cell matrix with stimulus name and value pairs (e.g. {'LPS',0.1})
% names        (n x 1) cell vector with output names to return
% params_mod   (n x 3) matrix that contains parameter index (rxn, subparm) and modified value
% init_mod     (n x 2) cell matrix with species name and initial value pairs (e.g. {'NFkB', 0.1})
% options      options structure - simulation-specific options (e.g. time step)
%
% OUTPUTS:
% t            time vector (0:v.DT:v.SIM_TIME)
% x            data matrix (1 column per entry in 'names')
% simdata      output data structure options, basal values; etc.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deleted


%global steadyyy_state i_ste

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization: set defaults, modify using inputs
if nargin <4
    init_mod = {};
    if nargin < 3
        params_mod = [];
        if nargin <2
            [~, tmp] = nfkbInitialize;
            names = tmp.NAMES;
            if nargin<1
                stimulus = {};
            end
        end
    end
end

% OPTIONS: (times and tolerances)
v.DT = 1;
v.SIM_TIME     = 18*60; %  length of stimulation phase (phase 2), in min

v.T_EQUILIBRATE   = -4000;
v.DEBUG = 0;
v.PULSE_TIME = inf;
% Defalut tolderances for ODE solver
v.RELTOL = 1e-6;
v.ABSTOL = 1e-8;



% Cycle input 'options' and modify
if nargin>4
    fields = fieldnames(options);
    for i = 1:length(fields)
        v.(fields{i}) = options.(fields{i});
    end
end
% If only length of simulation is specified, specify start time = 0
if length(v.SIM_TIME) == 1
    v.SIM_TIME = [0 v.SIM_TIME];
end
% Error checking: if an initial state isn't provided, start time MUST be 0
if (v.SIM_TIME(1) > 0) && ~isfield(v,'STEADY_STATE')
    error('Initial/steady state vals (e.g. options.STEADY_STATE) must be provided if not starting simulation from t=0')
end


ode_opt = odeset('RelTol',v.RELTOL,'AbsTol',v.ABSTOL);

% PARAMETERS AND INITIAL VALUES
if v.DEBUG;disp('INITIALIZATION: setting parameters and initial vals'); end
% [v.PARAMS, v.SPECIES] = nfkbInitialize();
v.PARAMS = options.v.PARAMS;
v.SPECIES = options.v.SPECIES;
if isfield(options.v,'init_val')
    v.init_val = options.v.init_val;
end


% Cycle input 'params_mod' and modify parameters
for i = 1:size(params_mod,1)
    v.PARAMS(params_mod(i,1),params_mod(i,2)) = params_mod(i,3);
    if v.DEBUG
        disp(['- param(', num2str(params_mod(i,1)),',',num2str(params_mod(i,2)),...
            ') set to ',num2str(params_mod(i,3)) ])
    end
end

% set the corelated parameter the same value
v.PARAMS(40,1) = v.PARAMS(36,1);
% v.PARAMS(61,1) = v.PARAMS(58,1);
% v.PARAMS(64,1) = v.PARAMS(58,1);
v.PARAMS(65,2) = v.PARAMS(52,2);

% Concentrations for fixed species (can cycle between forms, but aren't significantly synthesized or degraded)
% v.SPECIES.INIT  = zeros(length(v.SPECIES.NAMES),1); % all species start at zero by default
% v.SPECIES.INIT(strcmp('NFkB',v.SPECIES.NAMES))   = 0.08;
% v.SPECIES.INIT(strcmp('IKK_off',v.SPECIES.NAMES))   = 0.08;
% % TLR4-specific
% v.SPECIES.INIT(strcmp('MyD88_off',v.SPECIES.NAMES))   = 0.1;
% v.SPECIES.INIT(strcmp('TRIF_off',v.SPECIES.NAMES))   = 0.1;
% v.SPECIES.INIT(strcmp('TRAF6_off',v.SPECIES.NAMES))   = 0.1;
% % TNF-specific
% v.SPECIES.INIT(strcmp('TTR',v.SPECIES.NAMES))   = 0.035;
% v.SPECIES.INIT(strcmp('TAK1_off',v.SPECIES.NAMES))   = 0.1;

% Concentrations for fixed species (can cycle between forms, but aren't significantly synthesized or degraded)
v.SPECIES.INIT  = zeros(length(v.SPECIES.NAMES),1); % all species start at zero by default
if isfield(v,'init_val') && isfield(v.init_val,'NFkB')
    v.SPECIES.INIT(strcmp('NFkB',v.SPECIES.NAMES))   = v.init_val.NFkB;
else
    v.SPECIES.INIT(strcmp('NFkB',v.SPECIES.NAMES))   = 0.08;
end

if isfield(v,'init_val') && isfield(v.init_val,'IKK_off')
    v.SPECIES.INIT(strcmp('IKK_off',v.SPECIES.NAMES))   = v.init_val.IKK_off;
else
    v.SPECIES.INIT(strcmp('IKK_off',v.SPECIES.NAMES))   = 0.08;
end
% TLR4-specific
if isfield(v,'init_val') && isfield(v.init_val,'MyD88_off')
    v.SPECIES.INIT(strcmp('MyD88_off',v.SPECIES.NAMES))   = v.init_val.MyD88_off;
else
    v.SPECIES.INIT(strcmp('MyD88_off',v.SPECIES.NAMES))   = 0.1;
end

if isfield(v,'init_val') && isfield(v.init_val,'TRIF_off')
    v.SPECIES.INIT(strcmp('TRIF_off',v.SPECIES.NAMES))   = v.init_val.TRIF_off;
else
    v.SPECIES.INIT(strcmp('TRIF_off',v.SPECIES.NAMES))   = 0.1;
end

if isfield(v,'init_val') && isfield(v.init_val,'TRAF6_off')
    v.SPECIES.INIT(strcmp('TRAF6_off',v.SPECIES.NAMES))   = v.init_val.TRAF6_off;
else
    v.SPECIES.INIT(strcmp('TRAF6_off',v.SPECIES.NAMES))   = 0.1;
end
% TNF-specific

if isfield(v,'init_val') && isfield(v.init_val,'TTR')
    v.SPECIES.INIT(strcmp('TTR',v.SPECIES.NAMES))   = v.init_val.TTR;
else
    v.SPECIES.INIT(strcmp('TTR',v.SPECIES.NAMES))   = 0.035;
end

if isfield(v,'init_val') && isfield(v.init_val,'TAK1_off')
    v.SPECIES.INIT(strcmp('TAK1_off',v.SPECIES.NAMES))   = v.init_val.TAK1_off;
else
    v.SPECIES.INIT(strcmp('TAK1_off',v.SPECIES.NAMES))   = 0.1;
end

if isfield(v,'init_val') && isfield(v.init_val,'C1_off')
    v.SPECIES.INIT(strcmp('C1_off',v.SPECIES.NAMES))   = v.init_val.C1_off;
else
    v.SPECIES.INIT(strcmp('C1_off',v.SPECIES.NAMES))   = 0.1;
end

if isfield(v,'init_val') && isfield(v.init_val,'TRAF2_off')
    v.SPECIES.INIT(strcmp('TRAF2_off',v.SPECIES.NAMES))   = v.init_val.TRAF2_off;
else
    v.SPECIES.INIT(strcmp('TRAF2_off',v.SPECIES.NAMES))   = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added
%% v.PARAMS(strcmp('
% 
% if nargin>4 && isfield(options,'params_size') && isfield(options,'params_value')
%     for i = 1:size(options.params_size,1)
%         v.PARAMS(options.params_size(i,1),options.params_size(i,2)) = options.params_value(i);
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cycle input 'init_mod' and modify corresponding initial values for species
for i = 1:size(init_mod,1)
    idx = find(strcmp(init_mod{i,1},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
        v.SPECIES.INIT(idx)   = init_mod{i,2};
        if v.DEBUG
            disp(['- initial ',init_mod{i,1},' altered to ',num2str(init_mod{i,2})])
        end
    elseif v.DEBUG
        disp(['- species ''',init_mod{i,1},''' not found']);
    end
end


%% Simulations

% PHASE1: find basal steady state
max_attempts = 100;
max_pct = 0.0005;
min_small = 8e-7; % Ignore differences in smaller quantities than this (corresponds to ~1 molecule/cell)

if ~isfield(v,'STEADY_STATE')
    v.PHASE     = 1;
    [~, r1] = ode15s('nfkbOde_TNF_TRAF2', [v.T_EQUILIBRATE 0], v.SPECIES.INIT,ode_opt,v);
    starting_vals = r1(end,:);
    % Iterate until we reach steady state (convergence to < 0.01%)
    if v.DEBUG; disp('PHASE1: running model to basal steady state'); end
    converge = 0;
    iter = 0;
    %    nfkbOde([],[],[],v); % Create global delay variable
    while ~converge
        [~, r1] = ode15s('nfkbOde_TNF_TRAF2', [v.T_EQUILIBRATE 0], starting_vals,ode_opt,v);
        pct_diff = abs((starting_vals - r1(end,:))./r1(end,:));
        starting_vals = r1(end,:);
        pct_diff(starting_vals < min_small) = 0;
        if (max(pct_diff) < max_pct) || (iter>max_attempts)
            converge = 1;
            if (iter>max_attempts)
                idx = find(pct_diff==max(pct_diff),1,'first');
                disp(['Note: max steady state iterations reached (',num2str(max_attempts),'). Aborting with current state.'])
                disp(['[Non-convergent species: ', v.SPECIES.NAMES{idx},'(#',num2str(idx),')'...
                    '. Current state: ',num2str(starting_vals(idx)),' uM (',num2str(100*pct_diff(idx)),'% difference from prev.)]'])
            end
        end
        iter = iter+1;
    end
    v.STEADY_STATE = starting_vals;
    if v.DEBUG; disp(['- iterations required: ',num2str(iter),' (',num2str(-iter*v.T_EQUILIBRATE),' min)']); end
end
% steadyyy_state(i_ste,:)=v.STEADY_STATE;
% i_ste=i_ste+1;
% PHASE2: add stimulus, and simulate from basal levels.
if v.DEBUG; disp('PHASE2: simulating dynamic response to stimulus'); end
starting_vals = v.STEADY_STATE;

for i = 1:size(stimulus,1)
    idx = find(strcmp(stimulus{i,1},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
        starting_vals(idx)   = stimulus{i,2};
        if v.DEBUG
            disp(['- stimulus ''',stimulus{i,1},''' set to ',num2str(stimulus{i,2})])
        end
    elseif v.DEBUG
        disp(['- stimulus ',stimulus{i,1},' not found']);
    end
end

v.PHASE = 2;
[t_sim, x_sim] = ode15s('nfkbOde_TNF_TRAF2', v.SIM_TIME,starting_vals,ode_opt,v);

if max(abs(imag(x_sim))) > eps
    warn(['Some numerical solutions have imaginary components. (Max imag. value = ',num2str(max(abs(imag(x_sim)))),'i)'])
end
x_sim = real(x_sim);
x_end = x_sim(end,:);


%% Interpolate final outputs
x = [];
t = (v.SIM_TIME(1):v.DT:v.SIM_TIME(2))';

simdata = v;
simdata.STIMULUS = stimulus;

simdata.output = {};
for i = 1:length(names)
    idx = find(strcmp(names{i},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
        x = cat(2,x,interp1(t_sim,x_sim(:,idx),t));
        simdata.output = cat(1,simdata.output,names{i});
    elseif v.DEBUG
        disp(['Output species ',names{i}, ' not found - skipping'])
    end
end


if v.DEBUG; disp('- - - - - - - - - - - - - - - - - - - -'); end


