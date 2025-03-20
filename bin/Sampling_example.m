% firstDim = (1:5)'
% mu = repmat(firstDim,1,3)
% sigma = eye(3)
% rng('default') ; % For reproducibility
% R = mvnrnd(mu,sigma)

% mu = [2 3, 4];
% sigma = [1 1.5,0; 1.5 3,0;0,0,2];
rng('default')  % For reproducibility
% R = mvnrnd(mu,sigma,100);
% plot(R(:,1),R(:,2),'+')

% scatter3(R(:,1),R(:,2),R(:,3))

% TNF example: 

%correlation matrix:

%           p66,    p99,        p100,       p101,       p53,        p54
mu =        [0.676, 0.232,      0.0293,     0.245,      4.53e-6,    0.00423 ]; % mean value
sigma =     [1.83,  2.47,       8.19,       2.06,       1.92,       1.91]; %standard deviation of the ramdom effects
corr_mat = [1,      -0.116      -0.234      0.0497      0.381       0.558
            -0.116  1           0.929       0.852       -0.134      0.334
            -0.234  0.929       1           0.846       -0.154      0.182
            0.0497  0.852       0.846       1           -0.0229     0.217
            0.381   -0.134      -0.154      -0.0229     1           0.36
            0.558   0.334       0.182       0.217       0.36        1];
para_min = [0.05188,0.01,0.01,0.01,8.224e-7,0.00116];

para_max = [5.188,1,1,1,8.224e-5,0.116];

sigma_mat = diag(sigma)*corr_mat*diag(sigma);

mu_log = log((mu - para_min)./(para_max-mu));

N_cells = 500;

R = mvnrnd(mu_log,sigma_mat,N_cells);

para_val = (exp(R).*(ones(N_cells,1)*para_max) + ones(N_cells,1)*para_min)/(1+exp(R));