
    %% sampling
    % proj_num =18, TNF
    mu =        [0.676, 0.232,      0.0293,     0.245,      4.53e-6,    0.00423,     0.0455 ]; % mean value
    sigma =     [1.83,  2.47,       8.19,       2.06,       1.92,       1.91,        2.83]; %standard deviation of the ramdom effects
    corr_mat = [1,      -0.116      -0.234      0.0497      0.558       0.381       0.593
        -0.116  1           0.929       0.852       0.334       -0.134      -0.383
        -0.234  0.929       1           0.846       0.182       -0.154      -0.342
        0.0497  0.852       0.846       1           0.217       -0.0229     -0.155
        0.558   0.334       0.182       0.217       1           0.36        -0.0325
        0.381   -0.134      -0.154      -0.0229     0.36        1           0.216
        0.593   -0.383      -0.342      -0.155      -0.0325     0.216       1];
    para_min = [0.05188,0.01,       0.01,       0.01,       8.224e-7,   0.00116,    0.04];
    
    para_max = [5.188,  1,          1,          1,          8.224e-5,   0.116       0.3];
    
    sigma_mat = diag(sigma)*corr_mat*diag(sigma);
    
    mu_log = log((mu - para_min)./(para_max-mu));
    
    N_cells = Num_sample;
    
    R = mvnrnd(mu_log,sigma_mat,N_cells);
    
    para_val = (exp(R).*(ones(N_cells,1)*para_max) + ones(N_cells,1)*para_min)./(1+exp(R));