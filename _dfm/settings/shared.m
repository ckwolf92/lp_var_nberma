%% ENCOMPASSING DFM

DF_model.reorder    = [1:76, 87:94, 77:86, 95:171, 181:195, 172:180, 196:207]; % index to reorder data to match variable list in Stock-Watson (2016)

DF_model.levels     = 0; % =1: variables in levels, 0=: first differences
DF_model.n_fac      = 6; % number of factors
DF_model.coint_rank = 2; % cointegration rank in factor process (for levels specification); []: estimate using Johansen test
DF_model.n_lags_fac = 2; % lag order of factors
DF_model.n_lags_uar = 2; % lag order of measurement error

%% PREPARATIONS FOR STRUCTURAL ESTIMANDS

% selection of DGPs from encompassing model

settings.specifications.random_select         = 1; % randomly select variables from DFM list?
settings.specifications.random_n_spec         = 100; % number of random specifications
settings.specifications.random_n_var          = 5; % number of variables in each random specification
settings.specifications.random_category_range = [1 20; 21 31; 32 76; 77 86; 87 94; 95 131; 132 141;...
                                                 142 159; 160 171; 172 180; 181 207]; % ranges for Stock-Watson variable categories (see their Table 1)
settings.specifications.random_category_setup = {[1,2,3], 6}; % at least draw one from certain categories
settings.specifications.random_from_key_series = 0; % randomly select from some key series in DFM list?
settings.specifications.random_key_series = [1,2,6,12,56,95,97,121,132,142,147,151,172,181,193,196,202]; % list of key series

settings.specifications.manual_var_select     = [1 142; 1 97]; % if manual selection, then these sets of variables will be selected

% preliminary structural shock settings for observed shock

settings.est.estimate_shock_weight    = 1; % do we estimate the loading of the structural shock on the reduced-form shocks?
settings.est.manual_shock_pos         = 1; % if loading is not estimated, we just mechanically choose the xth shock in the state equation

% IRFs of interest

settings.est.IRF_hor    = 21; % maximal horizon (include contemporary)
settings.est.IRF_select = 1:settings.est.IRF_hor; % which IRFs to study
settings.est.n_IRF      = size(settings.est.IRF_select,2);

% compute roots and VAR(p) fit in population using truncated infinite-order VAR 

settings.est.VAR_infinity_truncate = 1000; % truncation horizon for VAR(infinity)
settings.est.VAR_infinity_truncate_comp = 200; % use only this many lags to construct companion matrix of VAR(infinity)
settings.est.VAR_root_quant = 0.75; % store this quantile of the roots of the companion matrix (in addition to largest)

% reference lags for long-lag statistics

settings.est.n_lag_large_ref = 4;

% number of Monte Carlo draws

settings.simul.n_mc    = 500; % number of Monte Carlo reps
settings.simul.seed    = (1:settings.simul.n_mc)*10 + randi([0,9],1,settings.simul.n_mc); % random seed for each Monte Carlo

% sample settings

settings.simul.T      = 240; % time periods for each simulation
settings.simul.T_burn = 100; % burn-in

%% ESTIMATION SETTINGS

% estimation methods

settings.est.methods{1} = {'estimator', 'var',...
            'bias_corr', true};

settings.est.methods{2} = {'estimator', 'lp',...
            'bias_corr', true};

settings.est.no_const  = false; % true: omit intercept
settings.est.se_homosk = true; % true: homoskedastic ses
settings.est.alpha     = 0.1; % significance level
settings.est.boot_num  = 1e3;  % number of bootstrap samples
settings.est.bootstrap = 'var';  % VAR bootstrap

settings.est.methods_shared = {'resp_ind',  [], ...
                               'innov_ind', [], ...
                               'alpha',     settings.est.alpha,...
                               'no_const',  settings.est.no_const,...
                               'se_homosk', settings.est.se_homosk,...
                               'boot_num',  settings.est.boot_num};

% lag length selection

settings.est.est_n_lag  = isnan(lag_type);
settings.est.n_lags_max = 10;
settings.est.n_lags_fix = lag_type;

%% PARALLELIZATION

settings.sim.num_workers = 8;