%% ENCOMPASSING DFM

DF_model.reorder    = [1:76, 87:94, 77:86, 95:171, 181:195, 172:180, 196:207]; % index to reorder data to match variable list in Stock-Watson (2016)

DF_model.levels     = 0; % =1: variables in levels, 0=: first differences
DF_model.n_fac      = 6; % number of factors
DF_model.coint_rank = 2; % cointegration rank in factor process (for levels specification); []: estimate using Johansen test
DF_model.n_lags_fac = 2; % lag order of factors
DF_model.n_lags_uar = 2; % lag order of measurement error

DF_model.censor_arch_uar = 0.7;  % cutoff for right-censoring ARCH parameter when shock_type='arch'

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

settings.simul.n_mc    = 1000; % number of Monte Carlo reps
settings.simul.seed    = (1:settings.simul.n_mc)*10 + randi([0,9],1,settings.simul.n_mc); % random seed for each Monte Carlo

% time periods for each simulation
switch sample_length
    case 'short'
        settings.simul.T = 100; 
    case 'medium'
        settings.simul.T = 240; 
    case 'long'
        settings.simul.T = 720; 
end
settings.simul.T_burn = 100; % burn-in

% degree of misspecification

settings.misspec.indic        = 1; % indicator for computing the degree of mis-specification
settings.misspec.lags         = [2 4 8]; % lags for mis-specification analysis
settings.misspec.n_lags       = length(settings.misspec.lags);
settings.misspec.VAR_infinity = 100; % VAR lag length
settings.misspec.VMA_hor      = 200; % maximal VMA horizon
settings.misspec.zeta         = 0.5; % local-to-VAR coefficient

%% ESTIMATION SETTINGS

% current list, for recursive:
% VAR AIC, VAR 4, VAR 8, VAR AIC (homoscedastic), LP AIC, LP 4,
% LP AIC (homoscedastic), BVAR, default SLP, under-smoothed SLP

% current list, for observed shock:
% VAR AIC, VAR 4, VAR 8, VAR AIC (homoscedastic), LP AIC, LP 4, 
% LP AIC (homoscedastic), BVAR, default SLP, under-smoothed SLP
% VAR AIC small, LP AIC small

% -------------------------------------------------------------------------
% Smoothed LP and BVAR settings
% -------------------------------------------------------------------------

% Smoothed LP defaults

opts_slp.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, ...
                          2:1:19, 20:20:100, 200:200:2000];  % CV grid, scaled by T
opts_slp.CV_folds      = 5;                                  % # CV folds
opts_slp.irfLimitOrder = 2;                                  % shrink towards polynomial of that order
opts_slp.undersmooth   = false;                              % multiply optimal lambda by 0.1? 

% (Under-)smoothed LP defaults

opts_slp_undersmooth             = opts_slp;
opts_slp_undersmooth.undersmooth = true;

% BVAR defaults

if mode_type == 5
    opts_bvar.RW = true;    % random walk prior
else
    opts_bvar.RW = false;    % WN prior
end
opts_bvar.ndraw = 500;  % Posterior draws

%----------------------------------------------------------------
% Add to methods
%----------------------------------------------------------------

% VAR AIC, VAR 4, VAR 8 (EHW, block bootstrap)
[settings.est.methods{1:3}] = deal({'resp_ind' , [], ...
    'innov_ind', [], ...
    'bootstrap', 'var',...
    'estimator', 'var', ...
    'se_homosk', false,...
    'boot_blocklength', ceil(5.03*settings.simul.T^(1/4)),...
    'shrinkage', false});

% VAR AIC (homoscedastic, recursive bootstrap)
[settings.est.methods{4}] = deal({'resp_ind' , [], ...
    'innov_ind', [], ...
    'bootstrap', 'var',...
    'estimator', 'var', ...
    'se_homosk', true,...
    'boot_blocklength', 1,...
    'shrinkage', false});

% LP AIC, LP 4 (EHW, block bootstrap)
[settings.est.methods{5:6}] = deal({'resp_ind'  , [], ...
    'innov_ind' , [], ...
    'bootstrap', 'var',...
    'estimator' , 'lp' , ...
    'se_homosk', false,...
    'boot_blocklength', ceil(5.03*settings.simul.T^(1/4)),...
    'shrinkage' , false});  % LP

% LP AIC (homoscedastic, recursive bootstrap)
[settings.est.methods{7}] = deal({'resp_ind'  , [], ...
    'innov_ind' , [], ...
    'bootstrap', 'var',...
    'estimator' , 'lp' , ...
    'se_homosk', true,...
    'boot_blocklength', 1,...
    'shrinkage' , false});  % LP


% BVAR
settings.est.methods{8} = {'resp_ind',  [], ...
    'innov_ind', [],...
    'estimator', 'var',...
    'bootstrap', [],...
    'shrinkage', true,...
    'opts_bvar', opts_bvar};

% SLP, default
settings.est.methods{9} = {'resp_ind',  [], ...
    'innov_ind', [],...
    'estimator', 'lp', ...
    'bootstrap', [],...
    'shrinkage', true, ...
    'opts_slp',  opts_slp};

% SLP, under-smoothed
settings.est.methods{10} = {'resp_ind',  [], ...
    'innov_ind', [],...
    'estimator', 'lp', ...
    'bootstrap', [],...
    'shrinkage', true, ...
    'opts_slp',  opts_slp_undersmooth};


% Set lags
settings.est.n_lags_fix          = [NaN; 4; 8; NaN; % VAR 
                                    NaN; 4; NaN;    % LP
                                    4; 4; 4];       % Shrinkiage
settings.est.system_type         =  repmat({'big'}, 10, 1);


if strcmp(estimand_type, 'obsshock')
    settings.est.methods{11} = settings.est.methods{1};  % Small VAR AIC
    settings.est.methods{12} = settings.est.methods{5};  % Small LP AIC
    settings.est.n_lags_fix  = [settings.est.n_lags_fix; NaN; NaN];
    settings.est.system_type = [settings.est.system_type; repmat({'small'}, 2, 1)];
end

settings.est.est_n_lag  = isnan(settings.est.n_lags_fix);  % indicator if lags are estimated
settings.est.n_lags_max = 10;
settings.est.n_methods  = length(settings.est.methods); % number of estimation methods

%----------------------------------------------------------------
% Shared settings
%----------------------------------------------------------------

settings.est.alpha            = 0.1; % significance level. 1-alpha credible interval for BVAR.
settings.est.no_const         = false; % true: omit intercept
settings.est.boot_num         = 500;   % number of bootstrap samples
settings.est.bias_corr_var    = true;  % Pope (1990) VAR bias correction
settings.est.bias_corr_lp     = true;  % HJ (2024) LP bias correction

settings.est.methods_shared = {
                               'alpha'        , settings.est.alpha,...
                               'no_const'     , settings.est.no_const,...
                               'boot_num'     , settings.est.boot_num,...
                               'bias_corr_var', settings.est.bias_corr_var,...
                               'bias_corr_lp' , settings.est.bias_corr_lp};


%% PARALLELIZATION

settings.sim.num_workers = 8;