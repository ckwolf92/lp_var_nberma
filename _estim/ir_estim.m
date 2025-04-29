function [irs, ses, cis, cis_boot, ses_bootstrap] = ir_estim(Y, p, horzs, varargin)

% Wrapper function for (Bayesian) VAR, (smoothed) LP estimation of impulse responses

% Inputs: see the below "PARSE INPUTS" code section

% Outputs:
% irs            1 x H       estimated impulse responses at select horizons
% ses            1 x H       s.e. for impulse responses. Delta method for
%                            VAR, LP, and smoothed LP. Standard deviation of 
%                            posterior distribution for BVAR.
% cis            2 x H       lower and upper limits of delta method confidence
%                            intervals (for VAR, LP, and smoothed LP). 
%                            Equal-tailed credible interval for BVAR.
% cis_boot       2 x H x 3   lower and upper limits of bootstrap confidence 
%                            intervals (3rd index: type of interval, either 
%                            Efron, Hall, or Hall percentile-t)
% ses_bootstrap  1 x H       Standard deviation of bootstrap distribution.

%% SHRINKAGE ESTIMATOR DEFAULT SETTINGS 
% BVAR and smoothed local projection

% Smoothed local projection 
opts_slp.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, ...
                          2:1:19, 20:20:100, 200:200:2000];  % CV grid, scaled by T
opts_slp.CV_folds      = 5;      % # CV folds
opts_slp.irfLimitOrder = 2;      % shrink towards polynomial of that order
opts_slp.undersmooth   = false;  % multiply CV lambda by .1? (default: no)

% Bayesian VAR 
opts_bvar.RW    = true;  % Random walk prior? (default: yes)
opts_bvar.ndraw = 500;   % Number of posterior draws

%% PARSE INPUTS

ip = inputParser;

% -------------------------------------------------------------------------
% Required inputs
% -------------------------------------------------------------------------
addRequired(ip, 'Y', @isnumeric);
% T x n     data vector
addRequired(ip, 'p', @isnumeric);
% 1 x 1     lag length (not counting any augmentation)
addRequired(ip, 'horzs', @isnumeric);
% 1 x H     impulse response horizons of interest

% -------------------------------------------------------------------------
% Optional inputs: General settings
% -------------------------------------------------------------------------
addParameter(ip, 'resp_ind', 1, @isnumeric);
% Index of response variable of interest (default: first variable)
addParameter(ip, 'innov_ind', 1, @isnumeric);
% Index of innovation of interest (default: first innovation)
addParameter(ip, 'estimator', 'lp', @ischar);
% Estimator type, either 'var' or 'lp' (default: local projection)
addParameter(ip, 'shrinkage', false, @islogical) 
% Shrink using smoothed LP or Bayesian VAR? (default: no)
addParameter(ip, 'se_homosk', false, @islogical);
% Homoskedastic standard errors/bootstrap? (default: no)
addParameter(ip, 'no_const', false, @islogical);
% Omit intercept? (default: false)
addParameter(ip, 'alpha', 0.05, @isnumeric);
% Significance level (default: 0.05)
addParameter(ip, 'bias_corr', true, @islogical);
% Bias-correct VAR or LP estimates? (default: yes)

% -------------------------------------------------------------------------
% Optional inputs: Shrinkage and bias correction settings
% -------------------------------------------------------------------------
addParameter(ip, 'opts_slp', opts_slp)
% Options for smoothed LP (default: see above)
addParameter(ip, 'opts_bvar', opts_bvar)
% Options for BVAR (default: see above)

% -------------------------------------------------------------------------
% Optional inputs: Bootstrap settings
% -------------------------------------------------------------------------
addParameter(ip, 'bootstrap', [], @(x) ischar(x) || isempty(x));
% 'var' for VAR bootstrap or empty if delta method (default: delta method)
addParameter(ip, 'boot_num', 1000, @isnumeric);
% Bootstrap iterations (default: 1000)
addParameter(ip, 'boot_blocklength', 1)
% >=1 for block bootstrap and 'wild' for wild (default 1)
addParameter(ip, 'boot_workers', 0, @isnumeric);
% Number of parallel workers used for bootstrapping (default: 0, meaning no parallel computation)
parse(ip, Y, p, horzs, varargin{:});


%% PRELIMINARIES

nh       = length(horzs); % Number of horizons
cvs      = repmat(norminv(1-ip.Results.alpha/2), 1, nh);  % Default critical values: normal
cis_boot = nan(2,nh,3);       % Initializes NaN array

%% POINT ESTIMATES AND VAR-COV

% -------------------------------------------------------------------------
% VAR, no shrinkage
% -------------------------------------------------------------------------
if strcmp(ip.Results.estimator, 'var') & ~ip.Results.shrinkage 

    % VAR impulse responses
    [irs_all, irs_all_varcov] = var_ir_estim(Y, ...
        ip.Results.innov_ind, ...
        p,...
        horzs, ...
        ip.Results.bias_corr,...
        ip.Results.se_homosk,...
        ip.Results.no_const);

    % Impulse responses of interest and s.e.
    irs = irs_all(ip.Results.resp_ind,:);
    ses = sqrt(reshape(irs_all_varcov(ip.Results.resp_ind,ip.Results.resp_ind,:),1,[]));

% -------------------------------------------------------------------------
% LP, no shrinkage
% -------------------------------------------------------------------------
elseif strcmp(ip.Results.estimator, 'lp') & ~ip.Results.shrinkage 

    irs     = zeros(1,nh);
    ses     = zeros(1,nh);
    betahat = cell(1,nh);
    res     = cell(1,nh);
    X       = cell(1,nh);

    for h=1:nh % For each horizon...

        the_horz = horzs(h); % Horizon

        % LP regression
        [the_irs_all, the_irs_all_varcov, betahat{h},...
            ~, res{h}, X{h}] = lp_ir_estim(Y, ...
            p,...
            the_horz,...
            ip.Results.resp_ind,...
            ip.Results.innov_ind,...
            ip.Results.se_homosk,...
            ip.Results.no_const);

        % Impulse responses of interest and s.e.
        irs(h) = the_irs_all;
        ses(h) = sqrt(the_irs_all_varcov);

    end

    if ip.Results.bias_corr  % Herbst & Johannsen (2024) bias correction
        % Include contemporaneous/lagged controls and exclude the
        % contemporaneous innovation
        w   = [Y(:, 1:ip.Results.innov_ind-1), lagmatrix(Y, 1:p)];
        w   = w(p+1:end, :);
        irs = lp_biascorr(irs, w);
    end

% -------------------------------------------------------------------------
% VAR, shrinkage (Bayesian VAR)
% -------------------------------------------------------------------------    
elseif strcmp(ip.Results.estimator, 'var') & ip.Results.shrinkage  
    
    n_Y = size(Y,2);

    if ip.Results.opts_bvar.RW  % Random walk prior
        r = bvarGLP(Y, p, 'MNpsi', 0, 'Fcast', 0);
    else  % White noise prior
        r = bvarGLP(Y, p, 'MNpsi', 0, 'Fcast', 0, ...
            'sur', 0, 'noc', 0, 'posi', 1:n_Y);
    end

    if ip.Results.opts_bvar.ndraw == 0 % only use posterior means

        beta_draws  = r.postmax.betahat;
        sigma_draws = r.postmax.sigmahat;

    else % draw from posterior

        beta_draws  = nan(1+n_Y*p, n_Y, ip.Results.opts_bvar.ndraw);
        sigma_draws = nan(n_Y, n_Y, ip.Results.opts_bvar.ndraw);

        for j=1:ip.Results.opts_bvar.ndraw
            [beta_draws(:,:,j),sigma_draws(:,:,j)] = ...
                post_draw(r.postmax.betahat, r.postmax.Sinv, ...
                r.postmax.cholZZinv,r.postmax.T);
        end        
    end

    % Compute IRFs
    irs_draws = nan(n_Y, length(horzs), ip.Results.opts_bvar.ndraw);

    for j=1:ip.Results.opts_bvar.ndraw
        G                = chol(sigma_draws(:,:,j), 'lower');
        ShockVector      = G(:,ip.Results.innov_ind);
        irs_draws(:,:,j) = var_ir(squeeze(beta_draws(2:end,:,j))', ShockVector, horzs);
        irs_draws(:,:,j) = irs_draws(:,:,j)./irs_draws(ip.Results.innov_ind,1,j);
    end

    % Organize output
    irs_all = mean(irs_draws, 3); % posterior mean of IRF draws (or just single IRF if ndraw=0)
    cis_all = quantile(irs_draws, [ip.Results.alpha/2 1-ip.Results.alpha/2], 3);  % Equal-tailed
    irs     = irs_all(ip.Results.resp_ind, :);
    cis     = squeeze(cis_all(ip.Results.resp_ind, :, :))';
    ses     = std(irs_draws(ip.Results.resp_ind,:,:), [], 3);

% -------------------------------------------------------------------------
% LP, shrinkage
% -------------------------------------------------------------------------
elseif strcmp(ip.Results.estimator, 'lp') & ip.Results.shrinkage 

    % Covariate matrix
    Y_lag = lagmatrix(Y, 0:p);
    Y_lag = Y_lag(p+1:end, :);    
    y     = Y_lag(:, ip.Results.resp_ind);
    x     = Y_lag(:, ip.Results.innov_ind);
    w     = Y_lag(:, [1:ip.Results.innov_ind-1, size(Y,2)+1:end]);  % Contemporaneous and lagged controls    

    % Setup
    lambdaRange = ip.Results.opts_slp.lambdaRange * size(Y,1); % scale the grid of lambda (penalty strength) by # time periods
    lambdaRange = [1e-4, lambdaRange, 1e10];                   % allow regular OLS or completely smoothed
    r           = opts_slp.irfLimitOrder+1;
    H_min       = 0;
    H_max       = max(ip.Results.horzs);
    
    % Leave-one-out CV 
    rss_cv             = locproj_cv(y, x, w, H_min, H_max, r, ...
        lambdaRange, ip.Results.opts_slp.CV_folds);
    [~,lambda_opt_loc] = min(rss_cv);
    lambda_opt         = lambdaRange(lambda_opt_loc); % optimally tuned lambda

    if ip.Results.opts_slp.undersmooth
        lambda_opt = lambda_opt*.1;
    end

    % Estimate
    [irs,~,~,ses] = locproj(y, x, w, H_min, H_max, r, lambda_opt);
    irs           = irs(:)';
    ses           = ses(:)';

end

% If only point estimates and standard errors are requested, stop
if nargout<=2
    return;
end

%% DELTA METHOD CONFIDENCE INTERVALS
% For (smoothed) LP and VAR

if ~(strcmp(ip.Results.estimator, 'var') & ip.Results.shrinkage)  
    cis = irs + [-1; 1]*(cvs.*ses);
end

%% BOOTSTRAP CONFIDENCE INTERVALS

if ~isempty(ip.Results.bootstrap)

    estims_boot = zeros(ip.Results.boot_num,nh);
    ses_boot    = estims_boot;

    if strcmp(ip.Results.bootstrap, 'var') % Recursive VAR bootstrap. Only option.

        % VAR coefficient estimates that define bootstrap DGP
        [irs_var, ~, Ahat_var, ~, res_var] = var_ir_estim(Y, ...
            ip.Results.innov_ind, ...
            p,...
            horzs, ...
            ip.Results.bias_corr,...
            ip.Results.se_homosk,...
            ip.Results.no_const);

        pseudo_truth = irs_var(ip.Results.resp_ind,:); % Pseudo-true impulse responses in bootstrap DGP

        parfor(b=1:ip.Results.boot_num, ip.Results.boot_workers)

            % Generate bootstrap sample based on (possibly lag-augmented) VAR estimates
            Y_boot = var_boot(Ahat_var, res_var, Y, p, ...
                ip.Results.boot_blocklength, ip.Results.no_const);

            % Estimate on bootstrap sample
            [estims_boot(b,:),ses_boot(b,:)] = ir_estim(Y_boot, p, horzs, varargin{:});

        end    
    end

    % Compute bootstrap confidence intervals
    cis_boot      = boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, ip.Results.alpha);
    ses_bootstrap = std(estims_boot);

end