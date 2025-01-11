function [irs, ses, cis_dm, cis_boot, ses_bootstrap] = ir_estim(Y, p, horzs, varargin)

    % Wrapper function for (V)AR or LP estimation of impulse responses
    % Delta method and bootstrap confidence intervals
    
    % Inputs: see below
    
    % Outputs:
    % irs       1 x H       estimated impulse responses at select horizons
    % ses       1 x H       s.e. for impulse responses
    % cis_dm    2 x H       lower and upper limits of delta method confidence intervals
    % cis_boot  2 x H x 3   lower and upper limits of bootstrap confidence intervals (3rd index: type of interval, either Efron, Hall, or Hall percentile-t)
   
    
    %% Parse inputs
    
    ip = inputParser;
    
    % Required inputs
    addRequired(ip, 'Y', @isnumeric);
        % T x n     data vector
    addRequired(ip, 'p', @isnumeric);
        % 1 x 1     lag length (not counting any augmentation)
    addRequired(ip, 'horzs', @isnumeric);
        % 1 x H     impulse response horizons of interest
    
    % Optional inputs
    addParameter(ip, 'resp_ind', 1, @isnumeric);
        % Index of response variable of interest (default: first variable)
    addParameter(ip, 'innov_ind', 1, @isnumeric);
        % Index of innovation of interest (default: first innovation)
    addParameter(ip, 'estimator', 'lp', @ischar);
        % Estimator type, either 'var' or 'lp' (default: local projection)
    addParameter(ip, 'alpha', 0.05, @isnumeric);
        % Significance level (default: 0.05)
    addParameter(ip, 'bias_corr', true, @islogical);
        % Bias-correct VAR estimates? (default: yes)
    addParameter(ip, 'bias_corr_lp', true, @islogical)
        % Bias-correct LP estimates? (default: yes)
    addParameter(ip, 'se_homosk', false, @islogical);
        % Homoskedastic standard errors/bootstrap? (default: no)
    addParameter(ip, 'no_const', false, @islogical);
        % Omit intercept? (default: no)
    addParameter(ip, 'bootstrap', [], @(x) ischar(x) || isempty(x));
        % Bootstrap type, either 'var' or empty if delta method (default: delta method)
    addParameter(ip, 'boot_num', 1000, @isnumeric);
        % Bootstrap iterations (default: 1000)
    addParameter(ip, 'boot_workers', 0, @isnumeric);
        % Number of parallel workers used for bootstrapping (default: 0, meaning no parallel computation)    
    parse(ip, Y, p, horzs, varargin{:});
    
    
    %% Preliminaries
    
    nh ...
      = length(horzs); % Number of horizons
    
    cvs ...
      = repmat(norminv(1-ip.Results.alpha/2),...
              1, nh);      % Default critical values: normal

    cis_boot ...
      = nan(2,nh,3);       % Initializes NaN array
    
    
    %% Point estimates and var-cov
    
    if strcmp(ip.Results.estimator, 'var') % VAR
        
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
        
    elseif strcmp(ip.Results.estimator, 'lp') % LP
        
        irs = zeros(1,nh);
        
        ses = zeros(1,nh);
        
        betahat ...
            = cell(1,nh);
        
        res = cell(1,nh);
        
        X   = cell(1,nh);
        
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

        if ip.Results.bias_corr_lp  % Herbst & Johanssen (2024) bias correction           
            % Include contemporaneous/lagged controls and exclude the
            % contemporaneous innovation
            w   = [Y(:, 1:ip.Results.innov_ind-1), lagmatrix(Y, 1:p)];
            w   = w(p+1:end, :);            
            irs = LP_CorrectBias(irs, w);
        end
        
    end
    
    % If only point estimates and standard errors are requested, stop
    if nargout<=2
        return;
    end
    
    
    %% Delta method confidence intervals

    cis_dm = irs + [-1; 1]*(cvs.*ses);

    %% Bootstrap confidence intervals
    
    if ~isempty(ip.Results.bootstrap)
        
        estims_boot = zeros(ip.Results.boot_num,nh);
        ses_boot = estims_boot;
        
        if strcmp(ip.Results.bootstrap, 'var') % Recursive VAR bootstrap (only option)
        
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
            Y_boot = var_boot(Ahat_var, res_var, Y, p, ip.Results.se_homosk, ip.Results.no_const);
            
            % Estimate on bootstrap sample
            [estims_boot(b,:),ses_boot(b,:)] = ir_estim(Y_boot, p, horzs, varargin{:});

        end
        
        % Compute bootstrap confidence intervals
        cis_boot = boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, ip.Results.alpha);
        ses_bootstrap = std( estims_boot);
        
        end

    end