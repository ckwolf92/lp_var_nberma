function [irs, irs_varcov, betahat, betahat_varcov, res, X] ...
            = lp_ir_estim(Y, num_lags, horz, resp_ind, innov_ind, se_setting, no_const)

    % Local projection with contemporaneous lagged controls
    
    % Inputs:
    % Y         T x n       data matrix
    % num_lags  1 x 1       number of lags of Y to control for (in addition to the contemporaneous regressors)
    % horz      1 x 1       horizon of interest
    % resp_ind  1 x 1       index of response variables Y_t of interest
    % innov_ind 1 x 1       index of shock variable of interest
    % se_setting            EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                       OR function handle: function that returns HAC/HAR sandwich matrix
    % no_const   bool       true: omit intercept
    
    % Outputs:
    % irs               1 x 1                                    estimated impulse responses at selected horizon
    % irs_varcov        1 x 1                                    var-cov of irs
    % betahat           (n*num_lags+innov_ind+~no_const)         full vector of estimated regression coefficients
    % betahat_varcov    ((n*num_lags+innov_ind+~no_const))x      var-cov of vec(betahat)
    %                   ((n*num_lags+innov_ind+~no_const))
    % res               (T-p) x 1                                residuals
    % X                 (T-p) x (n*num_lags+innov_ind+~no_const) covariate data matrix (expanded if intercept included)
    
    % System sizes
    n = size(Y,2);

    % Covariate matrix
    Y_lag = lagmatrix(Y,0:num_lags);
    X = [Y_lag(num_lags+1:end-horz,1:innov_ind), Y_lag(num_lags+1:end-horz,n+1:end)];
    
    % Run local projection
    [betahat, betahat_varcov, res, X] = linreg(Y(num_lags+horz+1:end,resp_ind),X,se_setting,no_const);
    
    irs = betahat(:,innov_ind);
    irs_varcov = betahat_varcov(innov_ind,innov_ind);

end