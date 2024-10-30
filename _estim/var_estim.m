function [betahat, betahat_varcov, Sigmahat, res] ...
            = var_estim(Y, num_lags, se_setting, no_const)

    % Reduced-form Vector Autoregression
    
    % Inputs:
    % Y           T x n      data matrix
    % num_lags    1 x 1      number of lags in VAR
    % se_setting             EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                        OR function handle: function that returns HAC/HAR sandwich matrix
    % no_const    bool       true: omit intercept
    
    % Outputs:
    % betahat           n x (n*num_lags+n+~no_const)      full vector of estimated regression coefficients
    % betahat_varcov    (n x (n*num_lags+n+~no_const))x   var-cov of vec(betahat)
    %                   (n x (n*num_lags+n+~no_const))
    % Sigmahat          n x n                             residual variance-covariance matrix
    % res               (T-p) x m                         residuals
    

    % Covariate matrix
    Y_lag = lagmatrix(Y,0:num_lags-1);
    X = Y_lag(num_lags:end-1,:);
    
    % VAR
    [betahat, betahat_varcov, res, ~] = linreg(Y(num_lags+1:end,:),X,se_setting,no_const);
    Sigmahat = (res'*res)/(size(res,1)-size(betahat,2));

end