function [irs, irs_varcov, Ahat_estim, Sigmahat_estim, res_estim] ...
            = var_ir_estim(Y, innov_ind, p, horzs, bias_corr, homosk, no_const)
    
    % SVAR(p) least-squares estimates and delta method s.e.
    % Cholesky identification
    
    % Inputs:
    % Y           T x n   data vector
    % innov_ind   1 x 1   index of Cholesky-orthogonalized innovation of interest
    % p           1 x 1   lag length used for impulse response computations
    % horzs       H x 1   horizons of interest
    % bias_corr   bool    true: apply analytical bias correction (Pope, 1990)
    % homosk      bool    true: homoskedastic s.e., false: EHW s.e.
    % no_const    bool    true: omit intercept
    % T_eff       scalar  Specify effective sample size. (=[] for default).

    % Outputs:
    % irs            n x H          estimated impulse responses Theta_h*nu at select horizons
    % irs_varcov     n x n x H      var-cov matrix of Theta_h*nu at select horizons
    % Ahat_estim     n x np         VAR coefficient estimates [A_1,...,A_p] (possibly bias-corrected, possibly including intercept as last column)
    % Sigmahat_estim n x n          VAR residual variance-covariance matrix
    % res_estim      (T-p) x n      estimation residuals
    
    [T,n] = size(Y);
    
    % One-step forecasting regression of Y_{t+1} on (Y_t, ..., Y_{t-p+1})
    [Ahat_estim, Ahat_estim_varcov, Sigmahat_estim, res_estim] ...
        = var_estim(Y, p, homosk, no_const);
    
    % Bias correction, if desired
    if bias_corr
        Ahat_estim(:,1:end-1+no_const) = var_biascorr(Ahat_estim(:,1:end-1+no_const), Sigmahat_estim, T);
    end
    
    % Only use first p VAR coefficient matrices to compute impulse responses
    Ahat = Ahat_estim(:,1:n*p);
    Ahat_varcov = Ahat_estim_varcov(1:n^2*p,1:n^2*p);

    % Cholesky shock vector (implemented via numerically equivalent regression)
    [coef,coef_varcov] = linreg(res_estim,res_estim(:,1:innov_ind),homosk,true);
    nu = coef(:,innov_ind);
    nu_varcov = coef_varcov(end-n+1:end,end-n+1:end);

    % Degrees-of-freedom correction to account for additional controls in original VAR regression
    nu_varcov = nu_varcov*(size(res_estim,1)-innov_ind)/(size(res_estim,1)-size(Ahat_estim,2)-innov_ind);
    
    if nargout==1
        irs = var_ir(Ahat,nu,horzs); % Compute impulse responses
    else
        [irs, jacob_a, jacob_nu] = var_ir(Ahat,nu,horzs); % Compute impulse responses and Jacobians
        nh = length(horzs);
        irs_varcov = zeros(n,n,nh);
        for h=1:nh
            irs_varcov(:,:,h) = jacob_a(:,:,h)*Ahat_varcov*jacob_a(:,:,h)' ...
                + jacob_nu(:,:,h)*nu_varcov*jacob_nu(:,:,h)';
        end
    end

end