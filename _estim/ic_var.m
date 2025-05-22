function p_select = ic_var(y,p_max,method)

    % VAR lag length selection
    
    % Inputs:
    % y         T x n_v       data vector
    % p_max     1 x 1         maximal lag length to consider
    % method    1 x 1         lag length selection method (1 = AIC, 2 = BIC)
    
    % Output:
    % p_select  1 x 1         selected lag length

    
    n_v = size(y,2);
    T   = size(y,1);
    bic = zeros(1,p_max);
    aic = zeros(1,p_max);
    
    % go through from one lag to max number of lags
    for p = 1:p_max
        [~, ~, Sigma] = var_estim(y, p, true, false);
        bic(1,p) = log(det(Sigma)) + (n_v^2 * p + n_v) * log(T - p_max) / (T - p_max);
        aic(1,p) = log(det(Sigma)) + (n_v^2 * p + n_v) * 2 / (T - p_max);
    end
    
    if method == 1
        [~,p_select] = min(aic);
    elseif method == 2
        [~,p_select] = min(bic);
    end

end