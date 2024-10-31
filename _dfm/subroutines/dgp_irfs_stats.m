function model = dgp_irfs_stats(model, settings, estimand_type)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

var_select = settings.specifications.var_select;
[n_spec,n_var] = size(var_select);

model.LRV_Cov_tr_ratio = nan(n_spec,1);
model.dLRV_dCov_tr_ratio = nan(n_spec,1);
model.VAR_largest_root = nan(n_spec,1);
model.VAR_quant_root = nan(n_spec,1);
model.frac_coef_for_large_lags = nan(n_spec,1);
model.R0_sq = nan(n_spec,1);

%----------------------------------------------------------------
% IRFs in DFM
%----------------------------------------------------------------

model.irf = compute_irfs(model.ABCD,settings.est.shock_weight,settings.est.IRF_hor);

if any(strcmp(estimand_type, {'obsshock'})) % observed shock
    model.normalized_irf = compute_normalized_irfs(model,settings);
    model.target_irf = model.normalized_irf(settings.est.IRF_select, :);
else % recursive
    model.VAR_irf = nan(length(settings.est.IRF_select),n_spec); % To be computed below
    model.target_irf = nan(length(settings.est.IRF_select),n_spec);
end

%----------------------------------------------------------------
% Statistics for Individual DGPs
%----------------------------------------------------------------

for i_spec = 1:n_spec

    ABCD_obs = model.ABCD;
    ABCD_obs.C = ABCD_obs.C(var_select(i_spec,:),:);
    ABCD_obs.D = ABCD_obs.D(var_select(i_spec,:),:);

    % compute reduced-form VAR(infinity) representation

    [ABCD_small, shock_select] = ABCD_reduce(ABCD_obs);

    red_form = reduced_form_VAR(ABCD_small,settings.est.VAR_infinity_truncate);

    % Cholesky impulse responses for recursive identification

    if strcmp(estimand_type, 'recursive')
        G = chol(red_form.innov_var, 'lower');
        chol_irfs = compute_irfs(red_form.innov_ABCD,G(:,settings.est.recursive_shock_pos),settings.est.IRF_hor);
        model.VAR_irf(:,i_spec) = chol_irfs(settings.est.IRF_select,settings.est.IRF_response_var_pos);
        model.target_irf(:,i_spec) = model.VAR_irf(:,i_spec)/chol_irfs(1,settings.est.est_normalize_var_pos);
    end

    % roots of reduced-form VAR

    nlag_comp = settings.est.VAR_infinity_truncate_comp; % # lags to include in companion matrix
    comp_form = [cell2mat(red_form.coef(1:nlag_comp)); eye(n_var*(nlag_comp-1)) zeros(n_var*(nlag_comp-1), n_var)]; % Companion matrix
    roots = abs(eig(comp_form));
    model.VAR_largest_root(i_spec) = max(roots); % Max
    model.VAR_quant_root(i_spec) = quantile(roots,settings.est.VAR_root_quant); % Percentile

    % tr(LRV)/tr(Var) ratio

    if model.VAR_largest_root(i_spec)<0.999
        [cov, lrv] = cov_lrv(ABCD_small);
        model.LRV_Cov_tr_ratio(i_spec) = trace(lrv)/trace(cov);
    else
        dABCD = ABCD_diff(ABCD_small); % ABCD representation of first differences
        [dcov, dlrv] = cov_lrv(dABCD);
        model.dLRV_dCov_tr_ratio(i_spec) = trace(dlrv)/trace(dcov);
    end

    % fraction of VAR coefficients at long lags

    norms = cellfun(@(x) norm(x,'fro'), red_form.coef); % Frobenius norm of each VAR coefficient matrix at lags 1,2,...
    model.frac_coef_for_large_lags(i_spec) = 1-sum(norms(1:settings.est.n_lag_large_ref))/sum(norms);

    % degree of invertibility
    
    if ~strcmp(estimand_type, 'recursive')
        model.R0_sq(i_spec) = degree_invertibility(ABCD_small.D, red_form.innov_var, settings.est.shock_weight(shock_select));
    else
        model.R0_sq(i_spec) = 1;
    end

end

end