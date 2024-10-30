function shock_weight = shock_weight_fun_DFM(model,settings);

shock_weight = zeros(model.n_fac+model.n_y,1);
if settings.est.estimate_shock_weight==1 % if want to estimate the optimal shock weight to maximize impact response
    shock_weight(1:model.n_fac) = model.ABCD.D(settings.est.shock_optimize_var_IRF,1:model.n_fac); % optimal weight is in the same direction as impulse response
    shock_weight = shock_weight/sqrt(shock_weight'*shock_weight); % normalize weight
else
    shock_weight = zeros(n_eps, 1); % manually choose shock
    shock_weight(settings.est.manual_shock_pos) = 1;
end

end