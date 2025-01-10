function [data_y,nlags] = estim_prep(settings,data_sim, i_method);

% unpack settings

n_lags_max = settings.est.n_lags_max;
with_shock = settings.est.with_shock;

% collect data

switch settings.est.system_type{i_method}
    case 'big'
        if with_shock == 1 % observe shock: w_t = (shock, \bar{w}_t)
            data_y = [data_sim.data_shock,data_sim.data_y]; % Warning: correspond to w_t in our paper
        else % recursive: w_t = \bar{w}_t
            data_y = data_sim.data_y;
        end

    case 'small'
        % For "small", only include shock and response variable
        data_y = [data_sim.data_shock, data_sim.data_y(:, settings.est.IRF_response_var_pos)];
end


if settings.est.est_n_lag(i_method)
    % estimate lag length
    nlags = ic_var(data_y,n_lags_max,1);
else
    % Set to fix
    nlags = settings.est.n_lags_fix(i_method);
end


end