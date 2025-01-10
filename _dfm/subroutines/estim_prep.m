function [data_y,nlags] = estim_prep(settings,data_sim);

% unpack settings

n_lags_max = settings.est.n_lags_max;
with_shock = settings.est.with_shock;

% collect data

if with_shock == 1 % observe shock: w_t = (shock, \bar{w}_t)
    data_y = [data_sim.data_shock,data_sim.data_y]; % Warning: correspond to w_t in our paper
else % recursive: w_t = \bar{w}_t
    data_y = data_sim.data_y;
end

% For "small", only include shock and response variable

if strcmp(settings.specifications.system_type, 'small')  
    data_y = [data_sim.data_shock, data_sim.data_y(:, settings.est.IRF_response_var_pos)];
end

% estimate lag length

n_lags_est = ic_var(data_y,n_lags_max,1);

% set lag length

nlags                         = settings.est.n_lags_fix;  % Fixed
nlags(settings.est.est_n_lag) = n_lags_est;               % Set to estimated


end