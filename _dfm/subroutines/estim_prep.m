function [data_y,nlags] = estim_prep(settings,data_sim);

% unpack settings

n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
n_lags_fix = settings.est.n_lags_fix;

with_shock      = settings.est.with_shock;

% collect data

if with_shock == 1 % observe shock: w_t = (shock, \bar{w}_t)
    data_y = [data_sim.data_shock,data_sim.data_y]; % Warning: correspond to w_t in our paper
else % recursive: w_t = \bar{w}_t
    data_y = data_sim.data_y;
end

% estimate lag length

n_lags_est = ic_var(data_y,n_lags_max,1);

% set lag length

if est_n_lag == 0 % fix lag order
    nlags = n_lags_fix;
else % use estimated lag order
    nlags = n_lags_est; 
end

end