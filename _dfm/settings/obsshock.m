%% SPECIFIC SETTINGS FOR OBSERVED SHOCK IDENTIFICATION

% indicate estimand

settings.est.with_shock      = 1; % shock is observed and ordered first in data?
settings.est.recursive_shock = 0; % use recursive shock

% indicate normalization scheme

settings.est.normalize_with_shock_std_dev = 1; % normalize IRF with one unit of shock (shock std-dev is 1)? Otherwise using normalization variable

% save variable indices

settings.est.methods_shared(2) = {settings.est.IRF_response_var_pos + 1}; % +1 since shock is included
settings.est.methods_shared(4) = {1};