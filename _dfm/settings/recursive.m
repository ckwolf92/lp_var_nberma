%% SPECIFIC SETTINGS FOR RECURSIVE SHOCK IDENTIFICATION

% indicate estimand

settings.est.with_shock      = 0; % shock is observed and ordered first in data?
settings.est.recursive_shock = 1; % use recursive shock

% indicate position of recursive shock of interest

settings.est.recursive_shock_pos = settings.specifications.random_fixed_pos; % which is recursively defined shock?

% save variable indices

settings.est.methods_shared(2) = {settings.est.IRF_response_var_pos};
settings.est.methods_shared(4) = {settings.est.recursive_shock_pos};