% Function for estimating IRFs using a Bayesian VAR approach
function [irs, cis_post, GLP_hyper_est] = bvarglp_ir_estim(Y, p, horzs, RW, varargin)
%% Parse inputs
ip = inputParser;

% Required inputs
addRequired(ip, 'Y', @isnumeric);
% T x n     data vector
addRequired(ip, 'p', @isnumeric);
% 1 x 1     lag length (not counting any augmentation)
addRequired(ip, 'horzs', @isnumeric);
% 1 x H     impulse response horizons of interest
addRequired(ip, 'RW', @islogical);
% 1 x 1     Random Walk prior? 


% Optional inputs
addParameter(ip, 'resp_ind', 1, @isnumeric);
% Index of response variable of interest (default: first variable)
addParameter(ip, 'innov_ind', 1, @isnumeric);
% Index of innovation of interest (default: first innovation)
addParameter(ip, 'alpha', 0.05, @isnumeric);
% Significance level (default: 0.05)
addParameter(ip, 'ndraw', 500, @isnumeric);
% Posterior draws (default: 500)

parse(ip, Y, p, horzs, RW, varargin{:});


%% Preliminaries

n_Y = size(Y,2);

%% Estimate BVAR

if RW % Random walk prior
    r = bvarGLP(Y, p, 'MNpsi', 0, 'Fcast', 0);
else % White noise prior
    r = bvarGLP(Y, p, 'MNpsi', 0, 'Fcast', 0, 'sur', 0, 'noc', 0, 'posi', 1:n_Y);
end

%% Posterior draws of parameters and IRFs

if ip.Results.ndraw == 0 % only use posterior means
    beta_draws  = r.postmax.betahat;
    sigma_draws = r.postmax.sigmahat;

else % draw from posterior

    beta_draws  = nan(1+n_Y*p, n_Y, ip.Results.ndraw);
    sigma_draws = nan(n_Y, n_Y, ip.Results.ndraw);
    
    for j=1:ip.Results.ndraw
        [beta_draws(:,:,j),sigma_draws(:,:,j)] = post_draw(r.postmax.betahat,r.postmax.Sinv,r.postmax.cholZZinv,r.postmax.T);
    end
end

% IRFs
IRF_draws = nan(n_Y, length(horzs), ip.Results.ndraw);

for j=1:ip.Results.ndraw
    G                = chol(sigma_draws(:,:,j),'lower');
    ShockVector      = G(:,ip.Results.innov_ind);
    IRF_draws(:,:,j) = var_ir(squeeze(beta_draws(2:end,:,j))', ShockVector, horzs);
end

% hyperparameters
GLP_hyper_est = [r.postmax.lambda;r.postmax.theta;r.postmax.miu];

% Organize output
irs_all      = mean(IRF_draws, 3); % posterior mean of IRF draws (or just single IRF if ndraw=0)
cis_all_post = quantile(IRF_draws, [ip.Results.alpha/2 1-ip.Results.alpha/2], 3);  % Equal-tailed
irs          = irs_all(ip.Results.resp_ind, :);
cis_post     = cis_all_post(ip.Results.resp_ind, :);

end