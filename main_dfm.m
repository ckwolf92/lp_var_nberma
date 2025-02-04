%% EXTENDED DFM SIMULATION STUDY
% this version: 10/24/2024

%% HOUSEKEEPING

clc
clear all
close all

path = cd;

addpath(genpath([path '/_auxiliary_functions']));
addpath(genpath([path '/_estim']));
addpath(genpath([path '/_dfm']));
cd([path]);

rng(1, 'twister');

%% SET EXPERIMENT

dgp_type      = 'mp'; % structural shock: either 'G' or 'MP'
estimand_type = 'obsshock'; % structural estimand: either 'obsshock' or 'recursive'
mode_type     = 3; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables), 5 (salient + persistent series)
sample_length = 'medium';  % short (T=100), medium (T=240), long (T=720)
shock_type    = 'arch';  % either 'arch' or 'iid' 

%% SETTINGS

% apply shared settings as well as settings specific to DGP and estimand type

shared;
run(fullfile('_dfm/settings', dgp_type));
run(fullfile('_dfm/settings', estimand_type));
set_mode;

% storage folder for results

save_pre    = '_results'; % destination to store the results
save_folder = fullfile(save_pre, save_mode_dir);

%% ENCOMPASSING DFM MODEL

%----------------------------------------------------------------
% Set up Encompassing Model
%----------------------------------------------------------------

% estimate DFM from dataset

DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar, ...
    DF_model.reorder, DF_model.levels, DF_model.coint_rank);
set_shocks;


% extract and store estimated DFM parameters

DF_model.Phi           = DFM_estimate.Phi;
DF_model.Sigma_eta     = DFM_estimate.Sigma_eta;

DF_model.Lambda        = DFM_estimate.Lambda;
DF_model.delta         = DFM_estimate.delta;
DF_model.sigma_v       = DFM_estimate.sigma_v;

[DF_model.n_y,DF_model.n_fac] = size(DF_model.Lambda);


%----------------------------------------------------------------
% Represent as Model in ABCD Form
%----------------------------------------------------------------

DF_model.ABCD  = ABCD_fun_DFM(DF_model);

%----------------------------------------------------------------
% Shock Weights
%----------------------------------------------------------------

settings.est.shock_weight = shock_weight_fun_DFM(DF_model,settings);

%% PREPARATION

%----------------------------------------------------------------
% Select Individual DGPs from Encompassing Model
%----------------------------------------------------------------

settings.specifications = pick_var_fn(DF_model, settings);

%----------------------------------------------------------------
% Create Placeholders for Results
%----------------------------------------------------------------

% point estimates

estims = zeros(settings.specifications.n_spec, settings.est.n_methods, settings.est.n_IRF, settings.simul.n_mc);

% CIs

ses    = estims;

cis_lower = zeros(settings.specifications.n_spec, settings.est.n_methods, settings.est.n_IRF, 4, settings.simul.n_mc);
cis_upper = cis_lower;

% lags

ps = zeros(settings.specifications.n_spec, settings.simul.n_mc);

%% STRUCTURAL ESTIMANDS AND DGP SUMMARIES

DF_model = dgp_irfs_stats(DF_model, settings, estimand_type);

%% MONTE CARLO ANALYSIS

if settings.sim.num_workers > 0
    poolobj = parpool(settings.sim.num_workers);
end

n_mc = settings.simul.n_mc;
seed = settings.simul.seed;

parfor (i_mc = 1:n_mc, settings.sim.num_workers)

    if mod(i_mc, ceil(n_mc/10)) == 0
        fprintf('%6d%s\n', round(i_mc/n_mc*100), '%');
    end

    %----------------------------------------------------------------
    % Generate Data From Encompassing DFM
    %----------------------------------------------------------------
    
    rng(seed(i_mc), 'twister');

    data_sim_all = generate_data(DF_model,settings);

    %----------------------------------------------------------------
    % Temporary Storage Folders (due to parfor)
    %----------------------------------------------------------------
    
    temp_estims    = NaN(settings.specifications.n_spec, settings.est.n_methods, settings.est.n_IRF);
    temp_ses       = temp_estims;
    temp_cis_lower = zeros(settings.specifications.n_spec, settings.est.n_methods, settings.est.n_IRF, 4);
    temp_cis_upper = temp_cis_lower;
    temp_ps        = zeros(settings.specifications.n_spec, 1);

    %----------------------------------------------------------------
    % Estimation for each DGP
    %----------------------------------------------------------------

    for i_spec = 1:settings.specifications.n_spec
        
        % select data
        
        data_sim = select_data_fn(data_sim_all,settings,i_spec);

        % IRF inference

        i_rep_estims    = nan(settings.est.n_methods, settings.est.n_IRF);
        i_rep_ses       = i_rep_estims;
        i_rep_cis_lower = nan(settings.est.n_methods, settings.est.n_IRF, 4);
        i_rep_cis_upper = i_rep_cis_lower;

        for i_method = 1:settings.est.n_methods

            % Prepare data and get lags
            [data_y,nlags] = estim_prep(settings, data_sim, i_method);

            [i_rep_estims(i_method,:),i_rep_ses(i_method,:),i_cis_dm,i_cis_boot] ...
                = ir_estim(data_y, nlags, 0:settings.est.IRF_hor-1, ...
                settings.est.methods_shared{:}, settings.est.methods{i_method}{:});

            i_rep_cis_lower(i_method,:,1)     = i_cis_dm(1,:);
            i_rep_cis_upper(i_method,:,1)     = i_cis_dm(2,:);
            i_rep_cis_lower(i_method,:,2:end) = i_cis_boot(1,:,:);
            i_rep_cis_upper(i_method,:,2:end) = i_cis_boot(2,:,:);

        end

        temp_estims(i_spec,:,:)      = i_rep_estims;
        temp_ses(i_spec,:,:)         = i_rep_ses;
        temp_cis_lower(i_spec,:,:,:) = i_rep_cis_lower;
        temp_cis_upper(i_spec,:,:,:) = i_rep_cis_upper;
        temp_ps(i_spec)              = nlags;
        
    end

    estims(:,:,:,i_mc)      = temp_estims;
    ses(:,:,:,i_mc)         = temp_ses;
    cis_lower(:,:,:,:,i_mc) = temp_cis_lower;
    cis_upper(:,:,:,:,i_mc) = temp_cis_upper;
    ps(:,i_mc)              = temp_ps;

end

if settings.sim.num_workers > 0
    delete(poolobj); % Stop parallel workers
end

% clear temporary storage

clear i_mc i_spec data_sim_all data_sim i_method IRF_hor n_lags_max est_n_lag n_lags_fix ...
    response_pos normalize_pos with_shock recursive_shock ...
    normalize_with_shock_std_dev recursive_shock_pos data_y n_lags_est nlags ...
    i_cis_boot i_cis_dm i_rep_cis_lower i_rep_cis_upper i_rep_estims i_rep_ses n_mc seed poolobj

%% SUMMARIZE RESULTS

%----------------------------------------------------------------
% Store Results
%----------------------------------------------------------------

results = struct;

results.estims    = estims;
results.ses       = ses;
results.cis_lower = cis_lower;
results.cis_upper = cis_upper;
results.ps        = ps;

clear estims ses cis_lower cis_upper ps

%----------------------------------------------------------------
% Post-Computations
%----------------------------------------------------------------

% coverage

irs_true_reshape = permute(DF_model.target_irf, [2 3 1]);

results.cover_inds = (irs_true_reshape >= results.cis_lower ...
    & irs_true_reshape <= results.cis_upper) + 0; % coverage indicator

results.coverage_prob = mean(results.cover_inds,5); % coverage probability

results.coverage_prob = mean(results.coverage_prob,1);

% median length

results.lengths       = results.cis_upper-results.cis_lower;
results.median_length = median(results.lengths,5);

results.median_length = mean(results.median_length,1);

% bias squared

results.bias2 = (mean(results.estims,4) - irs_true_reshape).^2;

clear irs_true_reshape

% variance

results.vce = var(results.estims,1,4);

%----------------------------------------------------------------
% Save Results
%----------------------------------------------------------------

mkdir(save_folder);
save(fullfile(save_folder, ...
    strcat('dfm_', dgp_type, '_', estimand_type, '_', sample_length, '_', shock_type)), ...
    'DFM_estimate','DF_model','settings','results',...
    'dgp_type','estimand_type','mode_type', 'sample_length', 'shock_type', '-v7.3'); % save results