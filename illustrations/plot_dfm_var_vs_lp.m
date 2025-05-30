%% VAR VS. LP DFM ILLUSTRATION
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 04/29/2025

%% HOUSEKEEPING

clear
clc
close all

path = cd;

addpath(genpath('../_auxiliary_functions'));
addpath(genpath('../_dfm'))
addpath(genpath('../_estim'))

rng(1, 'twister');

%% SET EXPERIMENT

dgp_type      = 'mp'; % structural shock: either 'G' or 'MP'
estimand_type = 'recursive'; % structural estimand: either 'obsshock' or 'recursive'
lag_type      = 1; % # of lags to impose in estimation, or NaN (= AIC)
mode_type     = 3; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables)
sample_length = 'medium';  % short (T=100), medium (T=240), long (T=720)
shock_type    = 'iid'; % 'iid' or 'arch'

lags_list = [2 6 12];
n_lags    = length(lags_list);

plot_folder = '_figures/lpvar_illustr';
mkdir(plot_folder)

%% SETTINGS

% shared

shared;

settings.est = rmfield(settings.est,'methods');

% VAR

[settings.est.methods{1}] = deal({'resp_ind' , [], ...
    'innov_ind', [], ...
    'bootstrap', [],...
    'estimator', 'var', ...
    'shrinkage', false});

% LP 

[settings.est.methods{2}] = deal({'resp_ind'  , [], ...
    'innov_ind' , [], ...
    'bootstrap', [],...
    'estimator' , 'lp' , ...
    'shrinkage' , false});  % LP

% lags

settings.est.n_lags_fix          = [4; 4];
settings.est.system_type         = repmat({'big'}, 2, 1);

settings.est.est_n_lag  = isnan(settings.est.n_lags_fix);  % indicator if lags are estimated
settings.est.n_lags_max = 10;
settings.est.n_methods  = length(settings.est.methods); % number of estimation methods

% rest

run(fullfile('../_dfm/settings', dgp_type));
run(fullfile('../_dfm/settings', estimand_type));
set_mode;

% adjust settings

settings.specifications.random_n_spec = 1;

settings.simul.n_mc = 1;

settings.simul.T      = 5e5;
settings.simul.T_burn = 1e2;
% settings.simul.T      = 200;
% settings.simul.T_burn = 100;

settings.est.n_methods = length(settings.est.methods);

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

%% STRUCTURAL ESTIMANDS AND DGP SUMMARIES

settings.specifications = pick_var_fn(DF_model, settings);

settings.specifications.var_select = [56 1 121 151 142]; % unemployment, GDP, CPI, spread, FFR

DF_model = dgp_irfs_stats(DF_model, settings, estimand_type);

%% MONTE CARLO ANALYSIS

%----------------------------------------------------------------
% Generate Data
%----------------------------------------------------------------

rng(settings.simul.seed(1), 'twister');

% simulate DFM

data_sim_all = generate_data(DF_model,settings);

% select specification

data_sim = select_data_fn(data_sim_all,settings,1);

%----------------------------------------------------------------
% Estimation
%----------------------------------------------------------------

% estimation preparation

data_y = estim_prep(settings,data_sim,1);

% IRF inference

estims    = nan(settings.est.n_methods, settings.est.n_IRF, n_lags);
ses       = estims;

for i_lags = 1:n_lags

    disp(i_lags)
    nlags = lags_list(i_lags);
        
    for i_method = 1:settings.est.n_methods
    
        [estims(i_method,:,i_lags),ses(i_method,:,i_lags)] ...
            = ir_estim(data_y, nlags, 0:settings.est.IRF_hor-1, settings.est.methods_shared{:}, settings.est.methods{i_method}{:});
    
    end

end

%% PLOT RESULTS

%----------------------------------------------------------------
% Color Preparation
%----------------------------------------------------------------

colors.blue  = [102/255 178/255 255/255];
colors.dblue = 0.7 * colors.blue;
colors.lblue = 0.4 * colors.blue + 0.6 * [1 1 1];
colors.red   = [204/255 0/255 0/255];
colors.dred  = 0.7 * colors.red;
colors.lred  = 0.4 * colors.red + 0.6 * [1 1 1];

colors.grey = [170/255 170/255 170/255];

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

IRF_plot = 21;

scale = 1/max(abs(DF_model.target_irf));

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

plotwidth = 0.267;
gapsize = 0.045;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

figure(1)

for i_lags = 1:n_lags

subplot(1,3,i_lags)
pos = get(gca, 'Position');
pos(1) = left_pos(i_lags);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot-1,scale*DF_model.target_irf(1:IRF_plot),'linewidth',4,'linestyle','-','color',colors.grey)
hold on
plot(0:1:IRF_plot-1,scale*estims(2,:,i_lags),'linewidth',4,'linestyle','-.','color',colors.blue)
hold on
plot(0:1:IRF_plot-1,scale*estims(1,:,i_lags),'linewidth',4,'linestyle',':','color',colors.red)
hold on
plot([lags_list(i_lags) lags_list(i_lags)],[-3 2],'Color',colors.grey,'LineWidth',2,'LineStyle','-')
% hold on
set(gcf,'color','w')
title(['$p = \; $' num2str(lags_list(i_lags))],'interpreter','latex','fontsize',21)
xlabel('horizon','interpreter','latex','FontSize',18)
if i_lags == 1
    ylabel('\%','interpreter','latex','FontSize',18)
end
% ylim([-1.5 0.5])
% yticks([-1.5:0.5:0.5])
ylim([-0.5 1.5])
% % yticks([-1.25:0.25:0.25])
yticks([-0.5:0.5:1.5])
% xlim([0 16])
if i_lags == 2
    legend({'VAR/LP($\infty$)','LP($p$)','VAR($p$)'},'Location','Northeast','fontsize',17,'interpreter','latex')
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.05*pos(3) 0.95*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'lp_vs_var_dfm.eps'))