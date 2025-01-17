%% EXTENDED DFM SIMULATION STUDY
% Christian Wolf
% this version: 10/24/2024

%% HOUSEKEEPING

clc
clear all
close all

path = cd;

addpath(genpath([path '/_auxiliary_functions']));
addpath(genpath('_dfm'))
addpath(genpath('_estim'))

rng(1, 'twister');

%% SET EXPERIMENT

dgp_type      = 'mp'; % structural shock: either 'G' or 'MP'
estimand_type = 'recursive'; % structural estimand: either 'obsshock' or 'recursive'
lag_type      = 1; % # of lags to impose in estimation, or NaN (= AIC)
mode_type     = 1; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables)

lags_list = [3 6 12];
n_lags    = length(lags_list);

plot_folder = '_figures/lpvar_illustr';
mkdir(plot_folder)

%% SETTINGS

% apply shared settings as well as settings specific to DGP and estimand type

shared;
run(fullfile('_dfm/settings', dgp_type));
run(fullfile('_dfm/settings', estimand_type));
set_mode;

% adjust settings

settings.specifications.random_n_spec = 1;

settings.simul.n_mc = 1;

settings.simul.T      = 5e5;
settings.simul.T_burn = 1e2;

settings.est.n_methods = length(settings.est.methods);

%% ENCOMPASSING DFM MODEL

%----------------------------------------------------------------
% Set up Encompassing Model
%----------------------------------------------------------------

% estimate DFM from dataset

DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar, ...
    DF_model.reorder, DF_model.levels, DF_model.coint_rank);

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

settings.specifications.var_select = [1 56 120 159 142]; % GDP, unemployment, CPI, GZ bond premium, FFR

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
plot(0:1:IRF_plot-1,scale*estims(4,:,i_lags),'linewidth',4,'linestyle','-.','color',colors.blue)
hold on
plot(0:1:IRF_plot-1,scale*estims(1,:,i_lags),'linewidth',4,'linestyle',':','color',colors.red)
hold on
plot([lags_list(i_lags) lags_list(i_lags)],[-3 2],'Color',settings.colors.grey,'LineWidth',2,'LineStyle','-')
% hold on
set(gcf,'color','w')
title(['$p = \; $' num2str(lags_list(i_lags))],'interpreter','latex','fontsize',18)
xlabel('horizon','interpreter','latex','FontSize',18)
if i_lags == 1
    ylabel('\%','interpreter','latex','FontSize',18)
end
% ylim([-1.5 0.5])
% yticks([-1.5:0.5:0.5])
ylim([-1.25 0.25])
% yticks([-1.25:0.25:0.25])
yticks([-1:0.5:0])
% xlim([0 16])
if i_lags == 2
    legend({'VAR/LP($\infty$)','LP($p$)','VAR($p$)'},'Location','Southeast','fontsize',17,'interpreter','latex')
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.05*pos(3) 0.95*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'lp_vs_var_dfm.eps'))