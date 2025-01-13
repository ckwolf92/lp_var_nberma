%% EXTENDED DFM SIMULATION STUDY: GENERATE FIGURES
% this version: 01/12/2025

%% HOUSEKEEPING

clc
clear all
close all

path = '/Users/christianwolf/Documents/GitHub/lp_var_nberma';
path_results = '/Users/christianwolf/Dropbox/Research/lp_var_nberma/codes';

addpath(genpath([path '/_auxiliary_functions']));
addpath(genpath([path '/_estim']));
addpath(genpath([path '/_dfm']));
addpath(genpath(path_results))
cd([path]);

%% LOAD RESULTS

%----------------------------------------------------------------
% Set Experiment
%----------------------------------------------------------------

dgp_type_plot = 'both'; % structural shock: either 'g', or 'mp', or 'both'
estimand_type = 'obsshock'; % structural estimand: either 'obsshock' or 'recursive'
mode_type     = 5; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables), 5 (salient + persistent series)

%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

mode_list   = {'baseline', 'persistent', 'salient', 'more', 'persistent_salient'};
load_mode_dir = mode_list{mode_type};

load_pre = '_results_small_2025_0113';
load_folder = fullfile(load_pre, load_mode_dir);

% this should probably be cleaned up

if strcmp(dgp_type_plot,'both')

    load(fullfile(load_folder, strcat('dfm_', 'g', '_', estimand_type)))
    results_g = results;
    results_g.coverage_prob = mean(results_g.coverage_prob,1);
    results_g.target_irf = DF_model.target_irf;
    results_g.M = DF_model.M;
    clear results DF_model DFM_estimate dgp_type

    load(fullfile(load_folder, strcat('dfm_', 'mp', '_', estimand_type)))
    results_mp = results;
    results_mp.coverage_prob = mean(results_mp.coverage_prob,1);
    results_mp.target_irf = DF_model.target_irf;
    results_mp.M = DF_model.M;
    clear results DF_model DFM_estimate dgp_type

    results.target_irf = [results_g.target_irf,results_mp.target_irf];
    results.bias2 = [results_g.bias2;results_mp.bias2];
    results.vce   = [results_g.vce;results_mp.vce];
    results.mse   = 0.5 * results.bias2 + 0.5 * results.vce;
    results.coverage_prob = 0.5 * results_g.coverage_prob + 0.5 * results_mp.coverage_prob;
    results.median_length = 0.5 * results_g.median_length + 0.5 * results_mp.median_length;  
    results.M = [results_g.M;results_mp.M];
    
else

    load(fullfile(load_folder, strcat('dfm_', dgp_type_plot, '_', estimand_type)))
    results.target_irf = DF_model.target_irf;
    results.mse        = results.bias2 + results.vce;
    results.coverage_prob = mean(results.coverage_prob,1);
    results.M = DF_model.M;
    clear DF_model DFM_estimate dgp_type

end

settings.est.names = {'VAR$_{AIC}$', 'VAR$_{4}$', 'VAR$_{8}$', 'LP$_{AIC}$', 'LP$_{4}$', ...
    'VAR$_{4, small}$', 'LP$_{4, small}$', 'VAR$_{b, AIC}$', 'LP$_{b, AIC}$'};

%% GENERATE FIGURES

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

% figures folder

plot_pre    = '_figures';
plot_folder = fullfile(plot_pre, load_mode_dir);
plot_folder = fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type));
mkdir(plot_folder)

% plot objects

horzs = settings.est.IRF_select - 1;

% lines

colors.blue  = [102/255 178/255 255/255];
colors.dblue = 0.7 * colors.blue;
colors.lblue = 0.4 * colors.blue + 0.6 * [1 1 1];
colors.red   = [204/255 0/255 0/255];
colors.dred  = 0.7 * colors.red;
colors.lred  = 0.4 * colors.red + 0.6 * [1 1 1];

line_colors = [colors.red; colors.dred; colors.dred; ...
                colors.blue; colors.dblue; ...
                colors.lred; colors.lblue; ...
                colors.red; colors.blue];

line_specs = {'-', ':', ':', '-', ':', ':', ':', '-.', '-.'};

line_width = 5 * ones(length(line_specs),1);

%----------------------------------------------------------------
% Estimation: Bias-Variance-MSE Plots
%----------------------------------------------------------------

% preparations

if strcmp(estimand_type,'obsshock')
    proc_estim_indic = [1 2 6 4 5 7];
elseif strcmp(estimand_type,'recursive')
    proc_estim_indic = [1 2 4 5];
end

proc_estim    = settings.est.names(proc_estim_indic);
numproc_estim = length(proc_estim);

the_rms_irf = sqrt(mean(results.target_irf.^2)); % root average squared true IRF across horizons
the_objects = {'bias2','vce','mse'}; % objects to plot

% figures

for j = 1:length(the_objects)

    the_result = sqrt(results.(the_objects{j}));
    the_result = permute(the_result,[3 1 2]);
    the_result = squeeze(quantile(the_result./the_rms_irf, 0.5, 2));

    figure(j)

    pos = get(gca, 'Position');
    set(gca,'FontSize',18)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Position', pos)
    hold on
    for i_proc = 1:numproc_estim
        plot(horzs, the_result(:,proc_estim_indic(i_proc)), line_specs{proc_estim_indic(i_proc)}, ...
            'Color', line_colors(proc_estim_indic(i_proc),:), 'LineWidth', line_width(proc_estim_indic(i_proc)));
        hold on
    end
    if strcmp(estimand_type,'obsshock')
        xlim([min(horzs) max(horzs)])
    elseif strcmp(estimand_type,'recursive')
        xlim([min(horzs)+1 max(horzs)])
    end
    xlabel('horizon','interpreter','latex');
    legend(proc_estim, 'Location', 'eastoutside', 'NumColumns', 1, 'interpreter', 'latex', 'FontSize', 18);
    grid on

    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    exportgraphics(gcf, fullfile(plot_folder, [the_objects{j} '.eps']))

end

%----------------------------------------------------------------
% Inference: Coverage-Width Plot
%----------------------------------------------------------------

% preparations

if strcmp(estimand_type,'obsshock')
    proc_inference_indic_1 = [1 8 2 6 4 9 5 7];
    proc_inference_indic_2 = [1 1 2 6 4 4 5 7; 1 2 1 1 1 3 1 1];
elseif strcmp(estimand_type,'recursive')
    proc_inference_indic_1 = [1 8 2 4 9 5];
    proc_inference_indic_2 = [1 1 2 4 4 5; 1 2 1 1 3 1];
end

proc_inference    = settings.est.names(proc_inference_indic_1);
numproc_inference = length(proc_inference);

% plot settings

plotwidth = 0.39;
gapsize = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

yticks_length = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

% figure

figure(length(the_objects)+1)

% coverage probability

subplot(1,2,1)
pos = get(gca, 'Position');
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
for i_proc = 1:numproc_inference
    plot(horzs, squeeze(results.coverage_prob(1,proc_inference_indic_2(1,i_proc),:,proc_inference_indic_2(2,i_proc))), ...
        line_specs{proc_inference_indic_1(i_proc)}, 'Color', line_colors(proc_inference_indic_1(i_proc),:),...
        'LineWidth',line_width(proc_inference_indic_1(i_proc)));
    hold on
end
plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 3.5); % Nominal confidence level
if strcmp(estimand_type,'obsshock')
    xlim([min(horzs) max(horzs)])
elseif strcmp(estimand_type,'recursive')
    xlim([min(horzs)+1 max(horzs)])
end
xlabel('horizon','interpreter','latex');
ylim([0 1]);
title('coverage probability','interpreter','latex');
grid on

% median length

subplot(1,2,2)
pos = get(gca, 'Position');
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
for i_proc = 1:numproc_inference
    plot(horzs, squeeze(log10(results.median_length(1,proc_inference_indic_2(1,i_proc),:,proc_inference_indic_2(2,i_proc)))), ...
        line_specs{proc_inference_indic_1(i_proc)}, 'Color', line_colors(proc_inference_indic_1(i_proc),:),...
        'LineWidth',line_width(proc_inference_indic_1(i_proc)));
    hold on
end
if strcmp(estimand_type,'obsshock')
    xlim([min(horzs) max(horzs)])
elseif strcmp(estimand_type,'recursive')
    xlim([min(horzs)+1 max(horzs)])
end
xlabel('horizon','interpreter','latex');
ylim([min(yticks_length) max(yticks_length)]);
yticks(yticks_length);
yticklabels(yticklabels_length);
title('median length, log scale','interpreter','latex');
legend(proc_inference, 'Location', 'NorthEast','NumColumns', 2, 'interpreter','latex');
grid on;

% size

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'covgwidth.eps'))

%----------------------------------------------------------------
% Degree of Mis-Specification
%----------------------------------------------------------------

% preparations

M_draws = results.M;

n_lags   = size(M_draws,2);
n_kernel = 1001;

% colors

colors.purple  = [160/255 32/255 240/255];
colors.lpurple = 0.5 * colors.purple + 0.5 * [1 1 1];
colors.blue    = [116/255 158/255 178/255];
colors.lblue   = 0.5 * colors.blue + 0.5 * [1 1 1];
colors.orange  = [204/255 102/255 0/255];
colors.lorange = 0.5 * colors.orange + 0.5 * [1 1 1];

colors.base  = [colors.blue; colors.purple; colors.orange];
colors.light = [colors.lblue; colors.lpurple; colors.lorange];

% get Kernel densities

M_lb = 0;
M_ub = 1.2 * max(max(M_draws));

M_density = NaN(2,n_kernel,n_lags);
M_density(1,:,:) = repmat(linspace(M_lb,M_ub,n_kernel),1,1,n_lags);

for i_lag = 1:n_lags

    M_density(2,:,i_lag) = kernel(M_density(1,:,i_lag),M_draws(:,i_lag),0.4);
    M_density(2,:,i_lag) = M_density(2,:,i_lag) * (M_density(1,2,i_lag) - M_density(1,1,i_lag));

end

% figure

figure(length(the_objects)+2)

pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',17)
set(gca,'TickLabelInterpreter','latex')
hold on
for i_lag = 1:n_lags
    jbfill(M_density(1,:,i_lag),0*M_density(2,:,i_lag),M_density(2,:,i_lag),...
        colors.light(i_lag,:),colors.light(i_lag,:),0,0.5);
    hold on
end
for i_lag = 1:n_lags
    plot(M_density(1,:,i_lag),M_density(2,:,i_lag),'linewidth',3,'linestyle','-','color',colors.base(i_lag,:));
    hold on
end
set(gcf,'color','w')
xlabel('Degree of Misspecification','interpreter','latex','FontSize',20)
xlim([M_lb M_ub]);
legend({'$p = 1$','$p = 4$','$p = 8$'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.8*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'misspec.eps'))