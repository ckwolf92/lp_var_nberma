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
estimand_type = 'recursive'; % structural estimand: either 'obsshock' or 'recursive'
mode_type     = 6; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables), 5 (salient + persistent series),
                   % 6 (combine 3 & 5)

%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

mode_list   = {'baseline', 'persistent', 'salient', 'more', 'persistent_salient', 'salient_all'};

load_pre = '_results_small_2025_0115';

covg_cutoff = 0.8; % cut-off for coverage indicator plots

if mode_type <= 5

    load_mode_dir = mode_list{mode_type};
    load_folder = fullfile(load_pre, load_mode_dir);
    
    load_results_aux

elseif mode_type == 6

    load_mode_dir = mode_list{3};
    load_folder = fullfile(load_pre, load_mode_dir);
    
    load_results_aux
    results_1 = results;

    load_mode_dir = mode_list{5};
    load_folder = fullfile(load_pre, load_mode_dir);
    
    load_results_aux
    results_2 = results;

    clear results

    results.target_irf = [results_1.target_irf,results_2.target_irf];

    results.bias2 = [results_1.bias2;results_2.bias2];
    results.vce   = [results_1.vce;results_2.vce];
    results.mse   = 0.5 * results.bias2 + 0.5 * results.vce;

    results.p = [results_1.p;results_2.p];

    results.coverage_prob  = [results_1.coverage_prob;results_2.coverage_prob];
    results.coverage_avg   = mean(results.coverage_prob,1);
    results.coverage_indic = mean(results.coverage_prob >= covg_cutoff,1);

    results.median_length = [results_1.median_length;results_2.median_length];  
    results.median_avg    = mean(results.median_length,1); 

    results.M = [results_1.M;results_2.M];

    mode_type     = 6;
    load_mode_dir = mode_list{mode_type};

end

clear covg_cutoff

% procedure names

settings.est.names = {'VAR$_{AIC}$', 'VAR$_{4}$', 'VAR$_{8}$', 'LP$_{AIC}$', 'LP$_{4}$', ...
    'VAR$_{AIC, small}$', 'LP$_{AIC, small}$', 'VAR$_{b, AIC}$', 'LP$_{b, AIC}$', 'VAR$_{b, 4}$', 'LP$_{b, 4}$'};

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
                colors.red; colors.blue; ...
                colors.dred; colors.dblue];

line_specs = {'-', ':', ':', '-', ':', ':', ':', '-.', '-.', '-.', '-.'};

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
    if mode_type == 5
        ylim([0 1])
        yticks([0:0.2:1])
    elseif mode_type == 6 && strcmp(dgp_type_plot,'mp') && j == 2
        ylim([0 1])
        yticks([0:0.2:1])
    end
    if strcmp(estimand_type,'obsshock')
        xlim([min(horzs) max(horzs)])
    elseif strcmp(estimand_type,'recursive')
        xlim([min(horzs)+1 max(horzs)])
    end
    xlabel('horizon','interpreter','latex');
    if j == 1
    else
        if mode_type == 5
            legend(proc_estim, 'Location', 'southeast', 'NumColumns', 1, 'interpreter', 'latex', 'FontSize', 18);
        else
            legend(proc_estim, 'Location', 'east', 'NumColumns', 1, 'interpreter', 'latex', 'FontSize', 18);
        end
    end
%     legend(proc_estim, 'Location', 'eastoutside', 'NumColumns', 1, 'interpreter', 'latex', 'FontSize', 18);
    grid on

    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1.2*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    exportgraphics(gcf, fullfile(plot_folder, [the_objects{j} '.eps']))

end

%----------------------------------------------------------------
% Inference: Coverage-Width Plots
%----------------------------------------------------------------

% preparations

if strcmp(estimand_type,'obsshock')
    proc_inference_indic_1 = [1 8 2 6 4 9 5 7];
    proc_inference_indic_2 = [1 1 2 6 4 4 5 7; 1 2 1 1 1 4 1 1];
elseif strcmp(estimand_type,'recursive')
    proc_inference_indic_1 = [1 2 8 10 4 5 9 11];
    proc_inference_indic_2 = [1 2 1 2 4 5 4 5; 1 1 2 2 1 1 4 4];
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

% average coverage and width

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
    plot(horzs, squeeze(results.coverage_avg(1,proc_inference_indic_2(1,i_proc),:,proc_inference_indic_2(2,i_proc))), ...
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
    plot(horzs, squeeze(log10(results.median_avg(1,proc_inference_indic_2(1,i_proc),:,proc_inference_indic_2(2,i_proc)))), ...
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
if mode_type == 3
    legend(proc_inference, 'Location', 'NorthEast','NumColumns', 2, 'interpreter','latex');
else
    legend(proc_inference, 'Location', 'SouthEast','NumColumns', 2, 'interpreter','latex');
end
grid on;

% size

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'covgwidth.eps'))

% coverage indicator

figure(length(the_objects)+2)

pos = get(gca, 'Position');
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Position', pos)
hold on
for i_proc = 1:numproc_inference
    plot(horzs, squeeze(results.coverage_indic(1,proc_inference_indic_2(1,i_proc),:,proc_inference_indic_2(2,i_proc))), ...
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
ylim([0 1]);
legend(proc_inference, 'Location', 'SouthEast','NumColumns', 2, 'interpreter','latex');
grid on

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'covgindic.eps'))

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

figure(length(the_objects)+3)

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
legend({'$p = 2$','$p = 4$','$p = 8$'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.8*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'misspec.eps'))