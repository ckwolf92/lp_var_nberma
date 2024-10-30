%% EXTENDED DFM SIMULATION STUDY: GENERATE FIGURES
% Christian Wolf
% this version: 10/24/2024

%% HOUSEKEEPING

clc
clear all
close all

path = '/Users/christianwolf/Documents/GitHub/lp_var_nberma';

addpath(genpath([path '/_auxiliary_functions']));
addpath(genpath([path '/_estim']));
addpath(genpath([path '/_dfm']));
addpath(genpath([path '/_results']));
cd([path]);

%% LOAD RESULTS

%----------------------------------------------------------------
% Set Experiment
%----------------------------------------------------------------

dgp_type      = 'g'; % structural shock: either 'G' or 'MP'
estimand_type = 'obsshock'; % structural estimand: either 'obsshock' or 'recursive'
lag_type      = 8; % # of lags to impose in estimation, or NaN (= AIC)

%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

load_pre = '_results'; % destination to store the results

if isnan(lag_type)
    load_suff = '_aic';
else
    load_suff = num2str(lag_type);
end
load_folder = fullfile(load_pre, strcat('lag', load_suff));

load(fullfile(load_folder, strcat('dfm_', dgp_type, '_', estimand_type)))

%% GENERATE FIGURES

%----------------------------------------------------------------
% File Name
%----------------------------------------------------------------

plot_pre    = '_figures';
plot_folder = fullfile(plot_pre, strcat('lag', load_suff));
mkdir(plot_folder);

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

proc_names = {'VAR', 'LP'};
procs = [1 1; % first index: inference procedure; second index: type of confidence interval
         2 1];
line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
line_specs = {'-', '-.'};

numproc = length(proc_names); % number of inference procedures

ylim_cover = [0 1]; % y-limits for coverage prob plot
yticks_length = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1', '10'}; % y-tick labels for median length plot

horzs   = settings.est.IRF_select;
numhorz = length(horzs); % no. of estimated impulse response horizons

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

% size

plotwidth = 0.39;
gapsize = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% plot figure

d = 1;

the_f = figure;

% coverage probability
subplot(1,2,1)
pos = get(gca, 'Position');
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on;
for j=1:numproc
    plot(horzs, squeeze(results.coverage_prob(d,procs(j,1),:,procs(j,2))), ...
        line_specs{j}, 'Color', line_colors(j,:),'LineWidth',5);
end
plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 3.5); % Nominal confidence level
hold off;
xlim([min(horzs) max(horzs)])
xlabel('horizon','interpreter','latex');
ylim(ylim_cover);
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
hold on;
for j=1:numproc
    plot(horzs, squeeze(log10(results.median_length(d,procs(j,1),:,procs(j,2)))), ...
        line_specs{j}, 'Color', line_colors(j,:),'LineWidth',5);
end
hold off;
xlim([min(horzs) max(horzs)])
xlabel('horizon','interpreter','latex');
ylim([min(yticks_length) max(yticks_length)]);
yticks(yticks_length);
yticklabels(yticklabels_length);
title('median length, log scale','interpreter','latex');
legend(proc_names, 'Location', 'SouthEast','interpreter','latex');
grid on;

% size
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
saveas(the_f,fullfile(plot_folder, strcat('dfm_', dgp_type, '_', estimand_type,'.png')))