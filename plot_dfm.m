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

dgp_type_plot = 'mp'; % structural shock: either 'G' or 'MP'
estimand_type = 'obsshock'; % structural estimand: either 'obsshock' or 'recursive'
lag_type      = 4; % # of lags to impose in estimation, or NaN (= AIC)
mode_type     = 1; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables)

%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

mode_list   = {'baseline', 'persistent', 'salient', 'more'};
load_mode_dir = mode_list{mode_type};

load_pre = '_results'; % destination to store the results

if isnan(lag_type)
    load_suff = '_aic';
else
    load_suff = num2str(lag_type);
end
load_folder = fullfile(load_pre, load_mode_dir, strcat('lag', load_suff));

if strcmp(dgp_type_plot,'both')
    load(fullfile(load_folder, strcat('dfm_', 'g', '_', estimand_type)))
    results_g = results;
    target_irf_g = DF_model.target_irf;
    load(fullfile(load_folder, strcat('dfm_', 'mp', '_', estimand_type)))
    results_mp = results;
    DF_model.target_irf_g = target_irf_g;
    clear results target_irf_g
    results.coverage_prob = 0.5 * results_g.coverage_prob + 0.5 * results_mp.coverage_prob;
    results.median_length = 0.5 * results_g.median_length + 0.5 * results_mp.median_length;
    results.bias2 = [results_g.bias2;results_mp.bias2];
    results.vce   = [results_g.vce;results_mp.vce];
else
    load(fullfile(load_folder, strcat('dfm_', dgp_type_plot, '_', estimand_type)))
end

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
% Coverage & Width
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
% saveas(the_f,fullfile(plot_folder, strcat('dfm_', dgp_type, '_', estimand_type,'.png')))

%----------------------------------------------------------------
% Bias & Variance
%----------------------------------------------------------------

methods_names_plot  = {'VAR','LP'};
font_size = 18;
exper_plotname = estimand_type;

the_objects = {'bias2','vce'}; % objects to plot
the_titles =  {'Bias','Std'};  % Plot titles/file names

if strcmp(dgp_type_plot,'both')
    the_rms_irf  = sqrt(mean([DF_model.target_irf,DF_model.target_irf_g].^2)); % Root average squared true IRF across horizons
else
    the_rms_irf  = sqrt(mean(DF_model.target_irf.^2)); % Root average squared true IRF across horizons
end

for j=1:length(the_objects)
    
    the_result = sqrt(results.(the_objects{j}));
    the_result = permute(the_result,[3 1 2]);

    % normalized losses
    
    plot_loss(horzs-1, squeeze(quantile(the_result./the_rms_irf, 0.5, 2)), [], ...
        strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size);
%     plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_reltruth', remark_loss_quant)), output_suffix);
    
end