%% EXTENDED DFM SIMULATION STUDY: GENERATE FIGURES
% Christian Wolf
% this version: 10/24/2024

%% HOUSEKEEPING

clc
clear all
close all

path = cd;

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
mode_type     = 5; % robustness check mode:
                   % 1 (baseline), 2 (persistent), 3 (salient series),
                   % 4 (more observables)
sample_length = 'medium';
shock_type    = 'arch';
load_pre = '_results_small_2025_0225'; % destination to store the results
plot_pre    = '_figures_2025_0225';


%%
for dgp = {'mp', 'g'}
    for estimand_type = {'obsshock', 'recursive'}
        for mode_type = [3 5]
            for shock_type = {'iid', 'arch'}
                plot_dfm_summary_main(dgp{1}, estimand_type{1},...
                    mode_type, 'medium', shock_type{1}, load_pre, plot_pre)
            end
        end
    end
end

%% MAIN PLOTTING FUNCTION
function plot_dfm_summary_main(dgp_type_plot, estimand_type, mode_type, sample_length, shock_type, load_pre, plot_pre)
%----------------------------------------------------------------
% Load
%----------------------------------------------------------------

mode_list   = {'baseline', 'persistent', 'salient', 'more', 'persistent_salient'};
load_mode_dir = mode_list{mode_type};


load_folder = fullfile(load_pre, load_mode_dir);

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
    results.bias2         = [results_g.bias2;results_mp.bias2];
    results.vce           = [results_g.vce;results_mp.vce];
else
   temp = load(fullfile(load_folder, strcat('dfm_', dgp_type_plot, '_', estimand_type, '_', sample_length, '_', shock_type)));
   DF_model = temp.DF_model;
   DFM_estimate = temp.DFM_estimate;
   dgp_type = temp.dgp_type;
   estimand_type = temp.estimand_type;
   mode_type = temp.mode_type;
   results = temp.results;
   settings = temp.settings;
   shock_type = temp.shock_type;
end



%% GENERATE FIGURES

%----------------------------------------------------------------
% File Name
%----------------------------------------------------------------

plot_folder = fullfile(plot_pre, mode_list{mode_type});
plot_folder = fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type));
mkdir(plot_folder)

%% COVERAGE AND WIDTH

close all

% Lag 1
procs = [1 1; % first index: inference procedure; second index: type of confidence interval
         1 2;
         5 1;
         5 4];

plot_coveragewidth(procs, settings, results, estimand_type)
exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type, '_', shock_type, '.pdf')), "Append", false)



% Lag 4
procs = [2 1; % first index: inference procedure; second index: type of confidence interval
         2 2;
         6 1;
         6 4];
plot_coveragewidth(procs, settings, results, estimand_type)
exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type,'_', shock_type, '.pdf')), "Append", true)


% Homosk vs heterosk
procs = [1 2;
         4 2;
         5 4;
         7 4];

plot_coveragewidth(procs, settings, results, estimand_type)
exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type, '_', shock_type,'.pdf')), "Append", true)


% Alternative estimators
procs = [1 2
         5 4;
         8 1;
         9 1;
         10 1];

plot_coveragewidth(procs, settings, results, estimand_type)
exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type, '_', shock_type,'.pdf')), "Append", true)




if strcmp(estimand_type, 'obsshock')
%    % Big vs small systems
    procs = [1 2;
        5 4;
        11 2;
        12 4];
    plot_coveragewidth(procs, settings, results, estimand_type)
    exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type,'_', shock_type, '.pdf')), "Append", true)
end

% VAR comparison
procs = [1 2; 
         2 2;
         3 2
         5 4];
plot_coveragewidth(procs, settings, results, estimand_type)
    exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type, '_', shock_type,'.pdf')), "Append", true)


%% Bias & Variance

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
    horzs   = 0:settings.est.IRF_hor-1;

    plot_loss(horzs, squeeze(quantile(the_result(:,:,[2 6])./the_rms_irf, 0.5, 2)), [], ...
        strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size);    
    title(the_titles{j})
    set(gcf, 'Units','inches', 'Position', [2 2 6 3])
    set(gca, 'FontSize', 11)
    

    exportgraphics(gcf, fullfile(plot_folder, strcat(dgp_type_plot, '_', estimand_type, '_', shock_type,'.pdf')), "Append", true)
    
end
end
%% PLOT FUNCTIONS

% size
function plot_coveragewidth(procs, settings, results, estimand_type)
close all
% Colors for each estimator
line_colors = [153,0,13;      % Dark red
               250 77 86;     % Medium red
              252,187,161;    % Light red
              252,187,161;    % Light red
               31,120,180;    % Dark blue
               166,206,227;   % Light blue
               166,206,227;   % Light blue
               255,127,0;      % Orange
               77,175,74;      % Green
                [77,175,74]/2 + [255 255 255]/2;  % Light green
                255,127,0;
               166,86,40]./255;
line_specs = {'-', '-.', '-.', '-.'};  % Pattern for bootstrap
plotwidth     = 0.39;
gapsize       = 0.1;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos      = [gapsize_edges, gapsize_edges + gapsize + plotwidth];
LW            = 5;
ylim_cover         = [0 1]; % y-limits for coverage prob plot
yticks_length      = -3:1:1; % y-ticks for median length plot (log10 scale)
yticklabels_length = {'0.001', '0.01', '0.1', '1'}; % y-tick labels for median length plot


numproc = length(procs); % number of inference procedures
proc_names = {'VAR AIC', 'VAR(4)', 'VAR(8)', 'VAR AIC homosk', 'LP AIC', 'LP(4)',...
              'LP AIC homosk', 'BVAR', 'SLP', 'SLP undersm', 'VAR AIC, small', 'LP AIC, small'};

covg         = squeeze(mean(min(results.coverage_prob, .9), 1));
covg_frac_70 = squeeze(mean(results.coverage_prob>=.7,1));
median_length = squeeze(mean(results.median_length, 1));

% plot figure

horzs   = 0:settings.est.IRF_hor-1;
numhorz = length(horzs); % no. of estimated impulse response horizons

%%
figure('Visible', 'off')
% coverage probability
subplot(1,3,1)
for j=1:numproc
    ind_esttype  = procs(j, 1);
    ind_citype   = procs(j, 2);
    plot(horzs, squeeze(covg(procs(j,1),:,procs(j,2))), ...
        line_specs{ind_citype}, 'Color', line_colors(ind_esttype,:),'LineWidth',LW);
    hold on
end
plot(horzs, (1-settings.est.alpha) * ones(1,length(horzs)), ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 3.5); % Nominal confidence level
hold off;
if strcmp(estimand_type, 'obsshock')
    xlim([min(horzs) max(horzs)])
else
    xlim([min(horzs)+1 max(horzs)])
end
xlabel('horizon','interpreter','latex');
ylim(ylim_cover);
title('Coverage shortfall','interpreter','latex');
grid on
set(gca, 'TickLabelInterpreter', 'latex')



% coverage frac >=70%
subplot(1,3,2)
for j=1:numproc
    ind_esttype  = procs(j, 1);
    ind_citype   = procs(j, 2);
    plot(horzs, squeeze(covg_frac_70(procs(j,1),:,procs(j,2))), ...
        line_specs{ind_citype}, 'Color', line_colors(ind_esttype,:),'LineWidth',LW);
    hold on
end
hold off;
if strcmp(estimand_type, 'obsshock')
    xlim([min(horzs) max(horzs)])
else
    xlim([min(horzs)+1 max(horzs)])
end
xlabel('horizon','interpreter','latex');
ylim(ylim_cover);
title('Coverage (frac DGP $>= 0.7$)','interpreter','latex');
grid on
set(gca, 'TickLabelInterpreter', 'latex')



% median length
subplot(1,3,3)
hold on;
for j=1:numproc
    ind_esttype  = procs(j, 1);
    ind_citype   = procs(j, 2);
    plot(horzs, squeeze(log10(median_length(procs(j,1),:,procs(j,2)))), ...
        line_specs{ind_citype}, 'Color', line_colors(ind_esttype,:),'LineWidth',LW);
end
hold off;
if strcmp(estimand_type, 'obsshock')
    xlim([min(horzs) max(horzs)])
else
    xlim([min(horzs)+1 max(horzs)])
end
xlabel('horizon','interpreter','latex');
ylim([min(yticks_length) max(yticks_length)]);
yticks(yticks_length);
yticklabels(yticklabels_length);
title('median length, log scale','interpreter','latex');

lgd_name = string(proc_names(procs(:, 1)));
lgd_name(procs(:, 2)> 1) = strcat(lgd_name(procs(:, 2)> 1), '$_b$');

n_columns = 1+(length(lgd_name)>4)*1;

legend(lgd_name, 'Location', 'southwest','interpreter','latex', ...
    'NumColumns', n_columns);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11)


% size
set(gcf, 'Units', 'inches', 'Position',[2 2 8 4])


end