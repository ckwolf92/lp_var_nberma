%% ILLUSTRATIVE ARMA
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 04/29/2025
clear
clc
close all

path        = cd;
plot_folder = '_figures/arma';
mkdir(plot_folder)

addpath(genpath('../_auxiliary_functions'));
addpath(genpath('../_estim'))


%% SETTINGS

%----------------------------------------------------------------
% DGP settings
%----------------------------------------------------------------

dgp        = struct;
dgp.T      = 240;
dgp.rhos   = [.85, .95];
dgp.alphas = [.1, .2];
dgp.sigma2 = 1;

%----------------------------------------------------------------
% Simulation settings
%----------------------------------------------------------------

sim           = struct;
sim.n_rep     = 1e4;
sim.rng_seed  = 2025;
sim.n_workers = 11;

%----------------------------------------------------------------
% Estimator settings
%----------------------------------------------------------------

settings.est.horzs = 0:12;
settings.est.methods_shared = {'alpha', 0.1,...
    'no_const',  false,...
    'se_homosk', true,...
    'boot_num',  0,...
    'bootstrap', [],...
    'resp_ind',  1,...
    'innov_ind', 1,...
    'bias_corr', false,...
    'shrinkage', false};
[settings.est.methods{1:3}] = deal({'estimator', 'var'});
[settings.est.methods{4:6}] = deal({'estimator', 'lp'});
settings.est.n_lags_fix     = [1; 4; 8; 1; 4; 8];

%% SETUP

rng(sim.rng_seed, 'twister');      

% Combinations of DGP parameters
aux1   = repmat(dgp.rhos,size(dgp.alphas));
aux2   = repmat(dgp.alphas',size(dgp.rhos))';
dgps   = [aux1;reshape(aux2,[1,size(aux1,2)])];

settings.est.n_methods = length(settings.est.methods);

n_dgp     = size(dgps,2);                % No. of DGPs
n_horz    = length(settings.est.horzs);  % No. of horizons
n_methods = settings.est.n_methods;      % No. of regression specifications
n_rep     = sim.n_rep;                   % No. of repetitions

irs_true = nan(n_dgp, n_horz);

% Initialize results arrays
estims    = zeros(n_dgp, n_methods, n_horz, n_rep);
ses       = estims;
cis_lower = zeros(n_dgp, n_methods, n_horz, n_rep);
cis_upper = cis_lower;

timer = tic; % Start timer

for i_dgp = 1:n_dgp
    rho_i        = dgps(1,i_dgp);
    alpha_i      = dgps(2, i_dgp);
    i_rand_seeds = randi(2^32-1,1,n_rep);
    
    % ARMA
    irs_true(i_dgp, :) = irf_arma11(settings.est.horzs, rho_i, alpha_i);

    parfor(i=1:n_rep, sim.n_workers)
        rng(i_rand_seeds(i), 'twister')
        y = sim_arma11(dgp.T, rho_i, alpha_i, dgp.sigma2);

        
        i_estims    = zeros(n_methods, n_horz);        
        i_ses       = zeros(n_methods, n_horz);        
        i_cis_lower = zeros(n_methods, n_horz);         
        i_cis_upper = zeros(n_methods, n_horz);

        for j = 1:n_methods
            p = settings.est.n_lags_fix(j);
            
            [i_estims(j,:), i_ses(j,:), i_cis_dm] = ...
                ir_estim(y, p, settings.est.horzs,...
             settings.est.methods_shared{:}, settings.est.methods{j}{:});
            i_cis_lower(j,:) = i_cis_dm(1,:);
            i_cis_upper(j,:) = i_cis_dm(2,:);
        end

        estims(i_dgp,:,:,i)    = i_estims;
        ses(i_dgp,:,:,i)       = i_ses;
        cis_lower(i_dgp,:,:,i) = i_cis_lower;
        cis_upper(i_dgp,:,:,i) = i_cis_upper;

    end
end
sim.elapsed_time = toc(timer);

%% ORGANIZE RESULTS

dgp.irs_true     = irs_true;
irs_true_reshape = permute(repmat(irs_true, [1 1 6]), [1 3 2]);


results           = struct;
results.estims    = estims;
results.ses       = ses;
results.cis_lower = cis_lower;
results.cis_upper = cis_upper;

results.cover_inds ...
    = (irs_true_reshape >= results.cis_lower ...
      & irs_true_reshape <= results.cis_upper) + 0; 

results.coverage_prob = mean(results.cover_inds, 4); 

results.bias2 = (mean(results.estims,4) - irs_true_reshape).^2;
results.vce   = var(results.estims,1,4);


%% FIGURE: VAR AND LP DENSITY

close all

% Set colors
colors.blue    = [102/255 178/255 255/255];
colors.lblue   = 0.4 * colors.blue + 0.6 * [1 1 1];
colors.red     = [204/255 0/255 0/255];
colors.lred    = 0.4 * colors.red + 0.6 * [1 1 1];

BW          = .02;  % Bin width
h           = 2;    % Horizon
i_dgp       = 1;    % DGP
estims_dist = squeeze(results.estims(i_dgp, [1,4], h+1, :));

figure('Visible', 'on', 'Units', 'inches', 'Position', [2 2 6 3])

% Histogram version
p_lp = histogram(estims_dist(2,:), 'EdgeColor', colors.blue, 'FaceAlpha', .3,...
    'BinWidth', BW, 'FaceColor',colors.blue, 'EdgeAlpha', .3,...
    'Normalization','probability');
hold on
p_var = histogram(estims_dist(1,:), 'EdgeColor', colors.red, 'FaceAlpha', .5, ...
    'BinWidth', BW, 'FaceColor', colors.red, 'EdgeAlpha', .5,...
    'Normalization','probability');

ylabel('Probability', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex', 'FontSize', 12)

plot(dgp.irs_true(i_dgp,h+1)*ones(2,1), ylim, 'Color', 'k', ...
    'LineStyle', '--', 'LineWidth', 1)
legend([p_lp, p_var], {'LP(1)', 'VAR(1)'}, 'Interpreter','latex', ...
    'Location', 'northeast')
grid on;
exportgraphics(gcf, fullfile(plot_folder, 'arma_lpvar_dens.eps'))

%% BIAS STD DEV PLOTS

%----------------------------------------------------------------
% Plot by lag
%----------------------------------------------------------------

bias2_var = squeeze(results.bias2(1, 1:3, :));
bias2_lp  = squeeze(results.bias2(1, 4:6, :));
vce_var   = squeeze(results.vce(1, 1:3, :));
vce_lp    = squeeze(results.vce(1, 4:6, :));

bias2        = [bias2_var; bias2_lp];
vce         = [vce_var; vce_lp];
line_styles = {'-', '--', ':', '-', '--', ':'};
Cols        = [colors.red; colors.red; colors.red; colors.blue;colors.blue; colors.blue];


lgd_lab = strcat({'VAR(1)', 'VAR(4)', 'VAR(8)', 'LP(1)', 'LP(4)', 'LP(8)'});
reorder = [1 4 2 5 3 6];  % For the legend

plot_biasstd(settings.est.horzs, bias2(reorder,:), vce(reorder,:), ...
    lgd_lab(reorder), line_styles(reorder), Cols(reorder,:),  [0, .23])
exportgraphics(gcf, fullfile(plot_folder, 'arma_biasvce_lags.eps'))

%----------------------------------------------------------------
% Alternative DGPs
%----------------------------------------------------------------

close all
bias2_var = squeeze(results.bias2(1:3, 1, :));
bias2_lp  = squeeze(results.bias2(1:3, 4, :));
vce_var   = squeeze(results.vce(1:3, 1, :));
vce_lp    = squeeze(results.vce(1:3, 4, :));

bias2        = [bias2_var; bias2_lp];
vce         = [vce_var; vce_lp];
line_styles = {'-', '--', ':','-', '--', ':'};
Cols        = [colors.red;  colors.red;  colors.red; 
               colors.blue; colors.blue; colors.blue];

lgd_lab = {'VAR, baseline', ...
          'VAR, higher $\rho$', ...
          'VAR, higher $\alpha$', ...
            'LP, baseline', ...
            'LP, higher $\rho$', ...
            'LP, higher $\alpha$'};


reorder = [1 4 2 5 3 6];  % For the legend

plot_biasstd(settings.est.horzs, bias2(reorder,:), vce(reorder,:), ...
    lgd_lab(reorder), line_styles(reorder), Cols(reorder,:), [0, .23])
exportgraphics(gcf, fullfile(plot_folder, 'arma_biasvce_combined.eps'))


%% FUNCTIONS

function plot_biasstd(horzs, bias2, vce, lgd_lab, line_styles, Cols, yl)

LW     = 2;  % Line width
n_spec = size(bias2, 1);
fs     = 14;  % Font size
figure('Visible','on', 'Units','inches', 'Position', [2 2 6 4])
p = [];
for i = 1:n_spec
    subplot(5, 2, [1 3 5 7])
    p_i = plot(horzs, sqrt(bias2(i,:)), ...
        'Color', Cols(i, :), 'LineStyle', line_styles{i}, 'LineWidth', LW);
    hold on;
    p = [p, p_i];
    set(gca,'TickLabelInterpreter','latex', 'FontSize', fs)
    grid on;
    ylim(yl)
    xlim(horzs([1 end]))
    xticks(horzs(1:4:end))
    title('Abs. Bias', 'Interpreter', 'latex')

    subplot(5, 2, [2 4 6 8])
    plot(horzs, sqrt(vce(i,:)), ...
        'Color', Cols(i, :), 'LineStyle', line_styles{i}, 'LineWidth', LW);
    hold on;
    set(gca,'TickLabelInterpreter','latex', 'FontSize', fs )

    grid on;
    ylim(yl)
    xlim(horzs([1 end]))
    xticks(horzs(1:4:end))
    title('Std. Dev.', 'Interpreter','latex')
end

lgd = legend(lgd_lab, 'Interpreter', 'latex', 'Location', 'northeast', ...
    'NumColumns', 3, 'Position',  [0.1272    0.0677    0.7482    0.1056]);
set(gca, 'FontSize', fs);
end