%% VAR ESTIMATE IN LP CI: GENERATE FIGURE
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 05/12/2025

%% HOUSEKEEPING
clear
clc
close all

estimand_type = 'obsshock';  % 'obsshock' or 'recursive'
dgp_type      = {'g', 'mp'};
mode_type     = {'salient', 'persistent_salient'};
sample_length = 'medium';  % 'short' (T=100), 'medium' (T=240), 'long' (T=720)
shock_type    = 'arch';    % 'iid' or 'arch'

plot_pre    = '_figures';

%% LOAD AND COMBINE RESULTS

ind = 1;
for j_mode = 1:length(mode_type)
    for j_dgp = 1:length(dgp_type)
        
        % Load results
        Res = load([fullfile('_results', mode_type{j_mode}, ...
            ['dfm_' dgp_type{j_dgp} '_' estimand_type '_' sample_length...
            '_' shock_type '.mat'])]);
        VARinLP_probs_j = Res.results.VARinLP_prob;

        % Store
        if ind == 1
            VARinLP_inds = nan(size(VARinLP_probs_j, 1), size(VARinLP_probs_j, 2), ...
                               length(mode_type)*length(dgp_type));
        end
        
        VARinLP_inds(:,:,ind) = VARinLP_probs_j;
        ind                     = ind+1;
    end
end


%% FIGURE: PROPORTION VAR IN LP CI BY HORIZON

disp(1-mean(VARinLP_inds(:)))

horzs = Res.settings.est.IRF_select - 1;

figure('Units', 'inches', 'Position', [2 2 6 3])
plot(horzs, mean(mean(mean(VARinLP_inds, 1), 3),4), 'LineWidth', 5)
ylim([0, 1])
grid on;
xlabel('horizon', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)
xticks(0:2:20)
exportgraphics(gcf, fullfile(plot_pre, 'var_in_lpci.eps'))
