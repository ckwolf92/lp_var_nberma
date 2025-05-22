%% LP-VAR EMPIRICAL COMPARISON
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 04/29/2025

%% HOUSEKEEPING

clear
clc
close all

path = cd;

addpath(genpath([path '../_auxiliary_functions']));
addpath(genpath([path '../_results']));

cd(path);

%% IMPORT RESULTS

load res_application

plot_folder = '_figures/lpvar_illustr';
mkdir(plot_folder)

%% SETTINGS

%----------------------------------------------------------------
% Variables and Horizons
%----------------------------------------------------------------

var_select = {[1 2 3 4]; ...
              [1 2 3]; ...
              [1 2 3]; ...
              [1 2 3]; ...
              [1 2 3 4]; ...
              [1 2 3]; ...
              [1 2 3]; ...
              [1 2 3]};

hor_select = {[1:1:12]; ...
              [1:1:4]; ...
              [1:1:4]; ...
              [1:1:4];...
              [13:1:49]; ...
              [5:1:21]; ...
              [5:1:21]; ...
              [5:1:21]};

n_appl = length(var_select);

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.lgrey  = [200/255 200/255 200/255];
settings.colors.grey   = [170/255 170/255 170/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.5 * settings.colors.blue + 0.5 * [1 1 1];

settings.colors.list = [settings.colors.grey; settings.colors.blue];

%% EXTRACT RESULTS

%----------------------------------------------------------------
% SE Ratios
%----------------------------------------------------------------

ses     = cell(n_appl,1);

for i_appl = 1:n_appl

    ses_tmp = NaN(length(hor_select{i_appl}),length(var_select{i_appl}));
    for i_hor = 1:length(hor_select{i_appl})
        for i_var = 1:length(var_select{i_appl})
            if i_appl <= 4
                ses_tmp(i_hor,i_var) = appl(i_appl).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    ./ appl(i_appl).results.ses_boot(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
            else
                ses_tmp(i_hor,i_var) = appl(i_appl-4).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    ./ appl(i_appl-4).results.ses_boot(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
            end
        end
    end
    ses_tmp = ses_tmp(:);
    ses{i_appl} = ses_tmp;
    clear ses_tmp

end

ses_all_1 = [ses{1};ses{2};ses{3};ses{4}]; % short horizons
ses_all_2 = [ses{5};ses{6};ses{7};ses{8}]; % long horizons

%----------------------------------------------------------------
% Point Estimate Gap
%----------------------------------------------------------------

estim_delta     = cell(n_appl,1);

for i_appl = 1:n_appl

    estim_delta_tmp = NaN(length(hor_select{i_appl}),length(var_select{i_appl}));
    for i_hor = 1:length(hor_select{i_appl})
        for i_var = 1:length(var_select{i_appl})
            if i_appl <= 4
%                 scale = 1/max(abs(squeeze(appl(i_appl).results.estims(1,var_select{i_appl}(i_var),:))));
                scale = 1/appl(i_appl).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
                estim_delta_tmp(i_hor,i_var) = scale * abs(appl(i_appl).results.estims(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    - appl(i_appl).results.estims(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)));
            else
%                 scale = 1/max(abs(squeeze(appl(i_appl-4).results.estims(1,var_select{i_appl}(i_var),:))));
                scale = 1/appl(i_appl-4).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
                estim_delta_tmp(i_hor,i_var) = scale * abs(appl(i_appl-4).results.estims(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    - appl(i_appl-4).results.estims(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)));
            end
        end
    end
    estim_delta_tmp = estim_delta_tmp(:);
    estim_delta{i_appl} = estim_delta_tmp;
    clear estim_delta_tmp scale

end

estim_delta_all_1 = [estim_delta{1};estim_delta{2};estim_delta{3};estim_delta{4}]; % short horizons
estim_delta_all_2 = [estim_delta{5};estim_delta{6};estim_delta{7};estim_delta{8}]; % long horizons

%% PLOT RESULTS

x = 1:2;
ses_box         = [[ses_all_1;NaN(length(ses_all_2)-length(ses_all_1),1)],ses_all_2];
estim_delta_box = [[estim_delta_all_1;NaN(length(estim_delta_all_2)-length(estim_delta_all_1),1)],estim_delta_all_2];

plotwidth = 0.4;
gapsize = 0.06;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

figure(1)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',17)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:2
    boxchart(x(i)*ones(size(ses_box(:,i))), ses_box(:,i), ...
        'MarkerStyle', 'none', 'BoxFaceColor', settings.colors.list(i,:), 'LineWidth', 2);
    hold on
end
set(gcf,'color','w')
title('Ratio of SEs (VARs/LPs)','interpreter','latex','FontSize',18)
yticks([0:0.2:1])
xticklabels({})
grid on
hold off

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',17)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:2
    boxchart(x(i)*ones(size(estim_delta_box(:,i))), estim_delta_box(:,i), ...
        'MarkerStyle', 'none', 'BoxFaceColor', settings.colors.list(i,:), 'LineWidth', 2);
    hold on
end
set(gcf,'color','w')
title('Scaled LP-VAR $\Delta$','interpreter','latex','FontSize',18)
ylim([0 12])
yticks([0:2:12])
xticklabels({})
legend({'Short Horizons','Long Horizons'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.8*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'lp_vs_var_ramey.eps'))