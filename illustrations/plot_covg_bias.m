%% COVERAGE COST OF BIAS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 04/29/2025

%% HOUSEKEEPING

clear
clc
close all

path = cd;

cd(path);

%% COMPUTE COVERAGE DISTORTION

b_over_sigma = 0:0.01:3; 
covg_targets = [0.68 0.9 0.95];
norm_cutoff    = NaN(length(covg_targets),1);
for i_covg = 1:length(covg_targets)
    norm_cutoff(i_covg) = norminv(1 - (1 - covg_targets(i_covg))/2);
end

coverage_68  = normcdf(norm_cutoff(1)-b_over_sigma,0,1) - normcdf(-norm_cutoff(1)-b_over_sigma,0,1);
coverage_90  = normcdf(norm_cutoff(2)-b_over_sigma,0,1) - normcdf(-norm_cutoff(2)-b_over_sigma,0,1);
coverage_95  = normcdf(norm_cutoff(3)-b_over_sigma,0,1) - normcdf(-norm_cutoff(3)-b_over_sigma,0,1);

%% PLOT RESULTS

colors.blue  = [102/255 178/255 255/255];
colors.dblue = 0.7 * colors.blue;
colors.lblue = 0.4 * colors.blue + 0.6 * [1 1 1];
colors.red   = [204/255 0/255 0/255];
colors.dred  = 0.7 * colors.red;
colors.lred  = 0.4 * colors.red + 0.6 * [1 1 1];

plot_folder = '_figures/uncertainty';
mkdir(plot_folder)

figure(1)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',17)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(b_over_sigma,coverage_90,'linewidth',4,'linestyle','-','color',colors.red);
hold on
plot(0,0.9,'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',10)
hold on
plot(b_over_sigma,coverage_68,'linewidth',4,'linestyle',':','color',colors.lred);
hold on
plot(0,0.68,'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',10)
hold on
plot(b_over_sigma,coverage_95,'linewidth',4,'linestyle','-.','color',colors.dred);
plot(0,0.95,'s','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',10)
hold on
hold off
set(gcf,'color','w')
xlabel('$|b_h(p)| / \tau_{h,VAR}(p)$','interpreter','latex','FontSize',20)
ylabel('CI Coverage','interpreter','latex','FontSize',20)
yticks([0:0.2:1])
% legend({'LPs','VARs'},'Location','Southwest','fontsize',20,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'covg_distortion.eps'))