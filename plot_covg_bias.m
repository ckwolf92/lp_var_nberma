%% COVERAGE COST OF BIAS
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 01/16/2025

%% HOUSEKEEPING

clear
clc
close all

path = cd;

cd([path]);

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

settings.colors.black  = [0 0 0];
settings.colors.lgrey  = [200/255 200/255 200/255];
settings.colors.grey   = [170/255 170/255 170/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.5 * settings.colors.blue + 0.5 * [1 1 1];
settings.colors.dblue  = 0.7 * settings.colors.blue + 0.3 * [0 0 0];

plot_folder = '_figures/uncertainty';
mkdir(plot_folder)

figure(1)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',17)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(b_over_sigma,coverage_68,'linewidth',4,'linestyle',':','color',settings.colors.lblue);
hold on
plot(b_over_sigma,coverage_90,'linewidth',4,'linestyle','-','color',settings.colors.blue);
hold on
plot(b_over_sigma,coverage_95,'linewidth',4,'linestyle','-.','color',settings.colors.dblue);
hold off
set(gcf,'color','w')
xlabel('$|b_h(p)| / \tau_{h,VAR}(p)$','interpreter','latex','FontSize',20)
ylabel('VAR CI Coverage','interpreter','latex','FontSize',20)
yticks([0:0.2:1])
legend({'68\% target','90\% target','95\% target'},'Location','Northeast','fontsize',20,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(plot_folder, 'covg_distortion.eps'))