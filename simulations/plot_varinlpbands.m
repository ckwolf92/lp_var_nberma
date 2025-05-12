%% PRELIMINARIES
clear
clc
close all

estimand_type = 'obsshock';
dgp_type      = {'g', 'mp'};
mode_type     = {'salient', 'persistent_salient'};

%% LOAD AND COMBINE RESULTS
% Add both objects to the code in the dev branch.
% Keep the VARinLP_inds in the big file.
% Add VARinLP that averages over the 1000 simulations.
% Keep only VARinLP in the small file.


ind = 1;
for j_mode = 1:length(mode_type)
    for j_dgp = 1:length(dgp_type)
        
        Res = load([fullfile('_results_combined', mode_type{j_mode}, ...
            ['dfm_' dgp_type{j_dgp} '_' estimand_type '_medium_arch.mat'])]);
        VARinLP_inds_ij = squeeze(Res.results.estims(:,1,:,:)) >= squeeze(Res.results.cis_lower(:,5,:,4,:)) &...
                          squeeze(Res.results.estims(:,1,:,:)) <= squeeze(Res.results.cis_upper(:,5,:,4,:));

        if ind == 1
            VARinLP_inds = nan(size(VARinLP_inds_ij, 1), size(VARinLP_inds_ij, 2), ...
                               size(VARinLP_inds_ij, 3),length(mode_type)*length(dgp_type));
        end
        
        VARinLP_inds(:,:,:,ind) = VARinLP_inds_ij;
        ind                     = ind+1;
    end
end


%% FIGURE: PROPORTION VAR IN LP CI BY HORIZON

disp(1-mean(VARinLP_inds(:)))

horzs = Res.settings.est.IRF_select - 1;

figure('Units','inches','Position', [2 2 6 3])
plot(horzs, mean(mean(mean(VARinLP_inds, 1), 3),4), 'LineWidth', 5)
ylim([0, 1])
grid on;
xlabel('horizon', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex', 'FontSize', 12)
xticks(0:2:20)
mean(1-VARinLP_inds(:))
exportgraphics(gcf, 'var_in_lpci.eps')
