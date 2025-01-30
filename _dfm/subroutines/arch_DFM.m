%% ORGANIZE
function DFM_estimate = arch_DFM(DFM_estimate)
%% SETUP

% Store residuals
uar_res = DFM_estimate.uar_res';
fac_res = DFM_estimate.fac_res;

% Univariate AR residual ARCH parameters
arch_uar     = nan(size(uar_res, 2), 1);
archcons_uar = arch_uar;
arch_uar_se  = arch_uar;

% Factor residual ARCH parameters
arch_fac     = nan(size(fac_res, 2), 1);
archcons_fac = arch_fac;
arch_fac_se  = arch_fac;

mdl_arch1        = garch(0, 1);  
mdl_arch1.Offset = NaN;  % Include an offset term (uar_res isn't necessarily mean zero)

%% ESTIMATE

for i = 1:length(arch_uar)  % Idiosyncratic
    est_mdl         = estimate(mdl_arch1, uar_res(:, i), "Display", "off");
    archcons_uar(i) = est_mdl.Constant;
    arch_uar(i)     = est_mdl.ARCH{1};
    arch_uar_se(i)  = est_mdl.summarize.Table{2,2};
end


for i = 1:length(arch_fac)  % Factors
    est_mdl         = estimate(mdl_arch1, fac_res(:, i), "Display", "off");
    archcons_fac(i) = est_mdl.Constant;
    arch_fac(i)     = est_mdl.ARCH{1};
    arch_fac_se(i)  = est_mdl.summarize.Table{2,2};
end

DFM_estimate.arch_uar = arch_uar;
DFM_estimate.arch_fac = arch_fac ;

end