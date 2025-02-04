% shock settings
if strcmp(shock_type, 'arch')
    DFM_estimate      = arch_DFM(DFM_estimate);
    DF_model.arch_uar = DFM_estimate.arch_uar;
    DF_model.arch_uar(DF_model.arch_uar>DF_model.censor_arch_uar) = DF_model.censor_arch_uar;
    DF_model.arch_fac = DFM_estimate.arch_fac*0 + mean(DFM_estimate.arch_fac);

else  % iid case (arch parameters are unused)
    DF_model.arch_uar = [];
    DF_model.arch_fac = [];    

end

DF_model.shock_type    = shock_type;
