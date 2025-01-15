if strcmp(dgp_type_plot,'both')

    load(fullfile(load_folder, strcat('dfm_', 'g', '_', estimand_type)))
    results_g = results;
    results_g.target_irf = DF_model.target_irf;
    results_g.M = DF_model.M;
    clear results DF_model DFM_estimate dgp_type

    load(fullfile(load_folder, strcat('dfm_', 'mp', '_', estimand_type)))
    results_mp = results;
    results_mp.target_irf = DF_model.target_irf;
    results_mp.M = DF_model.M;
    clear results DF_model DFM_estimate dgp_type

    results.target_irf = [results_g.target_irf,results_mp.target_irf];

    results.bias2 = [results_g.bias2;results_mp.bias2];
    results.vce   = [results_g.vce;results_mp.vce];
    results.mse   = 0.5 * results.bias2 + 0.5 * results.vce;

    results.coverage_prob  = [results_g.coverage_prob;results_mp.coverage_prob];
    results.coverage_avg   = mean(results.coverage_prob,1);
    results.coverage_indic = mean(results.coverage_prob >= covg_cutoff,1);

    results.median_length = [results_g.median_length;results_mp.median_length];  
    results.median_avg    = mean(results.median_length,1);

    results.M = [results_g.M;results_mp.M];

    clear results_g results_mp
    
else

    load(fullfile(load_folder, strcat('dfm_', dgp_type_plot, '_', estimand_type)))

    results.target_irf     = DF_model.target_irf;
    results.mse            = 0.5 * results.bias2 + 0.5 * results.vce;
    results.coverage_avg   = mean(results.coverage_prob,1);
    results.coverage_indic = mean(results.coverage_prob >= covg_cutoff,1);
    results.median_avg     = mean(results.median_length,1);
    results.M              = DF_model.M;
    clear DF_model DFM_estimate dgp_type

end