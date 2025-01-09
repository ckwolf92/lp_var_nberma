%% SPECIFIC SETTINGS FOR ALTERNATIVE EXPERIMENTS

% set up directory for robustness-check modes

mode_list   = {'baseline', 'persistent', 'salient', 'more'};
save_mode_dir = mode_list{mode_type};

% rewrite some baseline settings in "shared.m" for different robustness check modes

switch mode_type

    case 1 % baseline
        
        % rewrite nothing and use all the settings in "shared.m"

    case 2 % persistent DFM

        DF_model.levels = 1;
        DF_model.n_lags_fac = 4;
        DF_model.n_lags_uar = 4;

    case 3 % salient series

        % selection of DGPs from encompassing model

        settings.specifications.random_from_key_series = 1; % randomly select from some key series in DFM list?

    case 4 % more observables
        
        % selection of DGPs from encompassing model

        settings.specifications.random_n_var = 7; % larger number of variables in each specification

    case 5  % Salient persistent
        settings.specifications.random_from_key_series = 1; % randomly select from some key series in DFM list?
        DF_model.levels                                = 1;
        DF_model.n_lags_fac                            = 4;
        DF_model.n_lags_uar                            = 4;

end