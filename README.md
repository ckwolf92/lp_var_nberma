# lp_var_nberma

MATLAB code for large-scale simulation studies on impulse responses using Local Projections, VARs, and several of their variants.

**Reference**: [Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2025), ''Local Projections or VARs? A Primer for Macroeconomists''

Tested in: MATLAB R024a on MacBook Pro 2023 (M3 Pro)

## Contents

**[documents/](documents/)**: Paper and supplement XX to add XX
- Main paper
- Supplement

**[_estim/](_estim/)**: Impulse response estimation functions
- [ir_estim.m](_estim/ir_estim.m): Wrapper function for (Bayesian) VAR and (smoothed) LP impulse responses
- [GLP/](_estim/GLP/): Bayesian VAR ([Giannone, Lenza, and Primiceri, 2015](https://direct.mit.edu/rest/article-abstract/97/2/436/58236/Prior-Selection-for-Vector-Autoregressions?redirectedFrom=fulltext)) 
- [LP_Penalize/](_estim/LP_Penalize/): Penalized LP ([Barnichon and Brownlees, 2019](https://direct.mit.edu/rest/article-abstract/101/3/522/58522/Impulse-Response-Estimation-by-Smooth-Local?redirectedFrom=fulltext))
- [lp_biascorr.m](_estim/lp_biascorr.m): Local projection bias correction as in [Herbst and Johannsen (2024)](https://www.sciencedirect.com/science/article/abs/pii/S0304407624000010)
- [var_biascorr.m](_estim/var_biascorr.m): VAR bias correction as in [Pope (1990)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9892.1990.tb00056.x)

**[_dfm/](_dfm/)**: Functions and scripts used for the encompassing Dynamic Factor Model (DFM) used in the simulation study
- [settings/](_dfm/settings/): Folder with simulation settings
- [subroutines](_dfm/subroutines/): Folder with supporting subroutines used for the simulations
- [sw_estim](_dfm/sw_estim/): DFM code and data adapted from [Lazarus, Lewis, Stock and Watson (2018)](https://www.tandfonline.com/doi/full/10.1080/07350015.2018.1506926)

**[illustrations/](illustrations/)**: Illustrative figures featured in the main text
- [plot_dfm_var_vs_lp.m](illustrations/plot_dfm_var_vs_lp.m): Plot of LP and VAR impulse response estimands by lag length (Figure 2.1 of the main text)
- [plot_emp_var_vs_lp.m](illustrations/plot_emp_var_vs_lp.m): Box plots of VAR-to-LP standard error ratios and normalized point estimate differences based on applications of [Ramey (2016)](https://www.sciencedirect.com/science/article/abs/pii/S1574004816000045) (Figure 3.1 of the main text)
- [res_application.mat](illustrations/res_application.mat): Data used in [plot_emp_var_vs_lp.m](illustrations/plot_emp_var_vs_lp.m). These are taken from [ Montiel Olea, Plagborg-Møller, Qian, and Wolf (2024)](https://github.com/ckwolf92/lp_var_inference/tree/main/emp_ses).
- [plot_arma.m](illustrations/plot_arma.m): Plots of impulse response estimates, absolute bias, and standard deviation (Figures 3.2 and 3.3 of the main text)
- [plot_covg_bias.m](illustrations/plot_covg_bias.m): Asymptotic coverage of VAR confidence interval as a function of relative VAR bias (Figure 6.1 of the main text)


**[simulations/](simulations)**: Running and generating figures for the simulation study 
- [main_dfm.m](simulations/main_dfm.m): Main file for the simulation study
- [plot_dfm_paper.m](simulations/plot_dfm_paper.m): Plots of bias, standard deviation, MSE, and confidence interval coverage across estimators and DGPs (Figures 5.1, 5.2, and 6.2 in the main text, Figures D.1-D.7 in the supplement).
- [plot_varinlpbands.m](simulations/plot_varinlpbands.m): Plot of fraction of DGPs for which the VAR point estimates are contained inside the LP confidence interval (Figure E.1 of the supplement).

## Estimating IRFs using local projections

Using our code suite, our recommended procedure is implemented below.

Use the AIC to select the lag length `p` by estimating auxiliary VARs of lag length 1 through `p_max`.
```m
p_max  = 20;   % Maximum lag length
method = 1;    % AIC
p      = ic_var(data_y, p_max, method);
```

Below returns the local projection impulse response function `irs` of the second variable to the first innovation.
```m
horzs   = 0:20;  % Horizons of interest
opts_lp = {'resp_ind',  2, 'innov_ind', 1, 'estimator', 'lp',...
    'shrinkage', false,...
    'bias_corr', true,...  % Herbst and Johannsen (2024) bias correction
    'alpha',     0.05,...  % Significance level
    'se_homosk', false};   % Eicher-Huber-White standard errors};

[irs, ses, cis] = ir_estim(data_y, p, horzs, opts_lp{:});  
```
We recommend using the procedure of [Herbst and Johannsen (2024)](https://www.sciencedirect.com/science/article/abs/pii/S0304407624000010) to bias-correct the LP estimator with heteroskedasticity-robust standard errors. Delta method standard errors and confidence intervals are stored in ``ses`` and ``cis`` respectively.

We find further finite-sample performance gains with the residual block bootstrap of [Brüggemann,
Jentsch, and Trenkler (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0304407615002547) using the block length determined by [Jentsch and Lunsford (2019)](https://www.aeaweb.org/articles?id=10.1257/aer.20162011).
```m
block_length = ceil(5.03*size(data_y, 1)^(1/4));  % Jentsch and Lunsford (2019) rule of thumb
opts_lp_boot = {'bootstrap', 'var',...  % Residual bootstrap
                'boot_num', 1000,...    % Number of bootstrap draws
                'boot_blocklength', block_length};
[~, ~, ~, cis_boot] = ir_estim(data_y, p, horzs, opts_lp{:}, opts_lp_boot{:}); 
```
We recommend reporting percentile-t bootstrap confidence intervals `cis_boot(:,:,3)`.



## Detailed replication instructions

### Illustrative figures
- *Figure 2.1*. Run the file [plot_dfm_var_vs_lp.m](illustrations/plot_dfm_var_vs_lp.m).
- *Figure 3.1*. Run the file [plot_emp_var_vs_lp.m](illustrations/plot_emp_var_vs_lp.m).
- *Figures 3.2 and 3.3*. Run the file [plot_arma.m](illustrations/plot_arma.m).
- *Figure 6.1*: Run the file [plot_covg_bias.m](illustrations/plot_covg_bias.m). 

### Figures from the simulation studies

1. **Estimate IRFs from simulated data**: This produces raw IRF estimates and confidence intervals for each estimator.
    - In [_dfm/settings/shared.m](_dfm/settings/shared.m), set `settings.simul.n_mc` to `1000` and `settings.specifications.random_n_spec` to `100`.
    - Run the file `simulations/main_dfm.m` with the following settings:
        - Run combinations of `dgp_type` (set to `g` or `mp`) and `mode_type` (set to `3` or `5`).
        - To produce the appendix results for recursive identification, run the same combinations for `estimand_type='recursive'`.
    - Results will be saved in the directory `simulations/_results/`.
     
2. **Generate figures**: The following steps generate the figures included in the main text and supplement.
    - Run the file [plot_dfm_paper.m](simulations/plot_dfm_paper.m) with the following figure-specific header modifications:
        - *Figures 5.1, 5.2, 6.2, and D.7*: Set `dgp_type_plot='both'`, `estimand_type='obsshock'`, and `mode_type=6`.
        - *Figures D.1 and D.2*: Set `dgp_type_plot='both'`, `estimand_type='recursive'`, and `mode_type=6`.
        - *Figures D.3 and D.4*: Set `dgp_type_plot='g'`, `estimand_type='obsshock'`, and `mode_type=6`.
        - *Figures D.5 and D.6*: Set `dgp_type_plot='mp'`, `estimand_type='obsshock'`, and `mode_type=6`.
    - *Figure E.1*: Run the file [plot_varinlpbands.m](simulations/plot_varinlpbands.m).
## Acknowledgements
The simulation study extends that of [Li, Plagborg-Møller, and Wolf (2024)](https://github.com/dake-li/lp_var_simul) with estimation routines based on the replication materials for [Montiel Olea, Plagborg-Møller, Qian, and Wolf (2024)](https://github.com/ckwolf92/lp_var_inference).

We rely on the BVAR code of [Domenico Giannone, Michele Lenza, and Giorgio Primiceri](http://faculty.wcas.northwestern.edu/gep575/GLPreplicationWeb.zip) and the penalized LP code of [Regis Barnichon and Christian Brownlees](https://drive.google.com/drive/folders/1Fjzw-U3hjIl467KXywRqeQod2jdHOmDo). We also use the Dynamic Factor Model code and data of [Eben Lazarus, Daniel Lewis, Jim Stock, and Mark Watson](http://www.princeton.edu/~mwatson/ddisk/LLSW_ReplicationFiles_071418.zip), modified to allow for conditional heteroskedasticity of the shocks for both the factors and idiosyncratic disturbances.

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049) and by the Alfred P. Sloan Foundation, and Wolf does the same for NSF [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).

