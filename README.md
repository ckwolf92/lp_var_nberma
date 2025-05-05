# lp_var_nberma

MATLAB code for large-scale simulation studies on impulse responses using Local Projections, VARs, and several of their variants.

**Reference**: [Montiel Olea, José Luis](https://www.joseluismontielolea.com), [Mikkel Plagborg-Møller](https://www.mikkelpm.com), [Eric Qian](https://www.eric-qian.com) and [Christian K. Wolf](https://www.christiankwolf.com/) (2025), ''Local Projections or VARs? A Primer for Macroeconomists''

Tested in: MATLAB R024a on MacBook Pro 2023 (M3 Pro)

## Contents

## Estimating IRFs using local projections

Using our code suite, our recommended procedure is implemented below.

Use the AIC to select the lag length `p` through estimating auxiliary VARs of lag length 1 through `p_max`.
```m
p_max  = 20;   % Maximum lag length
method = 1;    % AIC
p      = ic_var(data_y, p_max, method);
```

Below returns the local projection impulse response function `irs` of the second variable to the first innovation.
```m
horzs        = 0:20;  % Horizons of interest
block_length = ceil(5.03*size(data_y, 1)^(1/4));  % Jentsch and Lunsford (2019) rule of thumb

[irs, ses, cis, cis_boot] = ir_estim(data_y, p, horzs, ...
    'resp_ind',  2,... 
    'innov_ind', 1,...
    'estimator', 'lp',...
    'shrinkage', false,...
    'bias_corr', true,...   % Herbst and Johannsen (2024) bias correction
    'alpha',      0.05,...  % Significance level
    'se_homosk', false,...  % Eicher-Huber-White standard errors
    'bootstrap', 'var',...  % Residual bootstrap
    'boot_num', 1000,...    % Number of bootstrap draws
    'boot_blocklength', block_length);  % Bootstrap block length 
```
We recommend using the procedure of [Herbst and Johannsen (2024)](https://www.sciencedirect.com/science/article/abs/pii/S0304407624000010) to bias-correct the LP estimator with heteroskedasticity-robust standard errors. Delta method standard errors and confidence intervals are stored in ``ses`` and ``cis`` respectively.

We find further gains in finite-sample performance with the residual block bootstrap of [Brüggemann,
Jentsch, and Trenkler (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0304407615002547) with percentile-t bootstrap confidence intervals (stored in `cis_boot(:,:,3)`). We implement using the block length determined by [Jentsch and Lunsford (2019)](https://www.aeaweb.org/articles?id=10.1257/aer.20162011).



## Detailed replication instructions

XX Set `settings.simul.n_mc` to `1000` and `settings.specifications.random_n_spec` to `100`.

## Acknowledgements
The simulation study extends that of [Li, Plagborg-Møller, and Wolf (2024)](https://github.com/dake-li/lp_var_simul) with estimation routines based on the replication materials for [Montiel Olea, Plagborg-Møller, Qian, and Wolf (2024)](https://github.com/ckwolf92/lp_var_inference).

We rely on the BVAR code of [Domenico Giannone, Michele Lenza, and Giorgio Primiceri](http://faculty.wcas.northwestern.edu/gep575/GLPreplicationWeb.zip) and the penalized LP code of [Regis Barnichon and Christian Brownlees](https://drive.google.com/drive/folders/1Fjzw-U3hjIl467KXywRqeQod2jdHOmDo). We also use the Dynamic Factor Model code and data of [Eben Lazarus, Daniel Lewis, Jim Stock, and Mark Watson](http://www.princeton.edu/~mwatson/ddisk/LLSW_ReplicationFiles_071418.zip), modified to allow for conditional heteroskedasticity of the shocks for both the factors and idiosyncratic disturbances.

Plagborg-Møller acknowledges that this material is based upon work supported by the NSF under [Grant #2238049](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2238049) and by the Alfred P. Sloan Foundation, and Wolf does the same for NSF [Grant #2314736](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2314736).

