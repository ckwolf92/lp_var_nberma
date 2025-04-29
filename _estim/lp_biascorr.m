function irs_corr = lp_biascorr(irs, w)
    % LP bias correction
    % "BCC" estimator in Herbst & Johanssen (2022)
    %
    % Inputs:
    % irs  1xH   impulse response. First element is response on impact.
    % w    Txk   Control variables  
    %
    % Outputs:
    % irs_corr 1xh  Bias-corrected impulse response.

    irs_hor = length(irs)-1;
    T       = size(w,1);

    % ACF term in bias correction
    acf_corr = nan(1,irs_hor);
    w = w - mean(w); % de-mean
    Sigma_0 = cov(w); % var-cov
    for j=1:irs_hor
        Sigma_j = (w(1:end-j,:)'*w(j+1:end,:))/(T-j-1); % j-th autocov
        acf_corr(j) = 1+trace(Sigma_0\Sigma_j);
    end

    % Iterate on bias correction
    irs_corr = irs;
    for h=1:irs_hor
        irs_corr(h+1) = irs(h+1) + (1/(T-h))*acf_corr(1:h)*irs_corr(h:-1:1)';
    end

end