function irf = irf_arma11(horzs, rho, alpha)
% irf_arma11  IRF of ARMA(1,1) process
% Model: y_t = rho y_t-1 + eps_t + alpha eps_t-1
% Input:
%   - horzs: Horizon of interest
%   - rho: AR parameter
%   - alpha: MA parameter
% Output:
%   - irf: Impulse response function

irf           = rho.^horzs + alpha * rho .^(horzs-1);
irf(horzs==0) = 1;

end