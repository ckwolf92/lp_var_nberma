function y = sim_arma11(T, rho, alpha, sigma2)
% gen_arma11  Generate ARMA(1,1) process with Gaussian shocks
%
% Model: y_t = rho y_t-1 + epsilon_t + alpha*epsilon_t-1
% Input:
%  - T:      Sample length
%  - rho:    AR parameter
%  - alpha:  MA parameter
%  - sigma2: Shock variance
%
% Output:
%  - y: simulated data

y      = nan(T+1, 1);
y(1)   = normrnd(0, 1)*sqrt((1+2*rho*alpha+alpha^2)*sigma2/(1-rho^2));
shocks = normrnd(0, 1, [T+1,1])*sqrt(sigma2);

for t = 2:length(y)
    y(t) = rho*y(t-1)+alpha*shocks(t-1)+shocks(t);
end

y = y(2:end);

end