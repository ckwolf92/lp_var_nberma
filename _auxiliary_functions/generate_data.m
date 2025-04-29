function data_sim = generate_data(model,settings)
% Function for generating simulated data
    % Use a general ABCD representation of the encompassing model (DFM, DSGE or others):
        % state transition:  s_t = A * s_{t-1} + B * e_t
        % measurement eq:    y_t = C * s_{t-1} + D * e_t

    % Warning: shock_weight' * e_t corresponds to true shock \epsilon_{1t} in our paper

% unpack settings

T      = settings.simul.T;
T_burn = settings.simul.T_burn;

A = model.ABCD.A;
B = model.ABCD.B;
C = model.ABCD.C;
D = model.ABCD.D;

[n_s,n_e] = size(B);

shock_weight = settings.est.shock_weight;

% draw shocks

if strcmp(model.shock_type, 'iid')
    data_e = randn(T_burn+T,n_e);
elseif strcmp(model.shock_type, 'arch')
    arch_all = [model.arch_fac;model.arch_uar];
    data_e   = NaN(T_burn+T,n_e);
    for i_e = 1:n_e
        arch_mdl = garch(Constant=1-arch_all(i_e),GARCH=0,ARCH=arch_all(i_e));
        [~,data_e(:,i_e)] = simulate(arch_mdl,T_burn+T);
    end
end

% simulate states & measurement error

s = zeros(n_s,1);
data_s = NaN(T_burn + T, n_s);
for t = 1:(T_burn + T)
    s = A * s + B * data_e(t,:)';
    data_s(t,:) = s';
end

% simulate observables

data_y = data_s(T_burn:end-1,:)*C' + data_e(T_burn+1:end,:)*D';

% collect results and shift timing

data_sim.data_y     = data_y;
data_sim.data_shock = data_e((T_burn+1):(T_burn+T), :) * shock_weight;

end