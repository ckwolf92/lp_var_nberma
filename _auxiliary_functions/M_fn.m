function M = M_fn(VAR_p,VMA,settings);

    % Get DGP for Simulations
    
    % Inputs:
    % VAR_p     struct          population VAR(p)
    % VMA       struct          residual VMA process
    % model     struct          ABCD representation of DGP
    % settings  struct          estimation settings
    
    % Output:
    % M         1x1             degree of mis-specification


    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------

    n_y     = size(VAR_p.Sigma_u,1);
    VMA_hor = settings.misspec.VMA_hor;
    T       = settings.simul.T;
    zeta    = settings.misspec.zeta;

    %----------------------------------------------------------------
    % Parameterization
    %----------------------------------------------------------------
    
    dgp.A  = VAR_p.VAR_coeff;
    dgp.HD = chol(squeeze(VMA.IRF_Wold(1,:,:)) * squeeze(VMA.IRF_Wold(1,:,:))','lower');
    dgp.H  = zeros(n_y,n_y);
    dgp.D  = zeros(n_y,n_y);
    for i_y = 1:n_y
        dgp.H(:,i_y)   = dgp.HD(:,i_y) / dgp.HD(i_y,i_y);
        dgp.D(i_y,i_y) = dgp.HD(i_y,i_y)^2;
    end
    
    dgp.alpha_tilde = zeros(VMA_hor,n_y,n_y);
    for i_l = 1:VMA_hor
        dgp.alpha_tilde(i_l,:,:) = dgp.H^(-1) * squeeze(VMA.IRF_Wold(i_l,:,:)) * squeeze(VMA.IRF_Wold(1,:,:))^(-1) * dgp.H;
    end
    
    %----------------------------------------------------------------
    % Magnitude of Mis-Specification
    %----------------------------------------------------------------
    
    dgp.M_tilde = 0;
    for i_l = 2:VMA_hor
        dgp.M_tilde = dgp.M_tilde + trace(dgp.D * squeeze(dgp.alpha_tilde(i_l,:,:))' * dgp.D^(-1) * squeeze(dgp.alpha_tilde(i_l,:,:)));
    end
    dgp.M_tilde = sqrt(dgp.M_tilde);
    
    clear i_l
    
    M = T^zeta * dgp.M_tilde;
    

end