function VMA = getresidVMA(VAR_infty,VAR_p,settings)

    % Get Residual VMA Process
    
    % Inputs:
    % VAR_infty struct          VAR(\infty) representation
    % VAR_p     struct          VAR(p) representation
    % settings  struct          estimation settings
    
    % Outputs:
    % VMA       struct          residual VMA process


    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    n_y           = size(VAR_infty.Sigma_u,1);
    VMA_hor       = settings.misspec.VMA_hor;
    VAR_laglength = settings.misspec.lags_temp;
    
    %----------------------------------------------------------------
    % VMA Wold Representation
    %----------------------------------------------------------------
    
    IRF_Wold = zeros(n_y,n_y,VMA_hor);
    
    for l = 1:VMA_hor
        for q = 0:min(l-1,VAR_laglength)
            IRF_Wold(:,:,l) = IRF_Wold(:,:,l) + VAR_p.A(:,:,q+1)' * squeeze(VAR_infty.IRF_Wold(l-q,:,:));
        end
    end
    
    VMA.IRF_Wold = permute(IRF_Wold,[3 1 2]);

end