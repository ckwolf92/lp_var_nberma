function VAR_p = popVARp(y_aux,settings)

    % Compute VAR(p) in ABCD system
    
    % Inputs:
    % y_aux     struct          auxiliary properties of observables y
    % settings  struct          estimation settings
    
    % Output:
    % VAR_p     struct          implied VAR(p) model

    
    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    n_y = size(y_aux.Sigma_u,1);
    
    Sigma_y_big = y_aux.Sigma_y_big;
    
    VAR_laglength   = settings.misspec.lags_temp;
    VMA_hor         = settings.misspec.VMA_hor;
    
    VAR_p.laglength = VAR_laglength;
    
    %----------------------------------------------------------------
    % VAR(p)
    %----------------------------------------------------------------
    
    VAR_p.VAR_coeff = Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,n_y+1:n_y+VAR_laglength*n_y)^(-1) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,1:n_y);
    VAR_p.Sigma_u   = Sigma_y_big(1:n_y,1:n_y) ...
                        - Sigma_y_big(1:n_y,n_y+1:n_y+VAR_laglength*n_y) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,n_y+1:n_y+VAR_laglength*n_y)^(-1) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,1:n_y);
    
    VAR_p.A = NaN(n_y,n_y,VAR_laglength+1);
    VAR_p.A(:,:,1) = eye(n_y);
    for l = 2:VAR_laglength+1
        VAR_p.A(:,:,l) = -VAR_p.VAR_coeff(1+(l-2)*n_y:(l-1)*n_y,:);
    end

end