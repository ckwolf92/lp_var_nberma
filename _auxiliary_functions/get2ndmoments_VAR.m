function y_aux = get2ndmoments_VAR(VAR,settings)

    % Get Second Moments from VAR
    
    % Inputs:
    % VAR       struct          VAR(p) representation
    % settings  struct          estimation settings
    
    % Output:
    % y_aux     struct          implied second-moment properties


    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    VAR_coeff = VAR.VAR_coeff;
    Sigma_u   = VAR.Sigma_u;
    n_y       = size(Sigma_u,1);
    p         = VAR.laglength;
    VMA_hor   = settings.misspec.VMA_hor;
    
    %----------------------------------------------------------------
    % Wold Representation
    %----------------------------------------------------------------
    
    Theta_Wold                = NaN(VMA_hor,n_y,n_y);
    VAR_coeff                 = [VAR_coeff(1:p*n_y,:);zeros(1+VMA_hor*n_y-p*n_y,n_y)];
    
    for i = 1:VMA_hor
        if i == 1
            Theta_Wold(i,:,:) = Sigma_u^(0.5);
        else
            Theta_Wold(i,:,:) = zeros(n_y,n_y);
            for j = 1:i-1
                Theta_Wold(i,:,:) = squeeze(Theta_Wold(i,:,:)) + VAR_coeff(1+(j-1)*n_y:j*n_y,:)' * squeeze(Theta_Wold(i-j,:,:)); % just compute Cov(x_t, eps_{t-h})
            end
        end
    end
    
    %----------------------------------------------------------------
    % Variances and Covariances
    %----------------------------------------------------------------
    
    Cov_y      = cell(VMA_hor,1);
    for i = 1:VMA_hor
        Cov_y{i} = 0;
        for j = 1:(VMA_hor-i+1)
            Cov_y{i} = Cov_y{i} + squeeze(Theta_Wold(j,:,:)) * squeeze(Theta_Wold(j+i-1,:,:))';
        end
    end
    
    Sigma_y_big = NaN(n_y*VMA_hor,n_y*VMA_hor);
    for i = 1:VMA_hor
        for j = 1:VMA_hor
            if i > j
                Sigma_y_big(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_y{1+abs(i-j)};
            else
                Sigma_y_big(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_y{1+abs(i-j)}';
            end
        end
    end
    
    %----------------------------------------------------------------
    % Collect Results
    %----------------------------------------------------------------
    
    y_aux.Theta_Wold        = Theta_Wold;
    y_aux.Cov_y             = Cov_y;
    y_aux.Sigma_u           = Sigma_u;
    y_aux.Sigma_y_big       = Sigma_y_big;

end