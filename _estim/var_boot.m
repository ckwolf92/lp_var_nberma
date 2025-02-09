function Y_boot = var_boot(A, res, Y, p, l, no_const)

    % VAR residual bootstrap, homoskedastic or wild
    
    % Inputs:
    % A         n x np      VAR(p) coefficient matrices [A_1,...,A_p]
    %                       OPTIONAL: column dimension could exceed n*p, in which case last column equals intercept
    %                       NOTE: only A(:,1:n*p) and A(:,end) will be used
    % res       T_res x n   residuals
    % Y         T x n       data vector
    % p         1 x 1       lag length
    % l                     l >=1 block bootstrap (=1 for homosk). 'wild' for wild bootstrap. 
    % no_const  bool        true: exclude intercept from bootstrap samples
    
    % Outputs:
    % Y_boot    T x n       bootstrap sample
    

    % Dimensions
    [T,n] = size(Y);
    T_res = size(res,1); % Effective sample size for residuals (need not equal T-p)
    
    % Draw block of initial T-T_res obs. from real data
    ind_init = randi(T_res+1);
    Y_init = Y(ind_init:ind_init+T-T_res-1,:);
    
    % Intercept (if applicable)
    c = zeros(n,1);
    if ~no_const
        c = A(:,end);
    end

    if ischar(l)
        if strcmp(l, 'wild')  % wild bootstrap
            res_boot = randn(T_res,1).*res;
        end

    else  % Block bootstrap

        % Setup
        n_blocks = ceil(T_res/l);
        res_boot = nan(size(res));
        
        % Draw blocks
        i_blocks            = randi([0, T_res-l], n_blocks, 1);  
        ind_blocks          = reshape(i_blocks + cumsum(ones(n_blocks, l), 2), [], 1);
        ind_blocks          = ind_blocks(1:T_res);
        res_boot_uncentered = res(ind_blocks,:);

        % Recenter, vectorized
        res_center = filter(ones(1, T_res-l+1)/(T_res-l+1), 1, res);
        res_center = res_center(end-l+1:end, :);
        res_center = reshape(permute(repmat(res_center, 1, 1,n_blocks), [1,3,2]), n_blocks*l, []);
        res_center = res_center(1:T_res, :);
        res_boot   = res_boot_uncentered - res_center;

        % Recenter, unvectorized
        % for s = 1:l
        %     for j=0:n_blocks-1
        %         if j*l+s <= T_res
        %             res_boot(j*l+s, :) = res_boot_uncentered(j*l+s, :) - ...
        %                                  mean(res(s+(0:T_res-l), :), 1);
        %         end
        %     end
        % end
        
    end
    
    % Generate VAR(p) data, with residuals and initial conditions as above
    Y_boot = [Y_init; var_sim(A(:,1:n*p), c, res_boot, Y_init)];

end