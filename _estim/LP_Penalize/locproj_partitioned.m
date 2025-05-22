function [IR, theta, gamma, IR_se] = locproj_partitioned(Y, B, Xb, w, Y_resw, Xb_resw, r, lambda, nlag)
% Auxiliary function for computing IRF and coefficients in penalized LP, exploiting partitioned formula

    % Uses results from  R. W. Farebrother (1978), Partitioned Ridge Regression, Technometrics, 20:2, 121-122
    % https://doi.org/10.1080/00401706.1978.10489635
    
    %%% Basic %%%
    % T:  sample size
    % K:  number of B-spline basis functions
    % HR: number of horizons
    % p:  number of controls
    
    %%% Input %%%
    % Y:       response variable and leads thereof
    % B:       B-spline
    % Xb:      impulse variable * B-spline
    % w:       controls
    % Y_resw:  Y residualized by w
    % Xb_resw: Xb residualized by w
    % r:       order of finite difference operator
    % lambda:  penalty strength
    % nlag:    Newey-West lag length
    
    %%% Output %%%
    % IR:    impulse response
    % theta: penalized coef, for Xb. Warning: correspond to b_k in our paper
    % gamma: unpenalized coef, for W. Warning: correspond to \zeta_h and \phi_{h,l} in our paper
    % IR_se: standard errors for IR
    
    % prepare
    [T,K,HR] = size(Xb);
    p = size(w,2);
    
    % r-th difference matrix
    D = eye(K);
    for k = 1:r 
        D = diff(D);
    end

    % First compute penalized coefficients
    Xb_resw_stack = reshape(permute(Xb_resw, [1 3 2]), HR*T, K);
    Y_resw_stack = Y_resw(:);
    select = ~isnan(Y_resw_stack);
    Xb_resw_stack_augment = [Xb_resw_stack(select,:); sqrt(lambda)*D];

    theta = Xb_resw_stack_augment\[Y_resw_stack(select); zeros(K-r,1)];

    IR = B * theta; % Impulse responses

    % Then compute unpenalized coefficients horizon by horizon
    gamma = nan(p,HR);
    res = nan(T,HR);
    if nargout>2
        for ih=1:HR
            the_select = ~isnan(Y(:,ih));
            gamma(:,ih) = w(the_select,:)\(Y(the_select,ih)-Xb(the_select,:,ih)*theta);
            res(the_select,ih) = Y(the_select,ih)-Xb(the_select,:,ih)*theta-w(the_select,:)*gamma(:,ih);
        end
    end

    % Standard errors (Barnichon & Brownlees, p. 525)
    if nargout>3

        % "Meat" in sandwich formula
        Xb_resw(isnan(Xb_resw))=0; % Set NaN to zero so they don't contribute to sums
        res(isnan(res))=0;
        meat = zeros(K);
        for l=0:nlag
            Gamma_l = zeros(K);
            for t=l+1:T-HR-1
                Gamma_l = Gamma_l + (permute(Xb_resw(t,:,:),[2 3 1])*res(t,:)')*(res(t-l,:)*permute(Xb_resw(t-l,:,:),[2 3 1])');
            end
            wgt = 1-l/(nlag+1)-0.5*(l==0);
            meat = meat + wgt*(Gamma_l + Gamma_l'); % Newey-West
        end

        % Standard errors
        bread = B/(Xb_resw_stack_augment'*Xb_resw_stack_augment);
        IR_varcov = bread*meat*bread';
        IR_se = sqrt(diag(IR_varcov));

    end
    
end
