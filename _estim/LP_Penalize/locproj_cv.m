function rss = locproj_cv(y, x, w, H_min, H_max, r, lambdaRange, K)
% Function for evaluating MSE of penalized LP using K-fold cross-validation


    % prepare
    w        = [ ones(size(w,1),1) w ]; % Add intercept to controls; to mirror function locproj()
    L        = length(lambdaRange);  
    rss_fold = zeros(L,K); % placeholder for MSE for each lambda and each fold
    
    HR = H_max-H_min+1;
    T = length(y);
    
    % Hold-out sample: design data matrix for penalized LP 
    B  = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
    Xb = nan(T,size(B, 2),HR);
    W  = nan(T,size(w, 2),HR);
    for ih = 1:HR
        h = H_min - 1+ih;
        lagm       = lagmatrix([x,w], h);
        Xb(:,:,ih) = lagm(:,1)*B(ih,:);
        W(:,:,ih)  = lagm(:,2:end); 
    end
    

    % Split sample into K equal-sized chunks
    chunks = ceil(K*(1:T)/T);
    
    for k=1:K % On each fold...

        % Hold-out sample
        the_select = (chunks==k);
        
        % Compute design matrix on estimation sample
        [the_B, the_Xb, the_Y, the_Y_resw, the_Xb_resw] = locproj_design(y(~the_select), x(~the_select), w(~the_select,:), H_min, H_max);


        % Loop over lambda
        for l=1:L
            
            % Estimate for given lambda
            [~, the_theta, the_gamma] = locproj_partitioned(the_Y, the_B, the_Xb, w(~the_select,:), the_Y_resw, the_Xb_resw, r, lambdaRange(l));

            % Compute prediction residuals on hold-out sample
            the_res = nan(sum(the_select),HR);
            for ih=1:HR
                the_res(:, ih) = y(the_select)-Xb(the_select,:,ih)*the_theta-W(the_select,:, ih)*the_gamma(:,ih);
            end
            the_res = the_res(~isnan(the_res)); % Remove NaN
            
            rss_fold(l,k) = (the_res'*the_res)/length(the_res); % MSE
            
        end

    end
    
    % Average MSE across folds
    rss = mean(rss_fold,2);

end