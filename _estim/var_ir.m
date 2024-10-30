function [irs, jacob_a, jacob_nu] = var_ir(A, nu, horzs)

    % SVAR(p) impulse responses and Jacobian wrt. parameters
    
    % Inputs:
    % A         n x np  VAR coefficient matrices (A_1, ..., A_p)
    % nu        n x 1   impact shock vector
    % horzs     1 x H   horizons of interest
    
    % Outputs:
    % irs       n x H               structural impulse responses Theta_h*nu at select horizons
    % jacob_a   n x (n^2*p) x H     Jacobian of Theta_h*nu at select horizons wrt. vec(A)
    % jacob_nu  n x n x H           Jacobian of Theta_h*nu at select horizons wrt. nu
    
    
    % Dimensions
    nh     = length(horzs);
    maxh   = max(horzs);
    [n,np] = size(A);
    p      = np/n;

    % Placeholders
    irs         = zeros(n,nh);
    jacob_a     = zeros(n,np*n,nh);
    jacob_nu    = zeros(n,n,nh);

    % Lags to aid computations
    ir_p = [eye(n); zeros(n*(p-1),n)]; % Will contain last p reduced-form impulse responses, stacked vertically
    jacob_a_p = zeros(n^2,n^2*p,p); % Will contain last p values of the Jacobian of vec(Theta_h) wrt. vec(A)

    % Loop through horizons
    for h=0:maxh
        
        % Compute impulse responses
        if h>0
            the_A       = A(:,1:n*min(h,p));
            the_past_ir = ir_p(1:n*min(h,p),:);
            the_ir      = the_A*the_past_ir; % Impulse response at horizon h
            ir_p        = [the_ir; ir_p(1:end-n,:)]; % Shift forward in time
        else
            the_ir = eye(n);
        end
        
        the_ind = find(horzs==h,1);
        if ~isempty(the_ind)
            irs(:,the_ind) = the_ir*nu; % Store impulse response at select horizons
        end
        
        if nargout>1
            
            % A Jacobian
            the_jacob_a_p = zeros(n^2,n^2*p);

            if h>0
                the_jacob_a_p(:,1:n^2*min(h,p)) = kron(the_past_ir',eye(n));
                for l=1:min(h,p)
                    the_jacob_a_p(:,1:n^2*min(h,p)) = the_jacob_a_p(:,1:n^2*min(h,p)) + kron_fast(A(:,(l-1)*n+1:l*n),jacob_a_p(:,1:n^2*min(h,p),l),1);
                end
                jacob_a_p(:,:,2:end) = jacob_a_p(:,:,1:end-1);
                jacob_a_p(:,:,1)     = the_jacob_a_p;
            end
            
            if ~isempty(the_ind)
                jacob_a(:,:,the_ind) = kron_fast(nu',the_jacob_a_p,0);
                jacob_nu(:,:,the_ind) = the_ir; % nu Jacobian
            end
            
        end
        
    end

end